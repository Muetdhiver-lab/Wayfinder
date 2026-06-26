# -*- coding: utf-8 -*-
"""Regression tests for Wayfinder's known MGA examples.

These tests intentionally avoid running the optimizer. They load completed
reference jobs from the bundled SQLite fixture and verify that the decoding
layer still reproduces the stored results.
"""

import json
import hashlib
import inspect
import math
import shutil
import sys
import tempfile
import unittest
from unittest.mock import patch
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
CORE_DIR = ROOT / "WayfinderCore"
REFERENCE_SOURCE_DB = ROOT / "Tests" / "wayfinder_reference.sqlite"
REFERENCE_TEMP_DIR = tempfile.TemporaryDirectory()
REFERENCE_DB = Path(REFERENCE_TEMP_DIR.name) / REFERENCE_SOURCE_DB.name
shutil.copy2(REFERENCE_SOURCE_DB, REFERENCE_DB)

sys.path.insert(0, str(CORE_DIR))

from _SQLiteStore import SQLiteJobStore  # noqa: E402
from _Wayfinder import Wayfinder, _make_mga_1dsm  # noqa: E402
from _Trajectory import decode_dV_tof as decode_dV_tof_func  # noqa: E402
from _Trajectory import decode_trajectory as decode_trajectory_func  # noqa: E402
from _Trajectory import EJ_angle_from_Pe  # noqa: E402
from _Trajectory import ejection_from_gene  # noqa: E402
from _Trajectory import fast_ejection_from_gene  # noqa: E402
from _Trajectory import _safe_acos  # noqa: E402
from _Optimization import WayfinderFitnessDecorator  # noqa: E402
from _OptimizationService import OptimizationService  # noqa: E402
import pykep as pk  # noqa: E402
import pygmo as pg  # noqa: E402
from pykep.trajopt import mga_1dsm  # noqa: E402


def load_wayfinder_reference(batch_name, planet_pack="Vanilla"):
    """Load a batch of completed reference jobs from the SQLite fixture."""
    plans = Wayfinder(planet_pack=planet_pack, datastore_name=batch_name)
    store = SQLiteJobStore(REFERENCE_DB)
    try:
        rows = store.best_results(
            planet_pack=planet_pack,
            batch_name=batch_name,
            limit=100,
        )
    finally:
        store.close()
    plans._reference_rows = {
        (row["sequence_short_name"], int(row["t0_min"]), int(row["tof_min"])): row
        for row in rows
    }
    return plans


def reference_udp(plans, row):
    bodies = json.loads(row["bodies_json"])
    return _make_mga_1dsm(
        seq=list(map(plans._fullname_dic.get, bodies)),
        t0=[row["t0_min"] / plans._Edy2Kdy, row["t0_max"] / plans._Edy2Kdy],
        tof=[row["tof_min"] / plans._Edy2Kdy, row["tof_max"] / plans._Edy2Kdy],
        vinf=[row["vinf_min"], row["vinf_max"]],
        tof_encoding=row["tof_encoding"],
        add_vinf_dep=bool(row["add_vinf_dep"]),
        add_vinf_arr=bool(row["add_vinf_arr"]),
        multi_objective=bool(row["multi_objective"]),
        orbit_insertion=bool(row["orbit_insertion"]),
        rp_target=row["rp_target"],
        e_target=row["e_target"],
        max_revs=int(row["lambert_max_revs"]),
    )


def decode_job(plans, job_index):
    """Decode a completed job into the main regression quantities."""
    row = plans._reference_rows[job_index]
    udp = reference_udp(plans, row)
    gene = json.loads(row["gene_json"])
    return decode_dV_tof_func(udp, gene, planet_pack=plans.planet_pack)


def decode_trajectory(plans, job_index):
    """Decode a completed job into detailed trajectory metrics."""
    row = plans._reference_rows[job_index]
    udp = reference_udp(plans, row)
    gene = json.loads(row["gene_json"])
    return decode_trajectory_func(udp, gene, planet_pack=plans.planet_pack)


def planet_semimajor_axis(planet):
    if hasattr(planet, "orbital_elements"):
        return planet.orbital_elements[0]
    return planet.elements(0, pk.el_type.KEP_M)[0]


def assert_decoded_matches_stored(
    test_case, plans, job_index, abs_tol=1.0, expected_objective=None,
):
    """Compare decoded values to the reference values stored in SQLite."""
    decoded = decode_job(plans, job_index)
    row = plans._reference_rows[job_index]
    stored = (
        row["objective_dv"] if expected_objective is None else expected_objective,
        row["result_t0"],
        row["result_tof"],
        row["ejection_vinf"],
    )

    for name, decoded_value, stored_value in zip(
        ["result_DV", "result_t0", "result_tof", "result_ej_vinf"],
        decoded,
        stored,
    ):
        test_case.assertAlmostEqual(
            decoded_value,
            stored_value,
            delta=abs_tol,
            msg=f"{job_index} {name}",
        )


class WayfinderRegressionTests(unittest.TestCase):
    def test_release_optimizer_default_is_explicit_three_stage_funnel(self):
        self.assertEqual(Wayfinder.DEFAULT_OPTIMIZER_STRATEGY, "funnel")
        self.assertEqual(
            inspect.signature(Wayfinder.optimize_sqlite)
            .parameters["optimizer_strategy"].default,
            "funnel",
        )

    def test_safe_acos_clamps_normalized_roundoff(self):
        self.assertEqual(_safe_acos(1.0 + 1e-12), 0.0)
        self.assertEqual(_safe_acos(-1.0 - 1e-12), math.pi)

    def test_core_is_pykep3_only(self):
        forbidden_patterns = [
            "pykep2",
            "pk.planet.keplerian",
            "hasattr(pk.planet",
            "fb_prop",
            "get_v1",
            "get_v2",
            "return propagate_lagrangian(r, v, tof, mu)",
            "return ic2par(r, v, mu)",
            "except TypeError as exc",
        ]
        core_files = [
            CORE_DIR / "_Wayfinder.py",
            CORE_DIR / "_Trajectory.py",
            CORE_DIR / "planet_packs" / "vanilla.py",
            CORE_DIR / "planet_packs" / "jnsq.py",
        ]

        offenders = []
        for path in core_files:
            text = path.read_text(encoding="utf-8")
            for pattern in forbidden_patterns:
                if pattern in text:
                    offenders.append(f"{path.name}: {pattern}")

        self.assertEqual(offenders, [])

    def test_wayfinder_does_not_monkey_patch_mga_1dsm(self):
        Wayfinder(planet_pack="Vanilla")

        self.assertFalse(hasattr(mga_1dsm, "decode_trajectory"))
        self.assertFalse(hasattr(mga_1dsm, "decode_dV_tof"))
        self.assertFalse(hasattr(mga_1dsm, "transx"))

    def test_wayfinder_accepts_planet_pack_as_first_argument(self):
        plans = Wayfinder("JNSQ")

        self.assertEqual(plans.planet_pack, "JNSQ")
        self.assertEqual(plans._Edy2Kdy, 2)

    def test_optimizer_island_count_can_use_cpu_count_minus_reserved_cores(self):
        with patch("os.cpu_count", return_value=8):
            self.assertEqual(Wayfinder._resolve_island_count(16), 6)
            self.assertEqual(Wayfinder._resolve_island_count(4), 4)

        with patch("os.cpu_count", return_value=16):
            self.assertEqual(Wayfinder._resolve_island_count(24), 14)

    def test_optimizer_island_count_can_keep_preset_value(self):
        with patch("os.cpu_count", return_value=8):
            self.assertEqual(Wayfinder._resolve_island_count(16, auto_workers=False), 16)

    def test_optimizer_facade_delegates_worker_policy_to_service(self):
        with patch.object(
            OptimizationService,
            "resolve_island_count",
            return_value=3,
        ) as resolve:
            result = Wayfinder._resolve_island_count(
                12, auto_workers=True, reserve_cores=4,
            )

        self.assertEqual(result, 3)
        resolve.assert_called_once_with(
            12, auto_workers=True, reserve_cores=4,
        )

    def test_archipelago_topology_factory_supports_benchmark_topologies(self):
        self.assertIsInstance(Wayfinder._make_archipelago_topology("unconnected", 4), pg.unconnected)
        ring = Wayfinder._make_archipelago_topology("ring", 4)
        full = Wayfinder._make_archipelago_topology("fully_connected", 4)
        self.assertIsInstance(ring, pg.ring)
        self.assertIsInstance(full, pg.fully_connected)
        self.assertEqual(pg.topology(ring).to_networkx().number_of_nodes(), 0)
        self.assertEqual(pg.topology(full).to_networkx().number_of_nodes(), 0)

        for name in ("ring", "fully_connected"):
            archi = pg.archipelago(
                algo=pg.de(gen=1), prob=pg.rosenbrock(2), n=4, pop_size=8,
                t=Wayfinder._make_archipelago_topology(name, 4),
            )
            self.assertEqual(len(archi), 4)
            self.assertEqual(
                archi.get_topology().to_networkx().number_of_nodes(), 4
            )
            before = [island.get_population() for island in archi]
            origins = Wayfinder._population_id_origins(before)
            self.assertTrue(all(isinstance(value, set) for value in origins.values()))
            archi.evolve(1)
            archi.wait_check()
            populations = [island.get_population() for island in archi]
            telemetry = Wayfinder._archipelago_telemetry(
                archi, populations, origins, name
            )
            self.assertEqual(telemetry["island_count"], 4)
            self.assertEqual(telemetry["topology_vertices"], 4)
            self.assertGreater(telemetry["topology_edges"], 0)
            self.assertEqual(telemetry["migrants_published"], 4)

        with self.assertRaises(ValueError):
            Wayfinder._make_archipelago_topology("unknown", 4)
        with self.assertRaises(ValueError):
            Wayfinder._make_archipelago_topology("ring", 0)

    def test_population_id_origins_retain_all_source_islands(self):
        class PopulationIds:
            def __init__(self, ids):
                self._ids = ids

            def get_ID(self):
                return self._ids

        origins = Wayfinder._population_id_origins(
            [PopulationIds([1, 2]), PopulationIds([1, 3])]
        )

        self.assertEqual(origins[1], {0, 1})
        self.assertEqual(origins[2], {0})
        self.assertEqual(origins[3], {1})

    def test_island_migration_policy_controls_published_individual_count(self):
        archi = pg.archipelago(t=pg.ring())
        for seed in range(4):
            archi.push_back(pg.island(
                algo=pg.de(gen=1, seed=seed),
                prob=pg.rosenbrock(2),
                size=8,
                seed=seed,
                s_pol=pg.select_best(2),
                r_pol=pg.fair_replace(2),
            ))

        archi.evolve(1)
        archi.wait_check()

        self.assertEqual(
            sum(len(entry[0]) for entry in archi.get_migrants_db()), 8
        )

    def test_funnel_stage_plan_strictly_reduces_small_archipelago(self):
        stages = Wayfinder._funnel_stage_plan(8, 32, 20, 20)

        self.assertEqual([stage["n_island"] for stage in stages], [8, 4, 2])
        self.assertEqual([stage["island_pop"] for stage in stages], [32, 14, 14])
        self.assertEqual([stage["evo_steps"] for stage in stages], [20, 40, 5])
        self.assertEqual(stages[-1]["ejection_model"], "vector_3d")
        self.assertEqual(stages[-1]["sade_gen"], 500)
        self.assertEqual(stages[-1]["initialization"], "random_plus_champion")
        self.assertEqual(stages[-1]["algorithm"], "hybrid")

    def test_optimization_service_preserves_all_funnel_plan_variants(self):
        expected_digests = {
            "legacy": "215c26160222c6a49320969561ee9288dc5d73ddcf44cca7080adc426e6a917c",
            "local": "db7f37234655ae43d6aed232ef0eb8d4ba167becd9655113e74b6f725e9f8f8a",
            "hybrid": "efe3a5db0cb574cc653554ce01bbd40776072082129151c137b270fe4711c4f4",
            "phase_elites_nm": "6497b3c946b5546601231c168028e9e9b38be7f59bc62626377c256fba5013b8",
            "phase_elites_nm_equal": "3ae7b2c97573b7110b6960eb1fe75a738159aea569f68be965912d256bae0236",
            "scout_archive_nm_32": "8c2276fdb73c63c4cedf9647c617949b650c5e9b8e7387817e9791aa219c45b8",
            "scout_archive_nm": "528ae24780546649e20ae5f4e53142422bec631b33a9b3ff443def3283bf5c3a",
            "scout_archive_nm_64": "528ae24780546649e20ae5f4e53142422bec631b33a9b3ff443def3283bf5c3a",
            "scout_archive_nm_128": "8d321031598d7dc5b0e33bd8f26d6b6053b4294b28feaa2798fc47425595e582",
        }

        for strategy, expected_digest in expected_digests.items():
            with self.subTest(strategy=strategy):
                plan = OptimizationService.funnel_stage_plan(
                    16, 32, 20, 20, exact_strategy=strategy,
                )
                payload = json.dumps(
                    plan, sort_keys=True, separators=(",", ":"),
                ).encode()
                self.assertEqual(
                    hashlib.sha256(payload).hexdigest(), expected_digest,
                )

    def test_funnel_stage_plan_enforces_sade_minimum_population(self):
        stages = Wayfinder._funnel_stage_plan(1, 1, 1, 1)

        self.assertEqual([stage["island_pop"] for stage in stages], [7, 7, 7])
        population = pg.population(pg.rosenbrock(3), size=stages[0]["island_pop"])
        pg.algorithm(pg.sade(gen=1)).evolve(population)

    def test_annealing_schedule_matches_sade_evaluation_budget(self):
        problem = pg.problem(pg.rosenbrock(18))
        schedule = Wayfinder._balanced_annealing_schedule(
            problem, population_size=14, sade_gen=20
        )

        self.assertEqual(schedule["target_evaluations"], 280)
        self.assertLessEqual(
            abs(schedule["nominal_evaluations"] - 280), problem.get_nx() * 2
        )

    def test_fast_exact_funnel_policies_are_explicit(self):
        local = Wayfinder._funnel_stage_plan(8, 32, 20, 20, exact_strategy="local")[-1]
        hybrid = Wayfinder._funnel_stage_plan(8, 32, 20, 20, exact_strategy="hybrid")[-1]

        self.assertEqual(local["initialization"], "local")
        self.assertEqual(hybrid["initialization"], "local_global")
        self.assertEqual(local["sade_gen"], 50)
        self.assertEqual(local["evo_steps"], 10)
        self.assertEqual(local["adaptive_stop"]["min_steps"], 2)

    def test_phase_elite_funnel_keeps_eight_islands_for_exact_portfolio(self):
        stages = Wayfinder._funnel_stage_plan(
            16, 32, 5, 20, exact_strategy="phase_elites_nm"
        )

        self.assertEqual([stage["n_island"] for stage in stages], [16, 8, 8])
        self.assertEqual(stages[1]["initialization"], "phase_elites_mixed")
        self.assertEqual(stages[1]["elite_fraction"], 0.35)
        self.assertEqual(stages[2]["algorithm"], "sade_nlopt")
        self.assertEqual(stages[2]["sade_gen"], 100)

    def test_equal_phase_elite_funnel_uses_same_stage_one_and_two_epochs(self):
        stages = Wayfinder._funnel_stage_plan(
            16, 32, 25, 50, exact_strategy="phase_elites_nm_equal"
        )

        self.assertEqual([stage["evo_steps"] for stage in stages], [25, 25, 5])

    def test_scout_archive_funnel_is_wide_budget_neutral_and_adaptive(self):
        stages = Wayfinder._funnel_stage_plan(
            16, 32, 20, 50, exact_strategy="scout_archive_nm"
        )

        self.assertEqual(
            [stage["name"] for stage in stages],
            ["scout_unconnected", "wide", "intermediate", "exact_ejection"],
        )
        self.assertEqual([stage["n_island"] for stage in stages], [64, 16, 8, 8])
        self.assertEqual([stage["island_pop"] for stage in stages], [8, 32, 14, 14])
        self.assertEqual([stage["evo_steps"] for stage in stages], [5, 15, 40, 10])
        self.assertEqual(stages[0]["topology"], "unconnected")
        self.assertEqual(stages[0]["migration_rate"], 0)
        self.assertEqual(stages[1]["initialization"], "scout_diverse")
        self.assertTrue(stages[1]["archive_exact"])
        self.assertTrue(stages[2]["archive_exact"])
        self.assertTrue(stages[3]["use_exact_archive"])
        self.assertEqual(stages[3]["adaptive_stop"]["min_steps"], 5)
        self.assertFalse(
            stages[3]["adaptive_stop"]["require_average_plateau"]
        )

        old_wide_evaluations = 16 * 32 * 20 * 50
        scout_plus_wide = (
            64 * 8 * 5 * 50 + 16 * 32 * 15 * 50
        )
        self.assertEqual(scout_plus_wide, old_wide_evaluations)

    def test_scout_width_variants_hold_level_zero_budget_constant(self):
        expected = {
            "scout_archive_nm_32": (32, 100),
            "scout_archive_nm_64": (64, 50),
            "scout_archive_nm_128": (128, 25),
        }
        budgets = []
        for strategy, (islands, generations) in expected.items():
            stages = Wayfinder._funnel_stage_plan(
                16, 32, 20, 50, exact_strategy=strategy
            )
            scout = stages[0]
            self.assertEqual(scout["n_island"], islands)
            self.assertEqual(scout["sade_gen"], generations)
            self.assertEqual(stages[1]["evo_steps"], 15)
            budgets.append(
                islands * scout["island_pop"]
                * scout["evo_steps"] * generations
            )
        self.assertEqual(budgets, [128000, 128000, 128000])

    def test_exact_phase_archive_is_bounded_and_keeps_exact_best(self):
        plans = Wayfinder(planet_pack="Vanilla")
        problem = pg.problem(pg.rosenbrock(2))
        genes = [[1.0, 1.0], [0.0, 0.0], [2.0, 2.0], [-1.0, 1.0]]

        with patch.object(
            plans, "_encounter_phase_embedding",
            side_effect=lambda context, problem, gene: np.asarray(gene),
        ):
            archive = plans._update_exact_phase_archive(
                {}, problem, [], genes, max_size=3
            )

        self.assertEqual(len(archive), 3)
        self.assertIn([1.0, 1.0], archive)

    def test_exact_diverse_seed_selection_is_preserved_by_service(self):
        problem = pg.problem(pg.rosenbrock(2))
        genes = [
            [1.0, 1.0], [0.0, 0.0], [2.0, 2.0],
            [-1.0, 1.0], [3.0, -2.0],
        ]

        selected = OptimizationService.select_exact_diverse_seeds(
            problem, genes, 3,
        )

        self.assertEqual(
            selected,
            [[1.0, 1.0], [3.0, -2.0], [-1.0, 1.0]],
        )
        self.assertEqual(
            Wayfinder._select_exact_diverse_seeds(problem, genes, 3),
            selected,
        )

    def test_exact_diverse_seed_selection_respects_zero_count(self):
        problem = pg.problem(pg.rosenbrock(2))
        genes = [[1.0, 1.0], [0.0, 0.0]]

        self.assertEqual(
            OptimizationService.select_exact_diverse_seeds(problem, genes, 0),
            [],
        )
        self.assertEqual(
            OptimizationService.select_exact_diverse_seeds(problem, genes, -1),
            [],
        )
        self.assertEqual(
            OptimizationService.select_exact_diverse_seeds(problem, genes, 1),
            [[1.0, 1.0]],
        )
        self.assertEqual(
            OptimizationService.select_exact_diverse_seeds(problem, [], 1),
            [],
        )

    def test_phase_embedding_receives_planet_context_explicitly(self):
        class LinearBody:
            def __init__(self, offset):
                self.offset = float(offset)

            def eph(self, epoch):
                value = float(epoch) + self.offset
                return [value, 1.0, 0.0], [0.0, value, 1.0]

        bodies = {"A": LinearBody(1.0), "B": LinearBody(2.0)}
        context = {"bodies_json": '["A", "B"]', "tof_encoding": "direct"}
        problem = pg.problem(pg.rosenbrock(6))
        gene = [10.0, 0.0, 0.0, 0.0, 0.0, 5.0]

        direct = OptimizationService.encounter_phase_embedding(
            bodies, context, problem, gene,
        )
        plans = Wayfinder(planet_pack="Vanilla")
        plans._fullname_dic = bodies

        np.testing.assert_allclose(
            plans._encounter_phase_embedding(context, problem, gene), direct,
        )
        self.assertEqual(direct.shape, (13,))

    def test_funnel_stage_plan_reduces_wide_archipelago_to_at_most_eight(self):
        stages = Wayfinder._funnel_stage_plan(48, 25, 10, 100)

        self.assertEqual([stage["n_island"] for stage in stages], [48, 24, 8])
        self.assertTrue(
            all(
                later["n_island"] < earlier["n_island"]
                for earlier, later in zip(stages, stages[1:])
            )
        )

    def test_adaptive_stop_continues_while_population_is_improving(self):
        history = [
            {"best": best, "average": average}
            for best, average in (
                (9000, 250000), (8500, 130000), (7500, 80000),
                (6000, 52000), (5200, 36000),
            )
        ]
        should_stop, plateau_windows = Wayfinder._adaptive_stop_decision(
            history, 0, {"min_steps": 5, "window": 3, "patience": 2}
        )
        self.assertFalse(should_stop)
        self.assertEqual(plateau_windows, 0)

    def test_adaptive_stop_can_ignore_population_average_for_exact_refinement(self):
        history = [
            {"best": best, "average": average}
            for best, average in (
                (2000.0, 5000.0), (1999.5, 4500.0), (1999.2, 4000.0),
                (1999.0, 3500.0), (1998.9, 3000.0), (1998.8, 2500.0),
            )
        ]
        options = {
            "min_steps": 5,
            "window": 2,
            "patience": 2,
            "best_relative_tolerance": 0.001,
            "average_relative_tolerance": 0.005,
            "require_average_plateau": False,
        }

        first_stop, windows = Wayfinder._adaptive_stop_decision(
            history[:-1], 0, options
        )
        second_stop, windows = Wayfinder._adaptive_stop_decision(
            history, windows, options
        )

        self.assertFalse(first_stop)
        self.assertTrue(second_stop)

    def test_adaptive_stop_requires_two_plateau_windows(self):
        options = {"min_steps": 5, "window": 3, "patience": 2}
        history = [{"best": 5000, "average": 10000} for _ in range(5)]
        should_stop, plateau_windows = Wayfinder._adaptive_stop_decision(
            history, 0, options
        )
        self.assertFalse(should_stop)
        history.append({"best": 5000, "average": 10000})
        should_stop, plateau_windows = Wayfinder._adaptive_stop_decision(
            history, plateau_windows, options
        )
        self.assertTrue(should_stop)

    def test_alpha_gene_to_direct_gene_preserves_keemo_trajectory(self):
        plans = Wayfinder(planet_pack="JNSQ")
        store = SQLiteJobStore(REFERENCE_DB)
        try:
            best = store.best_results(
                planet_pack="JNSQ",
                batch_name="JNSQ_KEEMo_5yr_high_vinf_pykep3",
                sequence_short_name="KEEMo",
                limit=1,
            )[0]
            context = store.run_context(best["run_id"])
        finally:
            store.close()

        alpha_gene = json.loads(context["gene_json"])
        alpha_udp = plans._mga_problem_from_sqlite_context(context)
        direct_udp, direct_gene, decoded_leg_tofs = plans._direct_local_problem_from_sqlite_context(
            context,
            alpha_gene,
            t0_wiggle_days=5,
            leg_tof_wiggle_days=5,
        )

        alpha_decoded = decode_dV_tof_func(alpha_udp, alpha_gene, planet_pack="JNSQ")
        direct_decoded = decode_dV_tof_func(direct_udp, direct_gene, planet_pack="JNSQ")

        self.assertEqual(len(decoded_leg_tofs), 3)
        self.assertEqual(len(direct_gene), len(alpha_gene) - 1)
        for alpha_value, direct_value in zip(alpha_decoded, direct_decoded):
            self.assertAlmostEqual(alpha_value, direct_value, places=6)

    def test_optional_sequence_generation_for_keemo(self):
        plans = Wayfinder(planet_pack="Vanilla")

        short_sequences = plans.generateShortSequences(
            [["Kerbin"], ["Eve"], ["Eve", "*"], ["Moho"]]
        )

        self.assertEqual(short_sequences, ["KEEMo", "KEMo"])

    def test_auto_tof_for_known_vanilla_sequences(self):
        plans = Wayfinder(planet_pack="Vanilla")

        self.assertEqual(plans.auto_tof(["Kerbin", "Eve", "Eve", "Moho"]), (400, 800))
        self.assertEqual(
            plans.auto_tof(["Kerbin", "Eve", "Kerbin", "Kerbin", "Jool"]),
            (1700, 3500),
        )

    def test_planner_style_direct_tof_bounds_exclude_known_kekkj_solution(self):
        plans = Wayfinder(planet_pack="Vanilla")
        estimate = plans.estimate_direct_tof_bounds(
            ["Kerbin", "Eve", "Kerbin", "Kerbin", "Jool"],
            profile="planner",
        )

        self.assertLess(estimate["total_upper_days"], 3469.52)
        self.assertAlmostEqual(estimate["total_upper_days"], 3047.279, places=3)

    def test_relaxed_direct_tof_bounds_include_known_kekkj_legs(self):
        plans = Wayfinder(planet_pack="Vanilla")
        estimate = plans.estimate_direct_tof_bounds(
            ["Kerbin", "Eve", "Kerbin", "Kerbin", "Jool"],
            profile="relaxed",
        )
        known_leg_tofs = [205.0734, 402.1447, 852.1840, 2010.1171]

        for known, (lower, upper) in zip(
            known_leg_tofs, estimate["direct_bounds_days"]
        ):
            self.assertLessEqual(lower, known)
            self.assertGreaterEqual(upper, known)
        self.assertAlmostEqual(estimate["total_lower_days"], 1833.863, places=3)
        self.assertAlmostEqual(estimate["total_upper_days"], 3991.993, places=3)

    def test_vanilla_keemo_reference_decodes_from_sqlite(self):
        plans = load_wayfinder_reference("Ex1_KEEMo_MGA", planet_pack="Vanilla")

        self.assertEqual(list(plans._reference_rows), [("KEEMo", 0, 400)])
        assert_decoded_matches_stored(
            self, plans, ("KEEMo", 0, 400),
            expected_objective=2628.874458141246,
        )

    def test_legacy_decode_tuple_matches_structured_decode(self):
        plans = load_wayfinder_reference("Ex1_KEEMo_MGA", planet_pack="Vanilla")
        job_index = ("KEEMo", 0, 400)

        legacy = decode_job(plans, job_index)
        decoded = decode_trajectory(plans, job_index)

        self.assertEqual(
            legacy,
            (
                decoded["objective_dv"],
                decoded["t0"],
                decoded["tof"],
                decoded["ejection_vinf"],
            ),
        )

    def test_circular_capture_metrics_are_available(self):
        plans = load_wayfinder_reference("Ex1_KEEMo_MGA", planet_pack="Vanilla")
        decoded = decode_trajectory(plans, ("KEEMo", 0, 400))

        self.assertIsNotNone(decoded["capture_dv"])
        self.assertAlmostEqual(decoded["objective_dv"], decoded["dv_with_capture"])
        self.assertGreater(decoded["dv_with_arrival_vinf"], decoded["dv_without_arrival"])

    def test_ejection_angle_from_periapsis_includes_soi_zenith_angle(self):
        plans = load_wayfinder_reference("Ex1_KEEMo_MGA", planet_pack="Vanilla")
        job_index = ("KEEMo", 0, 400)
        row = plans._reference_rows[job_index]
        udp = reference_udp(plans, row)
        decoded = decode_trajectory(plans, job_index)

        rsoi = planet_semimajor_axis(udp._seq[0]) * (
            udp._seq[0].mu_self / udp._seq[0].mu_central_body
        ) ** 0.4
        alt = row["ejection_altitude"]
        ejection_vinf = decoded["ejection_vinf"]
        orbital_speed = math.sqrt(udp._seq[0].mu_self / (udp._seq[0].radius + alt))
        v1 = math.sqrt(
            ejection_vinf * ejection_vinf
            + 2 * orbital_speed * orbital_speed
            - 2 * udp._seq[0].mu_self / rsoi
        )
        eccentricity = (udp._seq[0].radius + alt) * v1 * v1 / udp._seq[0].mu_self - 1
        semi_major_axis = (udp._seq[0].radius + alt) / (1 - eccentricity)
        true_anomaly_at_soi = math.acos(
            (semi_major_axis * (1 - eccentricity * eccentricity) - rsoi)
            / (eccentricity * rsoi)
        )
        zenith_angle_at_soi = math.asin(
            v1 * (udp._seq[0].radius + alt) / (ejection_vinf * rsoi)
        )

        self.assertAlmostEqual(
            EJ_angle_from_Pe(udp, ejection_vinf, rsoi, alt),
            true_anomaly_at_soi + zenith_angle_at_soi,
        )
        self.assertGreater(zenith_angle_at_soi, 0)

    def test_exact_fitness_and_decoder_share_ejection_calculation(self):
        plans = load_wayfinder_reference("Ex1_KEEMo_MGA", planet_pack="Vanilla")
        job_index = ("KEEMo", 0, 400)
        row = plans._reference_rows[job_index]
        udp = reference_udp(plans, row)
        gene = json.loads(row["gene_json"])
        decoded = decode_trajectory_func(udp, gene, planet_pack="Vanilla")
        rsoi = plans.soi_radius(udp._seq[0])

        ejection = ejection_from_gene(
            udp, gene, rsoi, row["ejection_altitude"]
        )
        self.assertAlmostEqual(ejection["dv"], decoded["ejection_dv"], places=9)

        decorator = WayfinderFitnessDecorator(
            planet_pack="Vanilla",
            bodies_by_name=plans._fullname_dic,
            soi_radius_by_name={
                name: plans.soi_radius(body)
                for name, body in plans._fullname_dic.items()
            },
            ejection_altitude=row["ejection_altitude"],
            tof_encoding=row["tof_encoding"],
            ejection_model="vector_3d",
        )
        problem = pg.problem(
            pg.decorator_problem(udp, fitness_decorator=decorator)
        )
        self.assertAlmostEqual(
            problem.fitness(gene)[0], decoded["objective_dv"], places=6
        )

    def test_fast_ejection_selects_minimum_of_direct_and_soi_split(self):
        plans = load_wayfinder_reference("Ex1_KEEMo_MGA", planet_pack="Vanilla")
        row = plans._reference_rows[("KEEMo", 0, 400)]
        udp = reference_udp(plans, row)
        gene = json.loads(row["gene_json"])
        gene[2] = 0.75
        gene[3] = 1800.0

        ejection = fast_ejection_from_gene(
            udp._seq[0], gene, row["ejection_altitude"]
        )
        self.assertAlmostEqual(
            ejection["split_dv"],
            ejection["planar_dv"] + ejection["soi_correction_dv"],
        )
        self.assertEqual(
            ejection["dv"], min(ejection["direct_dv"], ejection["split_dv"])
        )
        self.assertEqual(
            ejection["strategy"],
            "direct" if ejection["direct_dv"] <= ejection["split_dv"]
            else "split_soi",
        )

    def test_vanilla_kekkj_reference_decodes_from_sqlite(self):
        plans = load_wayfinder_reference("Ex2_KEKKJ_MGA", planet_pack="Vanilla")

        self.assertEqual(list(plans._reference_rows), [("KEKKJ", 600, 3000)])
        assert_decoded_matches_stored(
            self, plans, ("KEKKJ", 600, 3000),
            expected_objective=1801.401931428046,
        )

    def test_jnsq_moho_reference_decodes_with_jnsq_time_scale(self):
        plans = load_wayfinder_reference("JNSQ_Moho_reference_fixed", planet_pack="JNSQ")

        self.assertEqual(len(plans._reference_rows), 8)

        expected_objectives = {
            ("KEMo", 200, 450): 5865.630419341301,
            ("KEEMo", 200, 450): 5924.916955843217,
            ("KEMo", 200, 300): 6577.819321714192,
            ("KEEMo", 0, 300): 6749.088718787056,
            ("KEMo", 0, 450): 6919.519262156915,
            ("KEMo", 0, 300): 7118.4947561398185,
            ("KEEMo", 0, 450): 7357.605558481762,
            ("KEEMo", 200, 300): 7687.816056325191,
        }
        for job_index in plans._reference_rows:
            assert_decoded_matches_stored(
                self, plans, job_index,
                expected_objective=expected_objectives[job_index],
            )

    def test_vinf_metrics_are_split_for_jnsq_reference(self):
        plans = load_wayfinder_reference("JNSQ_Moho_reference_fixed", planet_pack="JNSQ")
        decoded = decode_trajectory(plans, ("KEEMo", 0, 300))

        self.assertIsNone(decoded["capture_dv"])
        self.assertIsNone(decoded["dv_with_capture"])
        self.assertAlmostEqual(
            decoded["objective_dv"],
            decoded["dv_without_arrival"] + decoded["arrival_vinf"],
        )
        self.assertAlmostEqual(decoded["objective_dv"], decoded["dv_with_arrival_vinf"])

    def test_direct_batch_bins_only_t0_and_persists_leg_bounds(self):
        plans = Wayfinder(planet_pack="Vanilla", datastore_name="direct-test")
        template = [["Kerbin"], ["Eve"], ["Kerbin"], ["Kerbin"], ["Jool"]]
        with tempfile.TemporaryDirectory() as directory:
            db_path = Path(directory) / "direct.sqlite"
            plans.add_direct_t0_batch_sqlite(
                template, db_path, batch_name="KEKKJ-direct",
                t0_min=500, t0_bin=50, n_t0_bins=3,
                tof_profile="relaxed", arrival_mode="flyby",
            )
            store = SQLiteJobStore(db_path)
            try:
                jobs = store.pending_jobs("Vanilla", limit=10, batch_name="KEKKJ-direct")
            finally:
                store.close()

        self.assertEqual(len(jobs), 3)
        self.assertEqual([job["t0_min"] for job in jobs], [500.0, 550.0, 600.0])
        self.assertTrue(all(job["tof_encoding"] == "direct" for job in jobs))
        self.assertTrue(all(job["arrival_mode"] == "flyby" for job in jobs))
        self.assertTrue(all(not job["add_vinf_arr"] for job in jobs))
        self.assertTrue(all(not job["orbit_insertion"] for job in jobs))
        self.assertEqual(len(json.loads(jobs[0]["leg_tof_bounds_json"])), 4)

        udp = plans._mga_problem_from_sqlite_context(jobs[0])
        self.assertEqual(udp._tof_encoding, "direct")


if __name__ == "__main__":
    unittest.main(verbosity=2)
