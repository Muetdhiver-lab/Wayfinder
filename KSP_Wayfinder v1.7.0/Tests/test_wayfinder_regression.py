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
TESTS_DIR = ROOT / "Tests"
REFERENCE_SOURCE_DB = ROOT / "Tests" / "wayfinder_reference.sqlite"
REFERENCE_TEMP_DIR = tempfile.TemporaryDirectory()
REFERENCE_DB = Path(REFERENCE_TEMP_DIR.name) / REFERENCE_SOURCE_DB.name
shutil.copy2(REFERENCE_SOURCE_DB, REFERENCE_DB)

sys.path.insert(0, str(CORE_DIR))
sys.path.insert(0, str(TESTS_DIR))

from _SQLiteStore import SQLiteJobStore  # noqa: E402
from _Wayfinder import Wayfinder, _make_mga_1dsm  # noqa: E402
from _Trajectory import decode_dV_tof as decode_dV_tof_func  # noqa: E402
from _Trajectory import decode_trajectory as decode_trajectory_func  # noqa: E402
from _Trajectory import EJ_angle_from_Pe  # noqa: E402
from _Trajectory import EJ_details  # noqa: E402
from _Trajectory import ejection_from_gene  # noqa: E402
from _Trajectory import fast_ejection_from_gene  # noqa: E402
from _Trajectory import fast_ejection_from_vinf  # noqa: E402
from _Trajectory import Kdate  # noqa: E402
from _Trajectory import _lambert_v1  # noqa: E402
from _Trajectory import _safe_acos  # noqa: E402
from _KSPTranslator import translate_first_leg  # noqa: E402
from _KSPTranslator import translate_mga_trajectory  # noqa: E402
from _KSPBridge import build_first_leg_flight_plan  # noqa: E402
from _KSPBridge import build_ksp_flight_plan  # noqa: E402
from pykep.core import DAY2SEC, lambert_problem  # noqa: E402
from _Optimization import WayfinderFitnessDecorator  # noqa: E402
from _OptimizationService import OptimizationService  # noqa: E402
from _OptimizationService import StageConfig, FunnelConfig  # noqa: E402
from _SequenceScout import (  # noqa: E402
    SequenceScout,
    TisserandScoutConfig,
    maximum_flyby_turn_angle,
    required_flyby_periapsis,
    tisserand_parameter,
)
from _LambertArcFilter import LambertArc1Config, LambertArc1Filter  # noqa: E402
from benchmark_analysis import (  # noqa: E402
    case_from_sequence,
    grouped_summary,
    summarize_values,
)
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


def ksp_plan_for_reference(plans, job_index, strategy="split_soi"):
    row = plans._reference_rows[job_index]
    udp = reference_udp(plans, row)
    gene = json.loads(row["gene_json"])
    arrival_mode = row.get("arrival_mode") or row["opt_injection"]
    if arrival_mode == "legacy":
        arrival_mode = row["opt_injection"]
    return build_ksp_flight_plan(
        udp,
        gene,
        planet_pack=plans.planet_pack,
        parking_altitude=row["ejection_altitude"],
        strategy=strategy,
        arrival_mode=arrival_mode,
    )


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

    def test_initial_seed_genes_is_forwarded_to_internal_optimizer(self):
        class FakeStore:
            def __init__(self):
                self.claims = 0

            def claim_pending_jobs(self, *args, **kwargs):
                if self.claims:
                    return []
                self.claims += 1
                return [{
                    "job_id": 1,
                    "worker_id": "test-worker",
                    "claimed_at": "2026-01-01T00:00:00",
                    "sequence_short_name": "KEMo",
                    "t0_min": 0.0,
                    "tof_min": 0.0,
                }]

            def close(self):
                pass

        seed_genes = [[1.0, 2.0, 3.0]]
        plans = Wayfinder(planet_pack="Vanilla")
        store = FakeStore()

        with patch("_SQLiteStore.SQLiteJobStore", return_value=store), \
             patch.object(Wayfinder, "_optimize_sqlite_job") as optimize_job:
            count = plans.optimize_sqlite(
                "unused.sqlite",
                n=1,
                batch_name="seeded",
                auto_workers=False,
                initial_seed_genes=seed_genes,
            )

        self.assertEqual(count, 1)
        self.assertIs(
            optimize_job.call_args.kwargs["initial_seed_genes"],
            seed_genes,
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

    def test_split_ring_topology_keeps_two_rings_and_one_weak_bridge(self):
        topology = Wayfinder._make_split_ring_topology(
            current_islands=16,
            alternative_islands=4,
            bridge_weight=0.25,
        )
        graph = pg.topology(topology).to_networkx()
        self.assertEqual(graph.number_of_nodes(), 20)
        self.assertEqual(graph.number_of_edges(), 42)
        self.assertEqual(graph[15][16]["weight"], 0.25)
        self.assertEqual(graph[16][15]["weight"], 0.25)
        self.assertFalse(graph.has_edge(0, 16))
        self.assertTrue(graph.has_edge(0, 1))
        self.assertTrue(graph.has_edge(16, 19))

        with self.assertRaisesRegex(ValueError, "at least three"):
            Wayfinder._make_split_ring_topology(16, 2, 0.25)

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

    def test_funnel_config_objects_round_trip_to_canonical_dicts(self):
        stages = OptimizationService.funnel_stage_plan(
            16, 32, 20, 20, exact_strategy="scout_archive_nm_32",
        )
        stage_round_trip = [StageConfig.from_dict(stage).to_dict() for stage in stages]
        self.assertEqual(stage_round_trip, stages)

        requested = {
            "n_islands": 16,
            "island_pop": 32,
            "evo_steps": 20,
            "sade_gen": 20,
        }
        funnel = FunnelConfig.from_stage_dicts(
            "funnel_scout_archive_32",
            "scout_archive_nm_32",
            requested,
            stages,
        )
        expected = {
            "kind": "funnel",
            "optimizer_strategy": "funnel_scout_archive_32",
            "exact_strategy": "scout_archive_nm_32",
            "requested": requested,
            "stages": stages,
        }

        self.assertEqual(funnel.to_dict(), expected)
        self.assertEqual(
            OptimizationService.funnel_run_config(
                "funnel_scout_archive_32", 16, 32, 20, 20,
            ),
            expected,
        )
        self.assertEqual(
            hashlib.sha256(
                json.dumps(
                    funnel.to_dict(), sort_keys=True, separators=(",", ":"),
                ).encode()
            ).hexdigest(),
            hashlib.sha256(
                json.dumps(expected, sort_keys=True, separators=(",", ":")).encode()
            ).hexdigest(),
        )

    def test_stage_config_rejects_invalid_runner_values(self):
        valid = {
            "name": "wide",
            "n_island": 4,
            "island_pop": 7,
            "evo_steps": 1,
            "sade_gen": 1,
            "ejection_model": "approximate",
            "initialization": "random",
            "algorithm": "sade",
            "topology": "ring",
            "migration_rate": 0,
        }

        for field, value, error in (
            ("n_island", 0, "n_island must be positive"),
            ("island_pop", 6, "island_pop must be at least 7"),
            ("evo_steps", 0, "evo_steps must be positive"),
            ("sade_gen", 0, "sade_gen must be positive"),
            ("algorithm", "scipy_nm", "algorithm must be one of"),
            ("ejection_model", "exactish", "ejection_model must be one of"),
            ("topology", "mesh", "topology must be one of"),
            ("migration_rate", -1, "migration_rate must be non-negative"),
            ("adaptive_stop", True, "adaptive_stop must be a dict"),
            ("annealing", True, "annealing must be a dict"),
        ):
            invalid = dict(valid)
            invalid[field] = value
            with self.subTest(field=field):
                with self.assertRaisesRegex(ValueError, error):
                    StageConfig.from_dict(invalid)

    def test_stage_config_reports_missing_required_fields(self):
        with self.assertRaisesRegex(ValueError, "Missing stage configuration"):
            StageConfig.from_dict({"name": "wide"})

    def test_stage_config_rejects_incomplete_split_ring_contracts(self):
        valid = {
            "name": "wide",
            "n_island": 20,
            "island_pop": 32,
            "evo_steps": 5,
            "sade_gen": 20,
            "ejection_model": "approximate",
            "initialization": "preselected",
            "algorithm": "sade",
            "topology": "split_ring_16_4",
            "split_ring": {
                "current_islands": 16,
                "alternative_islands": 4,
                "bridge_weight": 0.25,
            },
        }
        self.assertEqual(
            StageConfig.from_dict(valid).to_dict()["split_ring"],
            valid["split_ring"],
        )

        missing = dict(valid)
        missing.pop("split_ring")
        with self.assertRaisesRegex(ValueError, "requires split_ring options"):
            StageConfig.from_dict(missing)

        wrong_topology = dict(valid)
        wrong_topology["topology"] = "ring"
        with self.assertRaisesRegex(ValueError, "require split_ring_16_4"):
            StageConfig.from_dict(wrong_topology)

        too_small = dict(valid)
        too_small["n_island"] = 18
        too_small["split_ring"] = {
            "current_islands": 16,
            "alternative_islands": 2,
        }
        with self.assertRaisesRegex(ValueError, "at least three"):
            StageConfig.from_dict(too_small)

    def test_shared_benchmark_analysis_uses_consistent_statistics(self):
        self.assertEqual(
            case_from_sequence("Kerbin-Eve-Kerbin-Kerbin-Jool"),
            "KEKKJ",
        )
        summary = summarize_values([1.0, 3.0, 5.0], threshold=4.0)
        self.assertEqual(summary["median"], 3.0)
        self.assertEqual(summary["mean"], 3.0)
        self.assertEqual(summary["below_threshold"], 2)

        grouped = grouped_summary([
            {
                "strategy_label": "reference",
                "objective_dv": 1000.0,
                "runtime_seconds": 2.0,
            },
            {
                "strategy_label": "reference",
                "objective_dv": 1200.0,
                "runtime_seconds": 4.0,
            },
        ])
        self.assertEqual(grouped["reference"]["objective_dv"]["median"], 1100.0)
        self.assertEqual(grouped["reference"]["runtime_seconds"]["median"], 3.0)

    def test_ksp_calendar_dates_are_one_indexed(self):
        self.assertEqual(Kdate(0.0, "Vanilla"), "Y1 D1")
        self.assertEqual(Kdate(425.9, "Vanilla"), "Y1 D426")
        self.assertEqual(Kdate(426.0, "Vanilla"), "Y2 D1")
        self.assertEqual(Kdate(0.0, "JNSQ"), "Y1 D1")
        self.assertEqual(Kdate(365.0, "JNSQ"), "Y2 D1")

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

    def test_l0_recall_aliases_are_explicit_scout_archive_policies(self):
        default_config = OptimizationService.funnel_run_config(
            "funnel_l0_recall", 16, 32, 20, 50
        )
        legacy_config = OptimizationService.funnel_run_config(
            "funnel_scout_archive_64", 16, 32, 20, 50
        )

        self.assertEqual(default_config["optimizer_strategy"], "funnel_l0_recall")
        self.assertEqual(default_config["exact_strategy"], "scout_archive_nm_64")
        self.assertEqual(
            default_config["stages"], legacy_config["stages"]
        )
        self.assertEqual(default_config["stages"][0]["name"], "scout_unconnected")
        self.assertEqual(default_config["stages"][0]["topology"], "unconnected")
        self.assertEqual(default_config["stages"][0]["migration_rate"], 0)
        self.assertEqual(default_config["stages"][1]["initialization"], "scout_diverse")

        for suffix, expected_islands in (("32", 32), ("64", 64), ("128", 128)):
            with self.subTest(suffix=suffix):
                config = OptimizationService.funnel_run_config(
                    "funnel_l0_recall_" + suffix, 16, 32, 20, 50
                )
                self.assertEqual(
                    config["exact_strategy"], "scout_archive_nm_" + suffix
                )
                self.assertEqual(config["stages"][0]["n_island"], expected_islands)

    def test_requested_l0_strategy_name_is_preserved_in_results(self):
        self.assertEqual(
            Wayfinder._canonical_optimizer_strategy_from_exact_strategy(
                "scout_archive_nm_64"
            ),
            "funnel_scout_archive_64",
        )
        config = OptimizationService.funnel_run_config(
            "funnel_l0_recall", 16, 32, 20, 50
        )
        self.assertEqual(config["optimizer_strategy"], "funnel_l0_recall")

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

    def test_l0_recall_generation_sweep_strategies_override_scout_depth(self):
        for width in (32, 64, 128):
            for generations in (15, 30, 60):
                strategy = "funnel_l0_recall_{}_g{}".format(width, generations)
                exact_strategy = "scout_archive_nm_{}_g{}".format(
                    width, generations,
                )
                with self.subTest(strategy=strategy):
                    config = OptimizationService.funnel_run_config(
                        strategy, 16, 32, 20, 50,
                    )
                    scout = config["stages"][0]
                    self.assertEqual(config["exact_strategy"], exact_strategy)
                    self.assertEqual(scout["name"], "scout_unconnected")
                    self.assertEqual(scout["n_island"], width)
                    self.assertEqual(scout["sade_gen"], generations)
                    self.assertEqual(
                        Wayfinder._canonical_optimizer_strategy_from_exact_strategy(
                            exact_strategy,
                        ),
                        strategy,
                    )

    def test_l0_exact_selection_variants_rescore_scout_before_wide_stage(self):
        for strategy, exact_strategy in (
            ("funnel_l0_recall_128_exact", "scout_archive_nm_128_exact"),
            ("funnel_l0_recall_128_g30_exact", "scout_archive_nm_128_g30_exact"),
        ):
            with self.subTest(strategy=strategy):
                config = OptimizationService.funnel_run_config(
                    strategy, 16, 32, 20, 50,
                )
                self.assertEqual(config["exact_strategy"], exact_strategy)
                self.assertEqual(config["stages"][0]["name"], "scout_unconnected")
                self.assertEqual(config["stages"][1]["name"], "wide")
                self.assertEqual(config["stages"][1]["initialization"], "scout_diverse")
                self.assertEqual(config["stages"][1]["selection_problem"], "exact")
                self.assertEqual(
                    Wayfinder._canonical_optimizer_strategy_from_exact_strategy(
                        exact_strategy,
                    ),
                    strategy,
                )

    def test_l0_wide_override_reduces_scout_to_wide_compression(self):
        for strategy, exact_strategy in (
            ("funnel_l0_recall_128_w32", "scout_archive_nm_128_w32"),
            ("funnel_l0_recall_128_w32_exact", "scout_archive_nm_128_w32_exact"),
            ("funnel_l0_recall_128_g30_w32_exact", "scout_archive_nm_128_g30_w32_exact"),
        ):
            with self.subTest(strategy=strategy):
                config = OptimizationService.funnel_run_config(
                    strategy, 16, 32, 20, 50,
                )
                self.assertEqual(config["exact_strategy"], exact_strategy)
                self.assertEqual(
                    [stage["n_island"] for stage in config["stages"]],
                    [128, 32, 16, 8],
                )
                self.assertEqual(config["stages"][1]["initialization"], "scout_diverse")
                if strategy.endswith("_exact"):
                    self.assertEqual(config["stages"][1]["selection_problem"], "exact")
                self.assertEqual(
                    Wayfinder._canonical_optimizer_strategy_from_exact_strategy(
                        exact_strategy,
                    ),
                    strategy,
                )

    def test_l0_mbh_variants_insert_or_replace_intermediate_stage(self):
        between = OptimizationService.funnel_run_config(
            "funnel_l0_recall_64_mbh_between", 16, 32, 20, 50,
        )
        self.assertEqual(
            [stage["name"] for stage in between["stages"]],
            [
                "scout_unconnected", "wide", "mbh_refine",
                "intermediate", "exact_ejection",
            ],
        )
        self.assertEqual(between["stages"][2]["algorithm"], "mbh")
        self.assertEqual(between["stages"][2]["topology"], "unconnected")
        self.assertEqual(
            Wayfinder._canonical_optimizer_strategy_from_exact_strategy(
                "scout_archive_nm_64_mbh_between",
            ),
            "funnel_l0_recall_64_mbh_between",
        )

        replacement = OptimizationService.funnel_run_config(
            "funnel_l0_recall_64_mbh_l2", 16, 32, 20, 50,
        )
        self.assertEqual(
            [stage["name"] for stage in replacement["stages"]],
            ["scout_unconnected", "wide", "intermediate_mbh", "exact_ejection"],
        )
        self.assertEqual(replacement["stages"][2]["algorithm"], "mbh")
        self.assertEqual(replacement["stages"][2]["topology"], "unconnected")
        self.assertEqual(
            Wayfinder._canonical_optimizer_strategy_from_exact_strategy(
                "scout_archive_nm_64_mbh_l2",
            ),
            "funnel_l0_recall_64_mbh_l2",
        )

    def test_experimental_16_4_portfolio_is_named_and_serializable(self):
        strategy = (
            "funnel_l0_recall_64_mbh_between_"
            "portfolio_16_4_l0_pressure_cascade"
        )
        config = OptimizationService.funnel_run_config(
            strategy, 16, 32, 20, 50,
        )
        self.assertEqual(
            config["exact_strategy"],
            "scout_archive_nm_64_mbh_between_portfolio_16_4_l0",
        )
        self.assertEqual(config["optimizer_strategy"], strategy)
        self.assertTrue(config["pressure_cascade"]["enabled"])
        wide = config["stages"][1]
        self.assertEqual(wide["name"], "wide")
        self.assertEqual(wide["n_island"], 20)
        self.assertEqual(wide["initialization"], "preselected")
        self.assertEqual(wide["selection_policy"], "portfolio_16_4_l0")
        self.assertEqual(wide["topology"], "split_ring_16_4")
        self.assertEqual(wide["split_ring"], {
            "current_islands": 16,
            "alternative_islands": 4,
            "bridge_weight": 0.25,
        })
        self.assertEqual(
            Wayfinder._canonical_optimizer_strategy_from_exact_strategy(
                config["exact_strategy"],
            ),
            "funnel_l0_recall_64_mbh_between_portfolio_16_4_l0",
        )

    def test_pressure_cascade_strategy_keeps_base_funnel_and_adds_config(self):
        base = OptimizationService.funnel_run_config(
            "funnel_l0_recall_64_mbh_between", 16, 32, 20, 50,
        )
        cascade = OptimizationService.funnel_run_config(
            "funnel_l0_recall_64_mbh_between_pressure_cascade",
            16, 32, 20, 50,
        )

        self.assertEqual(
            cascade["exact_strategy"],
            "scout_archive_nm_64_mbh_between",
        )
        self.assertEqual(cascade["stages"], base["stages"])
        self.assertEqual(cascade["pressure_cascade"]["rescue_mode"], "cascade")
        self.assertEqual(cascade["pressure_cascade"]["branch_policy"], "l1_combined")
        self.assertEqual(cascade["pressure_cascade"]["top"], 32)
        self.assertNotIn("pressure_cascade", base)

    def test_pressure_helpers_detect_and_relax_high_tof_boundary(self):
        bounds = [[100.0, 200.0], [20.0, 60.0]]
        tof_vectors = [
            [198.0, 30.0],
            [199.0, 35.0],
            [197.5, 32.0],
            [196.0, 31.0],
        ]

        actions = OptimizationService.tof_boundary_pressure_actions(
            bounds, tof_vectors, near_fraction=0.03, min_pressure_count=2,
        )
        adjusted, deltas, modes = OptimizationService.adjust_leg_tof_bounds(
            bounds, actions, relax_fraction=0.20,
            min_relax_days=5.0, max_relax_days=60.0,
        )

        self.assertEqual(actions[0]["leg"], 1)
        self.assertEqual(actions[0]["side"], "high")
        self.assertEqual(adjusted[0], [100.0, 220.0])
        self.assertEqual(adjusted[1], [20.0, 60.0])
        self.assertEqual(deltas[(1, "high")], 20.0)
        self.assertEqual(modes[(1, "high")], "widen")

    def test_pressure_branch_and_retry_decisions_are_threshold_free(self):
        actions = [{"leg": 1, "side": "high"}]

        self.assertTrue(
            OptimizationService.pressure_branch_decision(
                actions,
                policy="l1_combined",
                l1_improvement=0.05,
                l1_best_to_l0_median=0.5,
            )
        )
        self.assertTrue(
            OptimizationService.pressure_branch_decision(
                actions,
                policy="l1_combined",
                l1_improvement=0.3,
                l1_best_to_l0_median=0.9,
            )
        )
        self.assertFalse(
            OptimizationService.pressure_branch_decision(
                actions,
                policy="l1_combined",
                l1_improvement=0.3,
                l1_best_to_l0_median=0.5,
            )
        )
        self.assertTrue(
            OptimizationService.pressure_retry_decision(
                baseline_dv=1000.0,
                current_best_dv=950.0,
                relaxed_dv=950.0,
                min_improvement=0.10,
            )
        )
        self.assertFalse(
            OptimizationService.pressure_retry_decision(
                baseline_dv=1000.0,
                current_best_dv=850.0,
                relaxed_dv=850.0,
                min_improvement=0.10,
            )
        )

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

    def test_handoff_embedding_balances_epoch_phase_tof_and_dsm(self):
        class LinearBody:
            def __init__(self, offset):
                self.offset = float(offset)

            def eph(self, epoch):
                value = float(epoch) + self.offset
                return [value, 1.0, 0.0], [0.0, value, 1.0]

        class BoundedProblem:
            @staticmethod
            def get_bounds():
                return (
                    [0.0, 0.0, 0.0, 0.0, 0.0, 1.0],
                    [20.0, 1.0, 1.0, 10.0, 1.0, 11.0],
                )

        bodies = {"A": LinearBody(1.0), "B": LinearBody(2.0)}
        context = {"bodies_json": '["A", "B"]', "tof_encoding": "direct"}
        gene = [10.0, 0.0, 0.0, 0.0, 0.25, 6.0]

        embedding = OptimizationService.handoff_diversity_embedding(
            bodies, context, BoundedProblem(), gene,
        )

        # t0 + 2 body phase states + one TOF + one DSM fraction.
        self.assertEqual(embedding.shape, (15,))
        self.assertAlmostEqual(embedding[0], 0.5)
        self.assertAlmostEqual(embedding[-2], 0.5)
        self.assertAlmostEqual(embedding[-1], 0.25)

    def test_pareto_handoff_keeps_champion_and_selects_distinct_basin(self):
        class FirstGeneFitness:
            @staticmethod
            def fitness(gene):
                return [float(gene[0])]

        genes = [
            [0.0, 0.0],
            [0.1, 0.01],
            [0.2, 1.0],
            [100.0, 100.0],
        ]
        selected = OptimizationService.select_pareto_diverse_elites(
            FirstGeneFitness(),
            genes,
            2,
            embedding_for_gene=lambda gene: np.asarray([gene[1]]),
            elite_fraction=0.75,
        )

        self.assertEqual(selected, [[0.0, 0.0], [0.2, 1.0]])
        self.assertNotIn([100.0, 100.0], selected)

    def test_pareto_handoff_strategy_is_explicit_and_composes_with_cascade(self):
        strategy = (
            "funnel_l0_recall_64_mbh_between_pareto_pressure_cascade"
        )
        config = OptimizationService.funnel_run_config(
            strategy, 16, 32, 20, 50,
        )

        self.assertEqual(
            config["exact_strategy"],
            "scout_archive_nm_64_mbh_between_pareto",
        )
        self.assertTrue(config["pressure_cascade"]["enabled"])
        handoff_stages = [
            stage for stage in config["stages"]
            if stage["initialization"] in (
                "scout_diverse", "phase_elites_mixed",
            )
        ]
        self.assertTrue(handoff_stages)
        self.assertTrue(all(
            stage["selection_policy"] == "pareto_all"
            for stage in handoff_stages
        ))
        self.assertEqual(
            Wayfinder._canonical_optimizer_strategy_from_exact_strategy(
                config["exact_strategy"],
            ),
            "funnel_l0_recall_64_mbh_between_pareto",
        )

    def test_pareto_l0_strategy_only_changes_scout_to_wide_handoff(self):
        strategy = (
            "funnel_l0_recall_64_mbh_between_pareto_l0_pressure_cascade"
        )
        config = OptimizationService.funnel_run_config(
            strategy, 16, 32, 20, 50,
        )

        self.assertEqual(
            config["exact_strategy"],
            "scout_archive_nm_64_mbh_between_pareto_l0",
        )
        wide = next(
            stage for stage in config["stages"] if stage["name"] == "wide"
        )
        self.assertEqual(wide["initialization"], "scout_diverse")
        self.assertEqual(wide["selection_policy"], "pareto_l0")
        later_handoffs = [
            stage for stage in config["stages"]
            if stage["initialization"] == "phase_elites_mixed"
        ]
        self.assertTrue(later_handoffs)
        self.assertTrue(all(
            "selection_policy" not in stage for stage in later_handoffs
        ))
        self.assertEqual(
            Wayfinder._canonical_optimizer_strategy_from_exact_strategy(
                config["exact_strategy"],
            ),
            "funnel_l0_recall_64_mbh_between_pareto_l0",
        )

    def test_basin_handoff_allocates_family_depth_and_keeps_each_basin(self):
        class FirstGeneFitness:
            @staticmethod
            def fitness(gene):
                return [float(gene[0])]

        genes = [
            [0.0, 0.00], [1.0, 0.01], [2.0, -0.01], [3.0, 0.02],
            [4.0, 10.0], [5.0, 10.1], [6.0, 9.9],
            [7.0, 20.0],
        ]
        selected = OptimizationService.select_basin_allocated_elites(
            FirstGeneFitness(),
            genes,
            6,
            embedding_for_gene=lambda gene: np.asarray([gene[1]]),
            elite_fraction=1.0,
            max_clusters=3,
        )

        near_zero = [gene for gene in selected if gene[1] < 5.0]
        near_ten = [gene for gene in selected if 5.0 <= gene[1] < 15.0]
        near_twenty = [gene for gene in selected if gene[1] >= 15.0]
        self.assertEqual(len(selected), 6)
        self.assertIn([0.0, 0.0], selected)
        self.assertGreater(len(near_zero), len(near_twenty))
        self.assertTrue(near_ten)
        self.assertTrue(near_twenty)

    def test_basin_l0_strategy_only_changes_scout_to_wide_handoff(self):
        strategy = (
            "funnel_l0_recall_64_mbh_between_basin_l0_pressure_cascade"
        )
        config = OptimizationService.funnel_run_config(
            strategy, 16, 32, 20, 50,
        )

        self.assertEqual(
            config["exact_strategy"],
            "scout_archive_nm_64_mbh_between_basin_l0",
        )
        wide = next(
            stage for stage in config["stages"] if stage["name"] == "wide"
        )
        self.assertEqual(wide["selection_policy"], "basin_l0")
        later_handoffs = [
            stage for stage in config["stages"]
            if stage["initialization"] == "phase_elites_mixed"
        ]
        self.assertTrue(all(
            "selection_policy" not in stage for stage in later_handoffs
        ))
        self.assertEqual(
            Wayfinder._canonical_optimizer_strategy_from_exact_strategy(
                config["exact_strategy"],
            ),
            "funnel_l0_recall_64_mbh_between_basin_l0",
        )

    def test_hill_valley_selection_separates_two_deep_valleys(self):
        class DoubleWell:
            @staticmethod
            def fitness(gene):
                value = float(gene[0])
                return [min((value + 1.0) ** 2, (value - 1.0) ** 2)]

            @staticmethod
            def get_bounds():
                return ([-2.0], [2.0])

            @staticmethod
            def get_nx():
                return 1

        selected, diagnostics = OptimizationService.select_hill_valley_elites(
            DoubleWell(),
            [[-1.0], [1.0], [-0.8], [0.8], [0.0]],
            2,
            elite_fraction=1.0,
            max_test_points=3,
            return_diagnostics=True,
        )

        self.assertEqual(selected, [[-1.0], [1.0]])
        self.assertGreaterEqual(diagnostics["cluster_count"], 2)
        self.assertGreater(diagnostics["hill_valley_evaluations"], 0)

    def test_hill_valley_periodic_interpolation_takes_short_arc(self):
        midpoint = OptimizationService._interpolate_gene(
            [3.0], [-3.0], 0.5,
            [-np.pi], [np.pi], periodic_indices=[0],
        )

        self.assertAlmostEqual(abs(midpoint[0]), np.pi)

    def test_hill_valley_l0_strategy_uses_quality_gated_full_l0_policy(self):
        strategy = (
            "funnel_l0_recall_64_mbh_between_hill_valley_l0_pressure_cascade"
        )
        config = OptimizationService.funnel_run_config(
            strategy, 16, 32, 20, 50,
        )

        self.assertEqual(
            config["exact_strategy"],
            "scout_archive_nm_64_mbh_between_hill_valley_l0",
        )
        wide = next(
            stage for stage in config["stages"] if stage["name"] == "wide"
        )
        self.assertEqual(wide["selection_policy"], "hill_valley_l0")
        self.assertEqual(wide["elite_fraction"], 0.35)
        self.assertEqual(
            Wayfinder._canonical_optimizer_strategy_from_exact_strategy(
                config["exact_strategy"],
            ),
            "funnel_l0_recall_64_mbh_between_hill_valley_l0",
        )

    def test_persistent_hill_valley_strategy_serializes_barrier_threshold(self):
        strategy = (
            "funnel_l0_recall_64_mbh_between_"
            "hill_valley_p2_l0_pressure_cascade"
        )
        config = OptimizationService.funnel_run_config(
            strategy, 16, 32, 20, 50,
        )

        self.assertEqual(
            config["exact_strategy"],
            "scout_archive_nm_64_mbh_between_hill_valley_p2_l0",
        )
        wide = next(
            stage for stage in config["stages"] if stage["name"] == "wide"
        )
        self.assertEqual(wide["selection_policy"], "hill_valley_p2_l0")
        self.assertEqual(wide["barrier_relative_tolerance"], 2.0)
        self.assertEqual(
            Wayfinder._canonical_optimizer_strategy_from_exact_strategy(
                config["exact_strategy"],
            ),
            "funnel_l0_recall_64_mbh_between_hill_valley_p2_l0",
        )

    def test_multiresolution_hill_valley_keeps_roots_and_family_depth(self):
        class DoubleWell:
            @staticmethod
            def fitness(gene):
                value = float(gene[0])
                return [min((value + 1.0) ** 2, (value - 1.0) ** 2)]

            @staticmethod
            def get_bounds():
                return ([-2.0], [2.0])

            @staticmethod
            def get_nx():
                return 1

        selected, diagnostics = OptimizationService.select_hill_valley_elites(
            DoubleWell(),
            [[-1.0], [1.0], [-0.8], [0.8], [0.0]],
            3,
            elite_fraction=1.0,
            valley_slot_fraction=0.5,
            return_diagnostics=True,
        )

        self.assertIn([-1.0], selected)
        self.assertIn([1.0], selected)
        self.assertEqual(diagnostics["valley_slots"], 2)
        self.assertEqual(diagnostics["family_slots"], 1)

    def test_multiresolution_hill_valley_strategy_serializes_slot_split(self):
        strategy = (
            "funnel_l0_recall_64_mbh_between_"
            "hill_valley_mr_l0_pressure_cascade"
        )
        config = OptimizationService.funnel_run_config(
            strategy, 16, 32, 20, 50,
        )

        self.assertEqual(
            config["exact_strategy"],
            "scout_archive_nm_64_mbh_between_hill_valley_mr_l0",
        )
        wide = next(
            stage for stage in config["stages"] if stage["name"] == "wide"
        )
        self.assertEqual(wide["selection_policy"], "hill_valley_mr_l0")
        self.assertEqual(wide["valley_slot_fraction"], 0.75)
        self.assertEqual(
            Wayfinder._canonical_optimizer_strategy_from_exact_strategy(
                config["exact_strategy"],
            ),
            "funnel_l0_recall_64_mbh_between_hill_valley_mr_l0",
        )

    def test_multiresolution_32_holds_wide_evaluation_budget_constant(self):
        current = OptimizationService.funnel_run_config(
            "funnel_l0_recall_64_mbh_between", 16, 32, 10, 20,
        )
        strategy = (
            "funnel_l0_recall_64_mbh_between_"
            "hill_valley_mr32_l0_pressure_cascade"
        )
        config = OptimizationService.funnel_run_config(
            strategy, 16, 32, 10, 20,
        )
        current_wide = next(
            stage for stage in current["stages"] if stage["name"] == "wide"
        )
        wide = next(
            stage for stage in config["stages"] if stage["name"] == "wide"
        )

        self.assertEqual(wide["selection_policy"], "hill_valley_mr32_l0")
        self.assertEqual((wide["n_island"], wide["island_pop"]), (32, 16))
        self.assertEqual(
            current_wide["n_island"] * current_wide["island_pop"]
            * current_wide["evo_steps"] * current_wide["sade_gen"],
            wide["n_island"] * wide["island_pop"]
            * wide["evo_steps"] * wide["sade_gen"],
        )
        self.assertEqual(
            Wayfinder._canonical_optimizer_strategy_from_exact_strategy(
                config["exact_strategy"],
            ),
            "funnel_l0_recall_64_mbh_between_hill_valley_mr32_l0",
        )

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

    def test_tisserand_helpers_match_circular_patched_conic_relations(self):
        plans = Wayfinder(planet_pack="Vanilla")
        scout = SequenceScout(plans._planet_pack_module)
        body = scout._body("Kerbin")
        vinf = 1200.0
        state = scout._orbit_from_relative_state(
            "Kerbin", vinf, math.radians(37.0), ("Kerbin",), vinf, 0.0, (),
        )
        self.assertIsNotNone(state)
        geometric = tisserand_parameter(
            state.semi_major_axis,
            state.eccentricity,
            body.orbit_radius,
        )
        velocity_form = 3.0 - (vinf / body.orbital_speed) ** 2
        self.assertAlmostEqual(geometric, velocity_form, places=10)

        maximum_turn = maximum_flyby_turn_angle(
            vinf, body.mu_self, body.safe_radius,
        )
        required_radius = required_flyby_periapsis(
            vinf, maximum_turn, body.mu_self,
        )
        self.assertAlmostEqual(required_radius, body.safe_radius, places=5)
        self.assertGreater(
            maximum_turn,
            maximum_flyby_turn_angle(
                2.0 * vinf, body.mu_self, body.safe_radius,
            ),
        )

    def test_tisserand_tree_recovers_known_kekkj_and_keemo_sequences(self):
        vanilla = Wayfinder(planet_pack="Vanilla")
        jnsq = Wayfinder(planet_pack="JNSQ")

        jool_candidates = vanilla.scout_sequences("Kerbin", "Jool")
        jool_sequences = [candidate.sequence for candidate in jool_candidates]
        self.assertIn(
            ("Kerbin", "Eve", "Kerbin", "Kerbin", "Jool"),
            jool_sequences,
        )
        self.assertNotIn(
            "Eeloo",
            {
                body
                for candidate in jool_candidates
                for body in candidate.sequence[1:-1]
            },
        )
        kekkj = next(
            candidate for candidate in jool_candidates
            if candidate.sequence
            == ("Kerbin", "Eve", "Kerbin", "Kerbin", "Jool")
        )
        self.assertTrue(all(
            flyby.operational_periapsis_radius_m
            == flyby.safe_periapsis_radius_m
            for flyby in kekkj.flybys
        ))
        self.assertTrue(all(
            flyby.required_periapsis_radius_m is None
            or flyby.required_periapsis_radius_m
            >= flyby.operational_periapsis_radius_m - 1e-6
            for flyby in kekkj.flybys
        ))

        moho_candidates = jnsq.scout_sequences(
            "Kerbin", "Moho", limit=10, as_dict=True,
        )
        self.assertIn(
            ["Kerbin", "Eve", "Eve", "Moho"],
            [candidate["sequence"] for candidate in moho_candidates],
        )
        self.assertTrue(all(
            candidate["model"] == "tisserand_circular_coplanar_unphased"
            for candidate in moho_candidates
        ))
        self.assertNotIn("Infinity", json.dumps(moho_candidates))

    def test_tisserand_scout_config_is_validated_and_serializable(self):
        config = TisserandScoutConfig(
            max_bodies=4,
            departure_vinf_samples=7,
            states_per_sequence=32,
        )
        self.assertEqual(config.to_dict()["states_per_sequence"], 32)
        with self.assertRaisesRegex(ValueError, "max_bodies"):
            TisserandScoutConfig(max_bodies=1)
        with self.assertRaisesRegex(ValueError, "bounds must be increasing"):
            TisserandScoutConfig(
                departure_vinf_min_mps=2000.0,
                departure_vinf_max_mps=1000.0,
            )
        with self.assertRaisesRegex(ValueError, "flyby_clearance_m"):
            TisserandScoutConfig(flyby_clearance_m=-1.0)

    def test_lambert_arc1_filter_keeps_known_kekkj_without_manual_dv_limit(self):
        plans = Wayfinder(planet_pack="Vanilla")
        candidates = plans.scout_sequences("Kerbin", "Jool")
        assessments = plans.filter_scout_sequences_arc1(
            candidates,
            t0_bounds_days=[0.0, 2000.0],
            accepted_only=False,
        )
        accepted = [assessment for assessment in assessments if assessment.accepted]
        kekkj = next(
            assessment for assessment in assessments
            if assessment.candidate.sequence
            == ("Kerbin", "Eve", "Kerbin", "Kerbin", "Jool")
        )

        self.assertGreater(len(assessments), len(accepted))
        self.assertGreater(len(accepted), 0)
        self.assertTrue(kekkj.accepted)
        self.assertEqual(kekkj.reason, "accepted")
        self.assertLess(kekkj.departure_ratio, 1.0)
        self.assertLess(
            kekkj.tisserand_vinf_relative_error,
            LambertArc1Config().maximum_tisserand_vinf_relative_error,
        )
        self.assertNotIn("Infinity", json.dumps(kekkj.to_dict()))

        direct = LambertArc1Filter(
            plans._planet_pack_module,
        )._direct_reference("Kerbin", "Jool", [0.0, 2000.0])
        self.assertAlmostEqual(
            kekkj.direct_reference_vinf_mps,
            direct.departure_vinf_mps,
            places=9,
        )
        self.assertAlmostEqual(
            kekkj.direct_reference_ejection_dv_mps,
            direct.ejection_dv_mps,
            places=9,
        )

        with self.assertRaisesRegex(ValueError, "t0_samples"):
            LambertArc1Config(t0_samples=1)

    def test_t0_bin_scanner_uses_global_direct_ejection_reference(self):
        plans = Wayfinder(planet_pack="Vanilla")
        candidates = plans.scout_sequences("Kerbin", "Jool")
        kekkj = next(
            candidate for candidate in candidates
            if candidate.sequence
            == ("Kerbin", "Eve", "Kerbin", "Kerbin", "Jool")
        )
        result = plans.scan_scout_sequence_bins(
            [kekkj],
            [500.0, 700.0],
            bin_width_days=100.0,
            t0_step_days=25.0,
            tof_samples=16,
            maximum_ejection_ratio=1.05,
        )

        self.assertEqual(result.bin_width_days, 100.0)
        self.assertGreater(result.direct_reference.ejection_dv_mps, 0.0)
        self.assertGreater(len(result.candidates), 0)
        self.assertTrue(all(
            row.direct_reference == result.direct_reference
            for row in result.candidates
        ))
        self.assertTrue(all(
            row.ejection_ratio <= 1.05 for row in result.candidates
        ))
        self.assertEqual(result.ranked_unique(1)[0].candidate.sequence, kekkj.sequence)
        self.assertNotIn("Infinity", json.dumps(result.to_dict()))

    def test_vector_ejection_api_matches_gene_wrapper(self):
        plans = Wayfinder(planet_pack="Vanilla")
        body = plans._fullname_dic["Kerbin"]
        gene = [0.0, 0.17, 0.63, 1250.0]
        from_gene = fast_ejection_from_gene(body, gene, 100000.0)
        from_vector = fast_ejection_from_vinf(
            body, from_gene["vinf"], 100000.0,
        )
        self.assertAlmostEqual(from_gene["dv"], from_vector["dv"], places=10)
        self.assertEqual(from_gene["strategy"], from_vector["strategy"])

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
            expected_objective=2616.6040426094487,
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

    def test_fast_ejection_uses_finite_soi_geometry(self):
        plans = load_wayfinder_reference("Ex1_KEEMo_MGA", planet_pack="Vanilla")
        job_index = ("KEEMo", 0, 400)
        row = plans._reference_rows[job_index]
        udp = reference_udp(plans, row)
        gene = json.loads(row["gene_json"])
        rsoi = plans.soi_radius(udp._seq[0])
        _, departure_velocity = udp._seq[0].eph(pk.epoch(gene[0]))
        prograde = departure_velocity / np.linalg.norm(departure_velocity)
        _, vinfx, vinfy, vinfz = udp._decode_times_and_vinf(gene)

        legacy_dv, legacy_inclination, _ = EJ_details(
            udp, [vinfx, vinfy, vinfz], rsoi, prograde, row["ejection_altitude"]
        )
        ejection = fast_ejection_from_gene(
            udp._seq[0], gene, row["ejection_altitude"], rsoi
        )

        self.assertAlmostEqual(ejection["direct_dv"], legacy_dv[0], places=9)
        self.assertAlmostEqual(
            ejection["direct_inclination"], legacy_inclination, places=9
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
        self.assertAlmostEqual(
            ejection["selected_normal_dv"],
            ejection["direct_normal_dv"]
            if ejection["strategy"] == "direct"
            else ejection["split_normal_dv"],
        )
        self.assertEqual(
            ejection["selected_normal_dv_location"],
            "parking_orbit" if ejection["strategy"] == "direct" else "soi",
        )

    def test_fast_ejection_reports_selected_z_component_source(self):
        plans = Wayfinder(planet_pack="Vanilla")
        kerbin = plans._fullname_dic["Kerbin"]
        altitude = 100000.0

        direct_gene = [0.0, 0.25, 0.05, 5000.0]
        direct = fast_ejection_from_gene(kerbin, direct_gene, altitude)
        self.assertEqual(direct["strategy"], "direct")
        self.assertEqual(direct["selected_normal_dv_location"], "parking_orbit")
        self.assertAlmostEqual(
            direct["selected_normal_dv"], direct["direct_normal_dv"]
        )
        self.assertNotAlmostEqual(
            abs(direct["selected_normal_dv"]),
            direct["soi_correction_dv"],
        )

        split_gene = [0.0, 0.0, 0.001, 8000.0]
        split = fast_ejection_from_gene(kerbin, split_gene, altitude)
        self.assertEqual(split["strategy"], "split_soi")
        self.assertEqual(split["selected_normal_dv_location"], "soi")
        self.assertAlmostEqual(
            split["selected_normal_dv"], split["split_normal_dv"]
        )
        self.assertAlmostEqual(
            abs(split["selected_normal_dv"]), split["soi_correction_dv"]
        )

    def test_lambert_kd_parking_node_components_match_ksp_validation(self):
        plans = Wayfinder(planet_pack="Vanilla")
        kerbin = plans._fullname_dic["Kerbin"]
        duna = plans._fullname_dic["Duna"]
        # Empirically validated in KSP after the finite-SOI ejection fix:
        # pure zero-rev Lambert, no DSM, 100 km Kerbin parking orbit.
        t0_kut = 228.77461294
        tof_kut = 239.85192938
        t0 = t0_kut / 4.0
        tof = tof_kut / 4.0
        r0, v0 = kerbin.eph(pk.epoch(t0))
        r1, _ = duna.eph(pk.epoch(t0 + tof))
        lambert = lambert_problem(
            r0, r1, tof * DAY2SEC, kerbin.mu_central_body,
            cw=False, multi_revs=0,
        )
        vinf = np.asarray(_lambert_v1(lambert)) - np.asarray(v0)
        vinf_norm = float(np.linalg.norm(vinf))
        theta = math.atan2(vinf[1], vinf[0]) % (2 * math.pi)
        gene = [
            t0,
            theta / (2 * math.pi),
            (1.0 - float(vinf[2] / vinf_norm)) / 2.0,
            vinf_norm,
            0.0,
            1.0,
        ]

        ejection = fast_ejection_from_gene(
            kerbin, gene, 100000.0, plans.soi_radius(kerbin)
        )

        self.assertEqual(ejection["strategy"], "direct")
        self.assertAlmostEqual(ejection["parking_node_dv"], 1038.315881650874, places=6)
        self.assertAlmostEqual(
            ejection["parking_node_prograde_dv"], 1031.932230848249, places=6
        )
        self.assertAlmostEqual(
            ejection["parking_node_normal_dv"], 114.95973653844, places=6
        )
        self.assertAlmostEqual(ejection["parking_node_radial_dv"], 0.0, places=6)

    def test_decode_reports_soi_timing_diagnostics(self):
        plans = load_wayfinder_reference("Ex1_KEEMo_MGA", planet_pack="Vanilla")
        decoded = decode_trajectory(plans, ("KEEMo", 0, 400))

        departure_offset = decoded["departure_soi_escape_time"]
        arrival_offset = decoded["arrival_soi_to_periapsis_time"]

        self.assertGreater(departure_offset, 0.0)
        self.assertGreater(arrival_offset, 0.0)
        self.assertGreater(decoded["heliocentric_arrival_met"], 0.0)
        self.assertIn("arrival_periapsis_altitude_estimate", decoded)

    def test_ksp_translation_rescores_kd_dsm_from_finite_soi_state(self):
        plans = Wayfinder(planet_pack="Vanilla")
        udp = _make_mga_1dsm(
            seq=[plans._fullname_dic["Kerbin"], plans._fullname_dic["Duna"]],
            t0=[50.0, 60.0],
            tof=[[50.0, 80.0]],
            vinf=[0.0, 2000.0],
            tof_encoding="direct",
            add_vinf_dep=False,
            add_vinf_arr=False,
            multi_objective=False,
            orbit_insertion=False,
            rp_target=0.0,
            e_target=0.0,
            max_revs=0,
        )
        gene = [
            57.06036947959551,
            0.28506896814935573,
            0.49717490891777705,
            853.582384564183,
            0.39830364563986964,
            64.58027340707079,
        ]

        translation = translate_first_leg(
            udp, gene, planet_pack="Vanilla", parking_altitude=100000.0
        )

        self.assertAlmostEqual(
            translation.arrival_met_days, 258.32109362828317, places=6
        )
        self.assertAlmostEqual(
            translation.arrival_epoch_days, 486.5625715466652, places=6
        )
        self.assertAlmostEqual(
            translation.dsm_met_days, 102.89023333782329, places=6
        )
        self.assertAlmostEqual(
            translation.ideal_dsm.magnitude, 5.160017441050398, places=6
        )
        self.assertAlmostEqual(
            translation.selected_case.corrected_dsm.magnitude,
            22.165309242180122,
            places=6,
        )
        self.assertAlmostEqual(
            translation.split_case.corrected_dsm.magnitude,
            22.096647817623474,
            places=6,
        )
        self.assertGreater(
            translation.selected_case.ideal_to_corrected_dsm_delta, 15.0
        )

        selected_plan = build_first_leg_flight_plan(
            udp, gene, planet_pack="Vanilla", strategy="selected",
            corrected_dsm=True,
        )
        self.assertEqual([node.label for node in selected_plan.nodes], [
            "departure", "dsm_corrected",
        ])
        self.assertAlmostEqual(
            selected_plan.nodes[0].ut_seconds, 228.24147791838204 * 21600,
            places=3,
        )
        self.assertAlmostEqual(
            selected_plan.nodes[1].prograde, 2.572997624565908, places=6
        )
        self.assertAlmostEqual(
            selected_plan.nodes[1].normal, 5.1394020922081, places=6
        )
        self.assertAlmostEqual(
            selected_plan.nodes[1].radial, -21.407175506359096, places=6
        )

        split_plan = build_first_leg_flight_plan(
            udp, gene, planet_pack="Vanilla", strategy="split_soi",
            corrected_dsm=True,
        )
        self.assertEqual([node.label for node in split_plan.nodes], [
            "departure_planar", "soi_normal_correction", "dsm_corrected",
        ])
        self.assertAlmostEqual(
            split_plan.nodes[1].met_days, 4.092827324854004, places=6
        )

    def test_vanilla_kekkj_reference_decodes_from_sqlite(self):
        plans = load_wayfinder_reference("Ex2_KEKKJ_MGA", planet_pack="Vanilla")

        self.assertEqual(list(plans._reference_rows), [("KEKKJ", 600, 3000)])
        assert_decoded_matches_stored(
            self, plans, ("KEKKJ", 600, 3000),
            expected_objective=1788.9615664299477,
        )

    def test_vanilla_kekkj_ksp_plan_exposes_full_timeline_with_fidelity_tags(self):
        plans = load_wayfinder_reference("Ex2_KEKKJ_MGA", planet_pack="Vanilla")
        job_index = ("KEKKJ", 600, 3000)
        row = plans._reference_rows[job_index]
        udp = reference_udp(plans, row)
        gene = json.loads(row["gene_json"])

        translation = translate_mga_trajectory(
            udp,
            gene,
            planet_pack="Vanilla",
            parking_altitude=row["ejection_altitude"],
            arrival_mode="elliptical",
        )
        plan = ksp_plan_for_reference(plans, job_index, strategy="split_soi")

        self.assertEqual(len(translation.flybys), 3)
        self.assertEqual([flyby.body for flyby in translation.flybys], [
            "Eve", "Kerbin", "Kerbin",
        ])
        self.assertTrue(
            all(flyby.model == "ksp_flyby_translated"
                for flyby in translation.flybys)
        )
        self.assertTrue(
            all(flyby.soi_entry_epoch_days < flyby.epoch_days
                < flyby.soi_exit_epoch_days
                for flyby in translation.flybys)
        )
        self.assertTrue(
            all(0.0 <= flyby.hyperbola_inclination_deg <= 180.0
                for flyby in translation.flybys)
        )
        self.assertTrue(
            all(len(flyby.flyby_plane_normal) == 3
                and len(flyby.periapsis_direction) == 3
                for flyby in translation.flybys)
        )
        self.assertEqual(len(plan.flybys), 3)
        self.assertTrue(
            all(0.0 <= flyby.hyperbola_inclination_deg <= 180.0
                for flyby in plan.flybys)
        )
        self.assertEqual(len(plan.nodes), 6)
        post_flyby_nodes = [node for node in plan.nodes if node.leg_index != 1]
        self.assertEqual(len(post_flyby_nodes), 3)
        self.assertTrue(
            all(node.model == "pykep_ideal_fallback"
                for node in post_flyby_nodes)
        )
        self.assertFalse(
            any(node.model.startswith("ksp_") for node in post_flyby_nodes)
        )
        self.assertEqual(plan.arrival.arrival_mode, "elliptical")
        self.assertEqual(plan.arrival.model, "pykep_arrival_fallback")
        self.assertGreater(plan.arrival.arrival_vinf, 1000.0)

    def test_jnsq_moho_reference_decodes_with_jnsq_time_scale(self):
        plans = load_wayfinder_reference("JNSQ_Moho_reference_fixed", planet_pack="JNSQ")

        self.assertEqual(len(plans._reference_rows), 8)

        expected_objectives = {
            ("KEMo", 200, 450): 5845.92533439078,
            ("KEEMo", 200, 450): 5910.310593318254,
            ("KEMo", 200, 300): 6558.101589910543,
            ("KEEMo", 0, 300): 6730.108872627146,
            ("KEMo", 0, 450): 6900.262227506135,
            ("KEMo", 0, 300): 7099.316352010953,
            ("KEEMo", 0, 450): 7338.7055346986,
            ("KEEMo", 200, 300): 7668.794071016291,
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
                tof_profile="relaxed",
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

    def test_direct_batch_accepts_custom_leg_tof_bounds(self):
        plans = Wayfinder(planet_pack="Vanilla", datastore_name="direct-test")
        template = [["Kerbin"], ["Duna"]]
        custom_bounds = [[220.0, 280.0]]
        with tempfile.TemporaryDirectory() as directory:
            db_path = Path(directory) / "direct_custom.sqlite"
            plans.add_direct_t0_batch_sqlite(
                template, db_path, batch_name="KD-custom",
                t0_min=200, t0_bin=20, n_t0_bins=1,
                tof_profile="relaxed",
                leg_tof_bounds_override=custom_bounds,
            )
            store = SQLiteJobStore(db_path)
            try:
                jobs = store.pending_jobs("Vanilla", limit=10, batch_name="KD-custom")
            finally:
                store.close()

        self.assertEqual(len(jobs), 1)
        self.assertEqual(json.loads(jobs[0]["leg_tof_bounds_json"]), custom_bounds)
        self.assertEqual(jobs[0]["tof_min"], 220.0)
        self.assertEqual(jobs[0]["tof_max"], 280.0)
        self.assertEqual(jobs[0]["tof_profile"], "relaxed_custom_bounds")

    def test_direct_batch_arrival_modes_have_distinct_objectives(self):
        plans = Wayfinder(planet_pack="Vanilla", datastore_name="arrival-mode-test")
        template = [["Kerbin"], ["Duna"]]
        expected = {
            "flyby": (False, False, 0.0),
            "none": (False, False, 0.0),
            "vinf": (True, False, 0.0),
            "circular": (True, True, 0.0),
            "elliptical": (True, True, 0.9),
        }
        with tempfile.TemporaryDirectory() as directory:
            db_path = Path(directory) / "arrival_modes.sqlite"
            for mode in expected:
                plans.add_direct_t0_batch_sqlite(
                    template, db_path, batch_name=f"mode-{mode}",
                    t0_min=0, t0_bin=50, n_t0_bins=1,
                    arrival_mode=mode,
                )

            store = SQLiteJobStore(db_path)
            try:
                jobs_by_mode = {}
                for mode in expected:
                    jobs = store.pending_jobs(
                        "Vanilla", limit=10, batch_name=f"mode-{mode}"
                    )
                    self.assertEqual(len(jobs), 1)
                    jobs_by_mode[mode] = jobs[0]
            finally:
                store.close()

        for mode, (add_vinf_arr, orbit_insertion, e_target) in expected.items():
            job = jobs_by_mode[mode]
            canonical_mode = "flyby" if mode == "none" else mode
            self.assertEqual(job["arrival_mode"], canonical_mode)
            self.assertEqual(job["opt_injection"], canonical_mode)
            self.assertEqual(bool(job["add_vinf_arr"]), add_vinf_arr)
            self.assertEqual(bool(job["orbit_insertion"]), orbit_insertion)
            self.assertEqual(float(job["e_target"]), e_target)

    def test_legacy_alpha_batch_persists_canonical_arrival_mode(self):
        plans = Wayfinder(planet_pack="Vanilla", datastore_name="alpha-arrival-test")
        template = [["Kerbin"], ["Duna"]]
        with tempfile.TemporaryDirectory() as directory:
            db_path = Path(directory) / "alpha_arrival.sqlite"
            plans.add_batch_sqlite(
                template,
                db_path=db_path,
                batch_name="alpha-flyby",
                t0_min=0,
                t0_bin=50,
                n_t0_bins=1,
                tof_min=200,
                tof_bin=50,
                n_tof_bins=1,
                opt_injection="flyby",
            )
            store = SQLiteJobStore(db_path)
            try:
                jobs = store.pending_jobs(
                    "Vanilla", limit=10, batch_name="alpha-flyby"
                )
            finally:
                store.close()

        self.assertEqual(len(jobs), 1)
        self.assertEqual(jobs[0]["opt_injection"], "flyby")
        self.assertEqual(jobs[0]["arrival_mode"], "flyby")
        self.assertFalse(bool(jobs[0]["add_vinf_arr"]))
        self.assertFalse(bool(jobs[0]["orbit_insertion"]))


if __name__ == "__main__":
    unittest.main(verbosity=2)
