# -*- coding: utf-8 -*-
"""SQLite datastore tests for cross-batch Wayfinder queries."""

import os
import sys
import contextlib
import io
import shutil
import unittest
import uuid
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt


ROOT = Path(__file__).resolve().parents[1]
CORE_DIR = ROOT / "WayfinderCore"
TESTS_DIR = ROOT / "Tests"
REFERENCE_DB = TESTS_DIR / "wayfinder_reference.sqlite"

sys.path.insert(0, str(CORE_DIR))

from _SQLiteStore import SQLiteJobStore  # noqa: E402
from _Wayfinder import Wayfinder  # noqa: E402


class FakePopulation:
    def __init__(self, genes, fitnesses):
        self._genes = genes
        self._fitnesses = fitnesses

    def get_x(self):
        return self._genes

    def get_f(self):
        return self._fitnesses


class SQLiteDatastoreTests(unittest.TestCase):
    def test_add_batch_sqlite_generates_wildcard_jobs(self):
        db_path = TESTS_DIR / f"wayfinder_sqlite_test_{uuid.uuid4().hex}.sqlite"
        plans = Wayfinder(planet_pack="Vanilla")

        plans.add_batch_sqlite(
            db_path=db_path,
            batch_name="wildcard_sqlite_only",
            swing_by_bodies=[["Kerbin"], ["Eve"], ["Eve", "*"], ["Moho"]],
            t0_min=0,
            t0_bin=100,
            n_t0_bins=2,
            auto_tof=True,
            opt_level="debug",
            opt_injection="vinf",
        )

        store = SQLiteJobStore(db_path)
        try:
            self.assertEqual(store.count_rows("batches"), 1)
            self.assertEqual(store.count_rows("sequences"), 2)
            self.assertEqual(store.count_rows("jobs"), 4)
            self.assertEqual(store.count_rows("results"), 0)

            route_jobs = store.best_results(
                planet_pack="Vanilla",
                start_body="Kerbin",
                target_body="Moho",
                limit=10,
            )
            self.assertEqual(route_jobs, [])
        finally:
            store.close()
            db_path.unlink(missing_ok=True)

    def test_sqlite_jobs_are_loaded_from_store_not_dataframe_staging(self):
        db_path = TESTS_DIR / f"wayfinder_sqlite_test_{uuid.uuid4().hex}.sqlite"
        plans = Wayfinder(planet_pack="Vanilla")
        plans.add_batch_sqlite(
            db_path=db_path,
            batch_name="sqlite_load_only",
            swing_by_bodies=[["Kerbin"], ["Eve"], ["Moho"]],
            t0_min=0,
            t0_bin=100,
            n_t0_bins=1,
            auto_tof=True,
            opt_level="debug",
            opt_injection="vinf",
        )

        try:
            with self.assertRaises(NotImplementedError):
                plans.load_sqlite_jobs(db_path, limit=1, batch_name="sqlite_load_only")

            store = SQLiteJobStore(db_path)
            try:
                jobs = store.pending_jobs("Vanilla", limit=1, batch_name="sqlite_load_only")
            finally:
                store.close()

            self.assertEqual(len(jobs), 1)
            self.assertEqual(jobs[0]["sequence_short_name"], "KEMo")
            self.assertEqual(jobs[0]["bodies_json"], '["Kerbin","Eve","Moho"]')
            self.assertEqual(jobs[0]["status"], "TODO")
        finally:
            db_path.unlink(missing_ok=True)

    def test_optimizer_snapshots_are_stored_for_running_run(self):
        db_path = TESTS_DIR / f"wayfinder_sqlite_test_{uuid.uuid4().hex}.sqlite"
        plans = Wayfinder(planet_pack="Vanilla")
        plans.add_batch_sqlite(
            db_path=db_path,
            batch_name="snapshot_store",
            swing_by_bodies=[["Kerbin"], ["Eve"], ["Moho"]],
            t0_min=0,
            t0_bin=100,
            n_t0_bins=1,
            auto_tof=True,
            opt_level="debug",
            opt_injection="vinf",
        )

        store = SQLiteJobStore(db_path)
        try:
            job = store.pending_jobs("Vanilla", limit=1, batch_name="snapshot_store")[0]
            run_id = store.start_run(job["job_id"], versions={"wayfinder": "test"})
            store.record_optimizer_snapshot(
                run_id,
                step=1,
                champion_fitness=[[123.4]],
                champion_genes=[[1.0, 2.0, 3.0]],
            )
            store.record_optimizer_population(
                run_id,
                step=1,
                populations=[
                    FakePopulation(
                        genes=[[1.0, 2.0], [3.0, 4.0]],
                        fitnesses=[[123.4], [456.7]],
                    )
                ],
                source="final",
            )
            self.assertEqual(store.count_rows("optimizer_snapshots"), 1)
            self.assertEqual(store.count_rows("optimizer_population_points"), 2)
            population_points = store.optimizer_population_points(run_id, source="final")
            self.assertEqual(population_points[1]["gene"], [3.0, 4.0])
        finally:
            store.close()
            db_path.unlink(missing_ok=True)

    def test_benchmark_batches_are_excluded_from_default_queries(self):
        db_path = TESTS_DIR / f"wayfinder_sqlite_test_{uuid.uuid4().hex}.sqlite"
        plans = Wayfinder(planet_pack="Vanilla")
        common_kwargs = {
            "db_path": db_path,
            "swing_by_bodies": [["Kerbin"], ["Eve"], ["Moho"]],
            "t0_bin": 100,
            "n_t0_bins": 1,
            "auto_tof": True,
            "opt_level": "debug",
            "opt_injection": "vinf",
        }
        plans.add_batch_sqlite(
            batch_name="production_shared",
            t0_min=0,
            purpose="production",
            **common_kwargs,
        )
        plans.add_batch_sqlite(
            batch_name="benchmark_shared",
            t0_min=0,
            purpose="benchmark",
            **common_kwargs,
        )
        plans.add_batch_sqlite(
            batch_name="benchmark_only",
            t0_min=100,
            purpose="benchmark",
            **common_kwargs,
        )

        store = SQLiteJobStore(db_path)
        try:
            rows = store.conn.execute(
                "SELECT name, purpose FROM batches ORDER BY name"
            ).fetchall()
            self.assertEqual(
                [(row["name"], row["purpose"]) for row in rows],
                [
                    ("benchmark_only", "benchmark"),
                    ("benchmark_shared", "benchmark"),
                    ("production_shared", "production"),
                ],
            )

            jobs = {
                float(job["t0_min"]): job
                for job in store.pending_jobs("Vanilla", limit=10)
            }
            self.assertEqual(len(jobs), 2)
            store.update_job_status(jobs[0.0]["job_id"], "DONE")
            store.upsert_result(
                jobs[0.0]["job_id"],
                objective_dv=6000.0,
                result_t0=10.0,
                result_tof=700.0,
                ejection_vinf=1.2,
                gene=[0.0, 1.0, 2.0],
            )
            store.update_job_status(jobs[100.0]["job_id"], "DONE")
            store.upsert_result(
                jobs[100.0]["job_id"],
                objective_dv=4000.0,
                result_t0=110.0,
                result_tof=650.0,
                ejection_vinf=1.1,
                gene=[3.0, 4.0, 5.0],
            )

            default_best = store.best_results(
                planet_pack="Vanilla",
                start_body="Kerbin",
                target_body="Moho",
                limit=1,
            )
            self.assertEqual(default_best[0]["t0_min"], 0.0)
            self.assertEqual(default_best[0]["objective_dv"], 6000.0)

            with_benchmarks = store.best_results(
                planet_pack="Vanilla",
                start_body="Kerbin",
                target_body="Moho",
                include_benchmarks=True,
                limit=1,
            )
            self.assertEqual(with_benchmarks[0]["t0_min"], 100.0)
            self.assertEqual(with_benchmarks[0]["objective_dv"], 4000.0)

            explicit_benchmark = store.best_results(
                planet_pack="Vanilla",
                batch_name="benchmark_only",
                limit=1,
            )
            self.assertEqual(explicit_benchmark[0]["t0_min"], 100.0)

            benchmark_rows = store.benchmark_results(
                planet_pack="Vanilla",
                benchmark_name="benchmark",
            )
            self.assertEqual(
                [(row["batch_name"], row["objective_dv"]) for row in benchmark_rows],
                [("benchmark_only", 4000.0), ("benchmark_shared", 6000.0)],
            )

            default_rows = store.result_rows(
                planet_pack="Vanilla",
                start_body="Kerbin",
                target_body="Moho",
            )
            self.assertEqual([row["t0_min"] for row in default_rows], [0.0])
        finally:
            store.close()
            db_path.unlink(missing_ok=True)

    def test_topology_benchmark_suite_creates_distinct_jobs(self):
        db_path = TESTS_DIR / f"wayfinder_sqlite_test_{uuid.uuid4().hex}.sqlite"
        plans = Wayfinder(planet_pack="Vanilla")

        created = plans.add_topology_benchmark_sqlite(
            db_path=db_path,
            benchmark_name="topology_smoke",
            swing_by_bodies=[["Kerbin"], ["Eve"], ["Moho"]],
            topologies=("unconnected", "ring"),
            seeds=[11, 12],
            t0_min=0,
            t0_bin=100,
            n_t0_bins=1,
            auto_tof=True,
            opt_level="debug",
            opt_injection="vinf",
        )

        store = SQLiteJobStore(db_path)
        try:
            self.assertEqual(len(created), 4)
            self.assertEqual(store.count_rows("batches"), 4)
            self.assertEqual(store.count_rows("jobs"), 4)

            rows = store.conn.execute(
                """
                SELECT b.name, b.purpose, j.optimizer_topology, j.optimizer_seed, j.status
                FROM batches b
                JOIN batch_jobs bj ON bj.batch_id = b.id
                JOIN jobs j ON j.id = bj.job_id
                ORDER BY b.name
                """
            ).fetchall()
            self.assertEqual({row["purpose"] for row in rows}, {"benchmark"})
            self.assertEqual(
                {(row["optimizer_topology"], row["optimizer_seed"]) for row in rows},
                {("unconnected", 11), ("unconnected", 12), ("ring", 11), ("ring", 12)},
            )

            job_id = rows[0]["status"] and store.conn.execute(
                "SELECT j.id FROM jobs j ORDER BY j.id LIMIT 1"
            ).fetchone()["id"]
            store.update_job_status(job_id, "DONE")
            plans.add_topology_benchmark_sqlite(
                db_path=db_path,
                benchmark_name="topology_smoke",
                swing_by_bodies=[["Kerbin"], ["Eve"], ["Moho"]],
                topologies=("unconnected", "ring"),
                seeds=[11, 12],
                t0_min=0,
                t0_bin=100,
                n_t0_bins=1,
                auto_tof=True,
                opt_level="debug",
                opt_injection="vinf",
            )
            self.assertEqual(
                store.conn.execute("SELECT status FROM jobs WHERE id = ?", (job_id,)).fetchone()["status"],
                "DONE",
            )
        finally:
            store.close()
            db_path.unlink(missing_ok=True)

    def test_run_metadata_records_benchmark_runtime_fields(self):
        db_path = TESTS_DIR / f"wayfinder_sqlite_test_{uuid.uuid4().hex}.sqlite"
        plans = Wayfinder(planet_pack="Vanilla")
        plans.add_batch_sqlite(
            db_path=db_path,
            batch_name="run_metadata",
            swing_by_bodies=[["Kerbin"], ["Eve"], ["Moho"]],
            t0_min=0,
            t0_bin=100,
            n_t0_bins=1,
            auto_tof=True,
            opt_level="debug",
            opt_injection="vinf",
            purpose="benchmark",
            optimizer_topology="ring",
            optimizer_seed=42,
        )

        store = SQLiteJobStore(db_path)
        try:
            job = store.pending_jobs("Vanilla", limit=1, batch_name="run_metadata")[0]
            run_id = store.start_run(
                job["job_id"],
                versions={"wayfinder": "test"},
                optimizer_metadata={
                    "optimizer_topology": job["optimizer_topology"],
                    "optimizer_seed": job["optimizer_seed"],
                    "requested_n_island": job["n_island"],
                    "actual_n_island": job["n_island"],
                    "island_pop": job["island_pop"],
                    "sade_gen": job["sade_gen"],
                    "n_evo_steps": job["n_evo_steps"],
                },
            )
            store.finish_run(run_id, runtime_seconds=12.5)
            row = store.conn.execute(
                """
                SELECT optimizer_topology, optimizer_seed, requested_n_island,
                       actual_n_island, island_pop, sade_gen, n_evo_steps,
                       runtime_seconds
                FROM runs WHERE id = ?
                """,
                (run_id,),
            ).fetchone()
            self.assertEqual(row["optimizer_topology"], "ring")
            self.assertEqual(row["optimizer_seed"], 42)
            self.assertEqual(row["runtime_seconds"], 12.5)
        finally:
            store.close()
            db_path.unlink(missing_ok=True)

    def test_query_best_known_route_from_sqlite_reference(self):
        db_path = TESTS_DIR / f"wayfinder_sqlite_test_{uuid.uuid4().hex}.sqlite"
        shutil.copyfile(REFERENCE_DB, db_path)

        store = SQLiteJobStore(db_path)
        try:
            best_route = store.best_results(
                planet_pack="JNSQ",
                batch_name="JNSQ_KEEMo_5yr_high_vinf_pykep3",
                start_body="Kerbin",
                target_body="Moho",
                t0_range=[0, 365 * 5],
                limit=1,
            )
            self.assertEqual(len(best_route), 1)
            self.assertEqual(best_route[0]["sequence_short_name"], "KEEMo")
            self.assertAlmostEqual(best_route[0]["objective_dv"], 4877.435209, places=3)

            best_keemo = store.best_results(
                planet_pack="JNSQ",
                batch_name="JNSQ_KEEMo_5yr_high_vinf_pykep3",
                sequence_short_name="KEEMo",
                contains_flyby="Eve",
                limit=1,
            )
            self.assertEqual(best_keemo[0]["body_path"], "Kerbin>Eve>Eve>Moho")
        finally:
            store.close()
            db_path.unlink(missing_ok=True)

    def test_find_best_known_plan_sqlite_displays_transx(self):
        db_path = TESTS_DIR / f"wayfinder_sqlite_test_{uuid.uuid4().hex}.sqlite"
        shutil.copyfile(REFERENCE_DB, db_path)
        try:
            plans = Wayfinder(planet_pack="JNSQ")
            output = io.StringIO()
            with contextlib.redirect_stdout(output):
                best = plans.find_best_known_plan_sqlite(
                    db_path=db_path,
                    batch_name="JNSQ_KEEMo_5yr_high_vinf_pykep3",
                    start_body="Kerbin",
                    target_body="Moho",
                    t0_range=[0, 365 * 5],
                )
        finally:
            db_path.unlink(missing_ok=True)

        self.assertEqual(best["sequence_short_name"], "KEEMo")
        self.assertAlmostEqual(best["objective_dv"], 4877.435209, places=3)
        self.assertIn("Arrival at Moho", output.getvalue())

    def test_plot_dv_vs_t0_sqlite_returns_best_points_per_bin(self):
        db_path = TESTS_DIR / f"wayfinder_sqlite_test_{uuid.uuid4().hex}.sqlite"
        shutil.copyfile(REFERENCE_DB, db_path)
        try:
            plans = Wayfinder(planet_pack="JNSQ")
            plotted = plans.plot_DVvsT0_sqlite(
                db_path=db_path,
                batch_name="JNSQ_KEEMo_5yr_high_vinf_pykep3",
                start_body="Kerbin",
                target_body="Moho",
                t0_range=[0, 365 * 5],
            )
        finally:
            db_path.unlink(missing_ok=True)

        self.assertEqual(len(plotted), 5)
        self.assertAlmostEqual(plotted["objective_dv"].min(), 4877.435209, places=3)

    def test_plot_by_sequences_sqlite_returns_best_per_sequence(self):
        db_path = TESTS_DIR / f"wayfinder_sqlite_test_{uuid.uuid4().hex}.sqlite"
        shutil.copyfile(REFERENCE_DB, db_path)
        try:
            plans = Wayfinder(planet_pack="JNSQ")
            plotted = plans.plot_by_sequences_sqlite(
                db_path=db_path,
                batch_name="JNSQ_KEEMo_5yr_high_vinf_pykep3",
                start_body="Kerbin",
                target_body="Moho",
                t0_range=[0, 365 * 5],
            )
        finally:
            db_path.unlink(missing_ok=True)

        self.assertEqual(len(plotted), 1)
        self.assertEqual(plotted.iloc[0]["sequence_short_name"], "KEEMo")
        self.assertAlmostEqual(plotted.iloc[0]["objective_dv"], 4877.435209, places=3)

    def test_porkchop_samples_store_metadata_and_plot_styles(self):
        db_path = TESTS_DIR / f"wayfinder_sqlite_test_{uuid.uuid4().hex}.sqlite"
        shutil.copyfile(REFERENCE_DB, db_path)
        cells_path = TESTS_DIR / f"wayfinder_sqlite_test_{uuid.uuid4().hex}_cells.png"
        continuous_path = TESTS_DIR / f"wayfinder_sqlite_test_{uuid.uuid4().hex}_continuous.png"

        store = SQLiteJobStore(db_path)
        try:
            best = store.best_results(
                planet_pack="JNSQ",
                batch_name="JNSQ_KEEMo_5yr_high_vinf_pykep3",
                sequence_short_name="KEEMo",
                limit=1,
            )[0]
            run_id = best["run_id"]
            rows = [
                {"t0": 100.1, "tof": 200.1, "metric": 5000.0, "fitness": 5000.0, "gene": [1.0]},
                {"t0": 101.1, "tof": 200.1, "metric": 5200.0, "fitness": 5200.0, "gene": [2.0]},
                {"t0": 100.1, "tof": 201.1, "metric": 5400.0, "fitness": 5400.0, "gene": [3.0]},
                {"t0": 101.1, "tof": 201.1, "metric": 5800.0, "fitness": 5800.0, "gene": [4.0]},
            ]
            store.replace_porkchop_samples(
                run_id,
                "unit_sampler",
                rows,
                metadata={"sampler_type": "unit", "sample_count": len(rows)},
            )
            metadata = store.porkchop_sampler_metadata(run_id, "unit_sampler")
        finally:
            store.close()

        try:
            self.assertEqual(metadata["sampler_type"], "unit")
            self.assertEqual(metadata["sample_count"], 4)

            plans = Wayfinder(planet_pack="JNSQ")
            cells = plans.plot_adaptive_binned_sampled_porkchop_sqlite(
                db_path,
                run_id,
                sampler_name="unit_sampler",
                coarse_bin_days=1.0,
                min_bin_days=1.0,
                style="cells",
                color_levels="low_detail",
                vmax="double_floor",
                output_path=cells_path,
            )
            continuous = plans.plot_adaptive_binned_sampled_porkchop_sqlite(
                db_path,
                run_id,
                sampler_name="unit_sampler",
                coarse_bin_days=1.0,
                min_bin_days=1.0,
                style="continuous",
                color_levels="low_detail",
                vmax="double_floor",
                output_path=continuous_path,
            )
            self.assertEqual(len(cells), 4)
            self.assertEqual(len(continuous), 4)
            self.assertTrue(cells_path.exists())
            self.assertTrue(continuous_path.exists())
        finally:
            plt.close("all")
            db_path.unlink(missing_ok=True)
            cells_path.unlink(missing_ok=True)
            continuous_path.unlink(missing_ok=True)


if __name__ == "__main__":
    unittest.main(verbosity=2)
