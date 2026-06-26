# -*- coding: utf-8 -*-
"""SQLite datastore tests for cross-batch Wayfinder queries."""

import os
import sys
import contextlib
import hashlib
import io
import json
import shutil
import sqlite3
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
from _OptimizationService import OptimizationService  # noqa: E402
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
    def test_jobs_are_claimed_atomically_and_expired_leases_are_recovered(self):
        db_path = TESTS_DIR / f"wayfinder_sqlite_test_{uuid.uuid4().hex}.sqlite"
        plans = Wayfinder(planet_pack="Vanilla")
        plans.add_batch_sqlite(
            db_path=db_path,
            batch_name="atomic_claim",
            swing_by_bodies=[["Kerbin"], ["Eve"], ["Moho"]],
            t0_min=0,
            t0_bin=100,
            n_t0_bins=1,
            auto_tof=True,
            opt_level="debug",
            opt_injection="vinf",
        )

        first = SQLiteJobStore(db_path)
        second = SQLiteJobStore(db_path)
        try:
            claimed = first.claim_pending_jobs(
                "Vanilla", limit=1, batch_name="atomic_claim",
                worker_id="worker-a", lease_seconds=3600,
            )
            self.assertEqual(len(claimed), 1)
            self.assertEqual(claimed[0]["status"], "RUNNING")
            self.assertEqual(claimed[0]["worker_id"], "worker-a")
            self.assertIsNotNone(claimed[0]["claimed_at"])
            self.assertIsNotNone(claimed[0]["claim_expires_at"])
            self.assertEqual(
                second.claim_pending_jobs(
                    "Vanilla", limit=1, batch_name="atomic_claim",
                    worker_id="worker-b",
                ),
                [],
            )

            expired_run_id = first.start_run(claimed[0]["job_id"])

            first.conn.execute(
                "UPDATE jobs SET claim_expires_at = ? WHERE id = ?",
                ("2000-01-01T00:00:00Z", claimed[0]["job_id"]),
            )
            first.conn.commit()
            recovered = second.claim_pending_jobs(
                "Vanilla", limit=1, batch_name="atomic_claim",
                worker_id="worker-b", lease_seconds=3600,
            )
            self.assertEqual(len(recovered), 1)
            self.assertEqual(recovered[0]["job_id"], claimed[0]["job_id"])
            self.assertEqual(recovered[0]["worker_id"], "worker-b")
            expired_run = second.conn.execute(
                "SELECT status, ended_at, stop_reason, error "
                "FROM runs WHERE id = ?",
                (expired_run_id,),
            ).fetchone()
            self.assertEqual(expired_run["status"], "FAILED")
            self.assertIsNotNone(expired_run["ended_at"])
            self.assertEqual(expired_run["stop_reason"], "claim_expired")
            self.assertIn("claim expired", expired_run["error"].lower())
            self.assertFalse(first.renew_job_claim(
                claimed[0]["job_id"], "worker-a",
                claimed_at=claimed[0]["claimed_at"],
            ))
            self.assertTrue(second.renew_job_claim(
                claimed[0]["job_id"], "worker-b",
                claimed_at=recovered[0]["claimed_at"],
            ))
            recovered_run_id = second.start_run(recovered[0]["job_id"])
            with self.assertRaisesRegex(RuntimeError, "stale result"):
                first.complete_claimed_run(
                    claimed[0]["job_id"], expired_run_id, "worker-a",
                    claimed[0]["claimed_at"], 999.0, 1.0, 2.0, 3.0,
                    [1.0, 2.0],
                )
            still_owned = second.conn.execute(
                "SELECT status, worker_id FROM jobs WHERE id = ?",
                (claimed[0]["job_id"],),
            ).fetchone()
            self.assertEqual(still_owned["status"], "RUNNING")
            self.assertEqual(still_owned["worker_id"], "worker-b")
            self.assertEqual(second.count_rows("results"), 0)
            second.complete_claimed_run(
                recovered[0]["job_id"], recovered_run_id, "worker-b",
                recovered[0]["claimed_at"], 123.0, 10.0, 20.0, 30.0,
                [4.0, 5.0], arrival_vinf=40.0, runtime_seconds=5.0,
                actual_n_evo_steps=3, stop_reason="test_complete",
            )
            row = second.conn.execute(
                "SELECT status, claimed_at, claim_expires_at, worker_id "
                "FROM jobs WHERE id = ?",
                (claimed[0]["job_id"],),
            ).fetchone()
            self.assertEqual(row["status"], "DONE")
            self.assertIsNone(row["claimed_at"])
            self.assertIsNone(row["claim_expires_at"])
            self.assertIsNone(row["worker_id"])
            self.assertEqual(
                second.conn.execute(
                    "SELECT value FROM metadata WHERE key = 'schema_version'"
                ).fetchone()["value"],
                "15",
            )
        finally:
            first.close()
            second.close()
            db_path.unlink(missing_ok=True)

    def test_schema_13_is_migrated_to_schema_15(self):
        db_path = TESTS_DIR / f"wayfinder_sqlite_test_{uuid.uuid4().hex}.sqlite"
        store = SQLiteJobStore(db_path)
        store.close()
        conn = sqlite3.connect(db_path)
        try:
            conn.execute("DROP INDEX idx_jobs_claim")
            for column in ("claimed_at", "claim_expires_at", "worker_id"):
                conn.execute(f"ALTER TABLE jobs DROP COLUMN {column}")
            for column in (
                "effective_optimizer_seed", "funnel_config_json",
                "funnel_config_hash", "code_revision", "planet_pack_hash",
            ):
                conn.execute(f"ALTER TABLE runs DROP COLUMN {column}")
            for column in (
                "topology_name", "migration_rate", "exact_archive_size",
            ):
                conn.execute(
                    f"ALTER TABLE optimizer_stages DROP COLUMN {column}"
                )
            conn.execute(
                "UPDATE metadata SET value = '13' WHERE key = 'schema_version'"
            )
            conn.commit()
        finally:
            conn.close()

        migrated = SQLiteJobStore(db_path)
        try:
            job_columns = {
                row["name"]
                for row in migrated.conn.execute("PRAGMA table_info(jobs)")
            }
            stage_columns = {
                row["name"]
                for row in migrated.conn.execute(
                    "PRAGMA table_info(optimizer_stages)"
                )
            }
            self.assertTrue(
                {"claimed_at", "claim_expires_at", "worker_id"}
                <= job_columns
            )
            run_columns = {
                row["name"]
                for row in migrated.conn.execute("PRAGMA table_info(runs)")
            }
            self.assertTrue(
                {
                    "effective_optimizer_seed", "funnel_config_json",
                    "funnel_config_hash", "code_revision", "planet_pack_hash",
                }
                <= run_columns
            )
            self.assertTrue(
                {"topology_name", "migration_rate", "exact_archive_size"}
                <= stage_columns
            )
            self.assertEqual(
                migrated.conn.execute(
                    "SELECT value FROM metadata WHERE key = 'schema_version'"
                ).fetchone()["value"],
                "15",
            )
        finally:
            migrated.close()
            db_path.unlink(missing_ok=True)

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

    def test_optimize_sqlite_claims_only_the_job_it_is_starting(self):
        db_path = TESTS_DIR / f"wayfinder_sqlite_test_{uuid.uuid4().hex}.sqlite"
        plans = Wayfinder(planet_pack="Vanilla")
        plans.add_batch_sqlite(
            db_path=db_path,
            batch_name="one_claim_at_a_time",
            swing_by_bodies=[["Kerbin"], ["Eve"], ["Moho"]],
            t0_min=0,
            t0_bin=100,
            n_t0_bins=2,
            auto_tof=True,
            opt_level="debug",
            opt_injection="vinf",
        )
        running_counts = []

        def fake_optimize(store, job, versions, **kwargs):
            running_counts.append(store.conn.execute(
                "SELECT COUNT(*) AS n FROM jobs WHERE status = 'RUNNING'"
            ).fetchone()["n"])
            store.update_job_status(job["job_id"], "DONE")

        plans._optimize_sqlite_job = fake_optimize
        try:
            self.assertEqual(
                plans.optimize_sqlite(
                    db_path, n=2, batch_name="one_claim_at_a_time",
                    worker_id="test-worker",
                ),
                2,
            )
            self.assertEqual(running_counts, [1, 1])
        finally:
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
                populations=[
                    FakePopulation(
                        genes=[[1.0, 2.0], [3.0, 4.0]],
                        fitnesses=[[123.4], [456.6]],
                    )
                ],
                telemetry={
                    "island_count": 8,
                    "topology_vertices": 8,
                    "topology_edges": 16,
                    "migrants_published": 16,
                    "migration_islands_active": 8,
                    "migrations_accepted": 5,
                },
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
            convergence = store.optimizer_convergence(run_id)
            self.assertEqual(convergence[0]["best_fitness"], 123.4)
            self.assertAlmostEqual(convergence[0]["average_fitness"], 290.0)
            self.assertEqual(convergence[0]["island_count"], 8)
            self.assertEqual(convergence[0]["topology_vertices"], 8)
            self.assertEqual(convergence[0]["topology_edges"], 16)
            self.assertEqual(convergence[0]["migrants_published"], 16)
            self.assertEqual(convergence[0]["migrations_accepted"], 5)
            population_points = store.optimizer_population_points(run_id, source="final")
            self.assertEqual(population_points[1]["gene"], [3.0, 4.0])
            store.record_optimizer_stage(
                run_id,
                {
                    "stage_index": 1,
                    "stage_name": "wide",
                    "n_island": 8,
                    "island_pop": 32,
                    "sade_gen": 20,
                    "planned_evo_steps": 20,
                    "actual_evo_steps": 12,
                    "ejection_model": "approximate",
                    "topology": "ring",
                    "migration_rate": 2,
                    "exact_archive_size": 17,
                    "algorithms": ["sade", "simulated_annealing"],
                    "topology_vertices": 8,
                    "topology_edges": 16,
                    "migrants_published": 160,
                    "migration_islands_active": 80,
                    "migrations_accepted": 42,
                    "best_fitness": 123.4,
                    "average_fitness": 290.0,
                    "runtime_seconds": 4.5,
                    "stop_reason": "convergence_plateau",
                },
            )
            stages = store.optimizer_stages(run_id)
            self.assertEqual(stages[0]["stage_name"], "wide")
            self.assertEqual(stages[0]["algorithms"], ["sade", "simulated_annealing"])
            self.assertEqual(stages[0]["actual_evo_steps"], 12)
            self.assertEqual(stages[0]["topology_name"], "ring")
            self.assertEqual(stages[0]["migration_rate"], 2)
            self.assertEqual(stages[0]["exact_archive_size"], 17)
            self.assertEqual(stages[0]["topology_vertices"], 8)
            self.assertEqual(stages[0]["topology_edges"], 16)
            self.assertEqual(stages[0]["migrations_accepted"], 42)
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
            funnel_config = OptimizationService.funnel_run_config(
                "funnel_scout_archive_32",
                job["n_island"],
                job["island_pop"],
                job["n_evo_steps"],
                job["sade_gen"],
            )
            run_id = store.start_run(
                job["job_id"],
                versions={"wayfinder": "test"},
                optimizer_metadata={
                    "optimizer_topology": job["optimizer_topology"],
                    "optimizer_seed": job["optimizer_seed"],
                    "effective_optimizer_seed": 42,
                    "requested_n_island": job["n_island"],
                    "actual_n_island": job["n_island"],
                    "island_pop": job["island_pop"],
                    "sade_gen": job["sade_gen"],
                    "n_evo_steps": job["n_evo_steps"],
                    "optimizer_strategy": "funnel_scout_archive_32",
                    "funnel_config": funnel_config,
                    "code_revision": "test-revision",
                    "planet_pack_hash": "test-planet-pack",
                },
            )
            store.finish_run(
                run_id, runtime_seconds=12.5, actual_n_evo_steps=7,
                stop_reason="convergence_plateau",
            )
            row = store.conn.execute(
                """
                SELECT optimizer_topology, optimizer_seed, requested_n_island,
                       actual_n_island, island_pop, sade_gen, n_evo_steps,
                       effective_optimizer_seed, actual_n_evo_steps,
                       stop_reason, runtime_seconds, optimizer_strategy,
                       funnel_config_json, funnel_config_hash, code_revision,
                       planet_pack_hash
                FROM runs WHERE id = ?
                """,
                (run_id,),
            ).fetchone()
            self.assertEqual(row["optimizer_topology"], "ring")
            self.assertEqual(row["optimizer_seed"], 42)
            self.assertEqual(row["effective_optimizer_seed"], 42)
            self.assertEqual(row["optimizer_strategy"], "funnel_scout_archive_32")
            self.assertEqual(json.loads(row["funnel_config_json"]), funnel_config)
            self.assertEqual(
                row["funnel_config_hash"],
                hashlib.sha256(row["funnel_config_json"].encode()).hexdigest(),
            )
            self.assertEqual(row["code_revision"], "test-revision")
            self.assertEqual(row["planet_pack_hash"], "test-planet-pack")
            self.assertEqual(row["runtime_seconds"], 12.5)
            self.assertEqual(row["actual_n_evo_steps"], 7)
            self.assertEqual(row["stop_reason"], "convergence_plateau")
        finally:
            store.close()
            db_path.unlink(missing_ok=True)

    def test_historical_run_replay_preserves_exact_funnel_config_and_seed(self):
        db_path = TESTS_DIR / f"wayfinder_sqlite_test_{uuid.uuid4().hex}.sqlite"
        plans = Wayfinder(planet_pack="Vanilla")
        plans.add_batch_sqlite(
            db_path=db_path,
            batch_name="historical_replay",
            swing_by_bodies=[["Kerbin"], ["Eve"], ["Moho"]],
            t0_min=0,
            t0_bin=100,
            n_t0_bins=1,
            auto_tof=True,
            opt_level="debug",
            opt_injection="vinf",
            purpose="benchmark",
            optimizer_topology="ring",
            optimizer_seed=None,
        )

        store = SQLiteJobStore(db_path)
        try:
            job = store.pending_jobs("Vanilla", limit=1, batch_name="historical_replay")[0]
            effective_seed = 987654321
            funnel_config = OptimizationService.funnel_run_config(
                "funnel_scout_archive_64",
                job["n_island"],
                job["island_pop"],
                job["n_evo_steps"],
                job["sade_gen"],
            )
            first_run_id = store.start_run(
                job["job_id"],
                versions={"wayfinder": "test"},
                optimizer_metadata={
                    "optimizer_topology": job["optimizer_topology"],
                    "optimizer_seed": job["optimizer_seed"],
                    "effective_optimizer_seed": effective_seed,
                    "requested_n_island": job["n_island"],
                    "actual_n_island": job["n_island"],
                    "island_pop": job["island_pop"],
                    "sade_gen": job["sade_gen"],
                    "n_evo_steps": job["n_evo_steps"],
                    "optimizer_strategy": "funnel_scout_archive_64",
                    "funnel_config": funnel_config,
                    "code_revision": "test-revision",
                    "planet_pack_hash": "test-planet-pack",
                },
            )
            first_context = store.run_context(first_run_id)
            replay_config = json.loads(first_context["funnel_config_json"])

            self.assertEqual(first_context["effective_optimizer_seed"], effective_seed)
            self.assertEqual(replay_config, funnel_config)
            self.assertEqual(
                first_context["funnel_config_hash"],
                hashlib.sha256(
                    json.dumps(
                        replay_config, sort_keys=True, separators=(",", ":")
                    ).encode()
                ).hexdigest(),
            )

            replay_run_id = store.start_run(
                first_context["job_id"],
                versions={"wayfinder": "test"},
                optimizer_metadata={
                    "optimizer_topology": first_context["run_optimizer_topology"],
                    "optimizer_seed": first_context["run_optimizer_seed"],
                    "effective_optimizer_seed": first_context[
                        "effective_optimizer_seed"
                    ],
                    "requested_n_island": first_context["requested_n_island"],
                    "actual_n_island": first_context["actual_n_island"],
                    "island_pop": first_context["run_island_pop"],
                    "sade_gen": first_context["run_sade_gen"],
                    "n_evo_steps": first_context["run_n_evo_steps"],
                    "optimizer_strategy": first_context["optimizer_strategy"],
                    "funnel_config": replay_config,
                    "code_revision": first_context["code_revision"],
                    "planet_pack_hash": first_context["planet_pack_hash"],
                },
            )
            replay_context = store.run_context(replay_run_id)

            self.assertEqual(
                replay_context["effective_optimizer_seed"],
                first_context["effective_optimizer_seed"],
            )
            self.assertEqual(
                replay_context["funnel_config_hash"],
                first_context["funnel_config_hash"],
            )
            self.assertEqual(
                json.loads(replay_context["funnel_config_json"]),
                replay_config,
            )
            self.assertEqual(
                replay_context["planet_pack_hash"],
                first_context["planet_pack_hash"],
            )
            self.assertEqual(
                replay_context["code_revision"],
                first_context["code_revision"],
            )
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
