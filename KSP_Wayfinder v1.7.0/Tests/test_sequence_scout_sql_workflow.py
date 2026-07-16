"""Tests for the persisted sequence-scout L0/funnel workflow."""

import sys
import unittest
import uuid
from pathlib import Path
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
CORE_DIR = ROOT / "WayfinderCore"
TESTS_DIR = ROOT / "Tests"
sys.path.insert(0, str(CORE_DIR))

from _LambertArcFilter import (  # noqa: E402
    LambertArcSolution,
    T0BinCandidate,
    T0BinScanResult,
)
from _OptimizationService import OptimizationService  # noqa: E402
from _SequenceScout import SequenceCandidate  # noqa: E402
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


class SequenceScoutSQLWorkflowTests(unittest.TestCase):
    def test_sql_workflow_strategies_are_explicit_continuations(self):
        l0 = OptimizationService.funnel_run_config(
            "funnel_l0_screen", 64, 8, 5, 20,
        )
        continuation = OptimizationService.funnel_run_config(
            "funnel_l0_continuation", 16, 32, 20, 20,
        )
        self.assertEqual(
            [stage["name"] for stage in l0["stages"]],
            ["scout_unconnected"],
        )
        self.assertEqual(
            [stage["name"] for stage in continuation["stages"]],
            ["wide", "mbh_refine", "intermediate", "exact_ejection"],
        )

    def test_prepare_promote_and_recover_parent_population(self):
        db_path = TESTS_DIR / (
            "wayfinder_scout_workflow_" + uuid.uuid4().hex + ".sqlite"
        )
        plans = Wayfinder(planet_pack="Vanilla")
        candidate = SequenceCandidate(
            sequence=("Kerbin", "Eve", "Jool"),
            proxy_cost_mps=1000.0,
            departure_vinf_mps=1000.0,
            arrival_vinf_mps=1000.0,
            max_turn_fraction=0.5,
            terminal_periapsis_m=1.0,
            terminal_apoapsis_m=2.0,
        )
        direct = LambertArcSolution(
            encounter_body="Jool", t0_days=50.0, tof_days=1000.0,
            departure_vinf_mps=3000.0, arrival_vinf_mps=2000.0,
            ejection_dv_mps=2000.0, ejection_strategy="direct",
        )
        first_arc = LambertArcSolution(
            encounter_body="Eve", t0_days=40.0, tof_days=160.0,
            departure_vinf_mps=1100.0, arrival_vinf_mps=900.0,
            ejection_dv_mps=1000.0, ejection_strategy="direct",
        )
        scanned = T0BinCandidate(
            candidate=candidate,
            bin_start_days=0.0,
            bin_end_days=100.0,
            solution=first_arc,
            direct_reference=direct,
            ejection_ratio=0.5,
            tisserand_vinf_relative_error=0.05,
            score=0.525,
        )
        scan = T0BinScanResult(
            t0_start_days=0.0,
            t0_end_days=100.0,
            bin_width_days=100.0,
            direct_reference=direct,
            candidates=(scanned,),
            maximum_ejection_ratio=1.05,
        )
        try:
            with patch.object(
                plans, "scout_sequences", return_value=[candidate],
            ), patch.object(
                plans, "scan_scout_sequence_bins", return_value=scan,
            ):
                prepared = plans.prepare_sequence_scout_sqlite(
                    db_path, "mock_scout", "Kerbin", "Jool",
                    config={
                        "t0_min_days": 0.0,
                        "t0_max_days": 100.0,
                        "bin_width_days": 100.0,
                        "top_sequences": 1,
                    },
                )
            self.assertEqual(prepared["l0_job_count"], 1)
            store = SQLiteJobStore(db_path)
            try:
                l0_job = store.pending_jobs(
                    "Vanilla", batch_name="mock_scout__l0", limit=1,
                )[0]
                self.assertEqual(
                    l0_job["optimizer_strategy"], "funnel_l0_screen",
                )
                run_id = store.upsert_result(
                    l0_job["job_id"], 1500.0, 20.0, 900.0, 1200.0,
                    [1.0, 2.0],
                )
                store.record_optimizer_population(
                    run_id, 5,
                    [FakePopulation(
                        [[1.0, 2.0], [3.0, 4.0]], [[1500.0], [1700.0]],
                    )],
                    source="stage_1_final",
                )
            finally:
                store.close()

            with self.assertRaisesRegex(RuntimeError, "still pending"):
                plans.promote_sequence_scout_sqlite(
                    db_path, prepared["scout_run_id"], per_bin=1,
                )
            promoted = plans.promote_sequence_scout_sqlite(
                db_path, prepared["scout_run_id"], per_bin=1,
                allow_partial=True,
            )
            self.assertEqual(promoted["promoted_job_count"], 1)
            self.assertEqual(promoted["newly_promoted_job_count"], 1)
            store = SQLiteJobStore(db_path)
            try:
                funnel_job = store.pending_jobs(
                    "Vanilla", batch_name="mock_scout__funnel", limit=1,
                )[0]
                self.assertEqual(
                    funnel_job["optimizer_strategy"],
                    "funnel_l0_continuation",
                )
                self.assertEqual(
                    store.promoted_seed_genes(funnel_job["job_id"]),
                    [[1.0, 2.0], [3.0, 4.0]],
                )
                lineage = store.conn.execute(
                    """
                    SELECT parent_job_id, parent_run_id, promotion_rank
                    FROM sequence_scout_jobs
                    WHERE job_id = ? AND workflow_stage = 'funnel'
                    """,
                    (funnel_job["job_id"],),
                ).fetchone()
                self.assertEqual(lineage["parent_job_id"], l0_job["job_id"])
                self.assertEqual(lineage["parent_run_id"], run_id)
                self.assertEqual(lineage["promotion_rank"], 1)
            finally:
                store.close()
        finally:
            db_path.unlink(missing_ok=True)


if __name__ == "__main__":
    unittest.main()
