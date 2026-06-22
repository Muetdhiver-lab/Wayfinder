# -*- coding: utf-8 -*-
"""Smoke-test alpha to direct local re-optimization on a known KEEMo run."""

import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "WayfinderCore"))

from _SQLiteStore import SQLiteJobStore  # noqa: E402
from _Wayfinder import Wayfinder  # noqa: E402


REFERENCE_DB = ROOT / "Tests" / "wayfinder_reference.sqlite"
BATCH_NAME = "JNSQ_KEEMo_5yr_high_vinf_pykep3"


def run_smoke_test():
    plans = Wayfinder(planet_pack="JNSQ")
    store = SQLiteJobStore(REFERENCE_DB)
    try:
        best = store.best_results(
            planet_pack="JNSQ",
            batch_name=BATCH_NAME,
            sequence_short_name="KEEMo",
            limit=1,
        )[0]
    finally:
        store.close()

    result = plans.reoptimize_direct_near_run_sqlite(
        REFERENCE_DB,
        run_id=best["run_id"],
        t0_wiggle_days=5,
        leg_tof_wiggle_days=5,
        sade_gen=40,
        pop_size=64,
    )
    print(json.dumps(result, indent=2))
    assert abs(result["seed_objective"] - best["objective_dv"]) < 1e-6
    assert result["champion_objective"] < 6000

    cloud = plans.sample_direct_local_cloud_sqlite(
        REFERENCE_DB,
        run_id=best["run_id"],
        sampler_name="direct_reopt_smoke_cloud",
        t0_wiggle_days=2,
        leg_tof_wiggle_days=2,
        sade_gen=10,
        n_island=3,
        island_pop=12,
    )
    valid = cloud.dropna(subset=["result_DV", "result_t0", "result_tof"])
    assert len(valid) > 1
    assert valid["result_DV"].min() < 6000


if __name__ == "__main__":
    run_smoke_test()
