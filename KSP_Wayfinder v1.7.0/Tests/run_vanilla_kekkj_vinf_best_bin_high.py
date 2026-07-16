# -*- coding: utf-8 -*-
"""Run the best coarse Vanilla KEKKJ vinf bin at high optimization level."""

import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "WayfinderCore"))

from _Wayfinder import Wayfinder  # noqa: E402


DB_PATH = ROOT / "Tests" / "vanilla_kekkj_vinf_best_bin_high.sqlite"
BATCH_NAME = "vanilla_kekkj_vinf_best_bin_fully_connected_high"


def main():
    plans = Wayfinder(planet_pack="Vanilla")
    plans.add_topology_benchmark_sqlite(
        db_path=DB_PATH,
        benchmark_name=BATCH_NAME,
        swing_by_bodies=[["Kerbin"], ["Eve"], ["Kerbin"], ["Kerbin"], ["Jool"]],
        topologies=("fully_connected",),
        seeds=(0,),
        t0_min=400,
        t0_bin=200,
        n_t0_bins=1,
        tof_min=3000,
        tof_bin=500,
        n_tof_bins=1,
        auto_tof=False,
        opt_level="high",
        opt_injection="vinf",
        injection_altitude=100000,
    )
    batch_name = f"{BATCH_NAME}_fully_connected_seed0"
    optimized = plans.optimize_sqlite(DB_PATH, n=1, batch_name=batch_name)
    print(f"Optimized {optimized} KEKKJ high job")
    results = plans.benchmark_results_sqlite(DB_PATH, benchmark_name=BATCH_NAME)
    columns = [
        "sequence_short_name",
        "optimizer_topology",
        "optimizer_seed",
        "t0_min",
        "t0_max",
        "tof_min",
        "tof_max",
        "objective_dv",
        "result_t0",
        "result_tof",
        "runtime_seconds",
        "actual_n_island",
    ]
    print(results.sort_values("objective_dv")[columns].to_string(index=False))


if __name__ == "__main__":
    main()
