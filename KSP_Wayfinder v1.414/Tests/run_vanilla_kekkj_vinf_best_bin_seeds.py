# -*- coding: utf-8 -*-
"""Re-run the best coarse Vanilla KEKKJ vinf bin with several seeds."""

import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "WayfinderCore"))

from _Wayfinder import Wayfinder  # noqa: E402


DB_PATH = ROOT / "Tests" / "vanilla_kekkj_vinf_best_bin_seeds.sqlite"
BATCH_NAME = "vanilla_kekkj_vinf_best_bin_fully_connected_moderate"
SEEDS = (0, 1, 2, 3, 4)


def main():
    plans = Wayfinder(planet_pack="Vanilla")
    plans.add_topology_benchmark_sqlite(
        db_path=DB_PATH,
        benchmark_name=BATCH_NAME,
        swing_by_bodies=[["Kerbin"], ["Eve"], ["Kerbin"], ["Kerbin"], ["Jool"]],
        topologies=("fully_connected",),
        seeds=SEEDS,
        t0_min=400,
        t0_bin=200,
        n_t0_bins=1,
        tof_min=3000,
        tof_bin=500,
        n_tof_bins=1,
        auto_tof=False,
        opt_level="moderate",
        opt_injection="vinf",
        injection_altitude=100000,
    )
    optimized = 0
    for seed in SEEDS:
        batch_name = f"{BATCH_NAME}_fully_connected_seed{seed}"
        optimized += plans.optimize_sqlite(DB_PATH, n=1, batch_name=batch_name)
    print(f"Optimized {optimized} KEKKJ best-bin seed jobs")
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
    ]
    print(results.sort_values("objective_dv")[columns].to_string(index=False))


if __name__ == "__main__":
    main()
