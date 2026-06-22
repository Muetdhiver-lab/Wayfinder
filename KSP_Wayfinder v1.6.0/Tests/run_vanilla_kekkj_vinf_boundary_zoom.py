# -*- coding: utf-8 -*-
"""Zoom around the Vanilla KEKKJ vinf bin boundary found in the coarse scan."""

import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "WayfinderCore"))

from _Wayfinder import Wayfinder  # noqa: E402


DB_PATH = ROOT / "Tests" / "vanilla_kekkj_vinf_boundary_zoom.sqlite"
BATCH_NAME = "vanilla_kekkj_vinf_boundary_zoom_moderate"


def main():
    plans = Wayfinder(planet_pack="Vanilla")
    plans.add_batch_sqlite(
        db_path=DB_PATH,
        batch_name=BATCH_NAME,
        swing_by_bodies=[["Kerbin"], ["Eve"], ["Kerbin"], ["Kerbin"], ["Jool"]],
        t0_min=550,
        t0_bin=150,
        n_t0_bins=2,
        tof_min=3250,
        tof_bin=300,
        n_tof_bins=2,
        auto_tof=False,
        opt_level="moderate",
        opt_injection="vinf",
        injection_altitude=100000,
        purpose="benchmark",
        optimizer_topology="fully_connected",
        optimizer_seed=0,
    )
    optimized = plans.optimize_sqlite(DB_PATH, n=4, batch_name=BATCH_NAME)
    print(f"Optimized {optimized} KEKKJ boundary zoom jobs")
    results = plans.benchmark_results_sqlite(DB_PATH, benchmark_name=BATCH_NAME)
    columns = [
        "sequence_short_name",
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
