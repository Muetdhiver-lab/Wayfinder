# -*- coding: utf-8 -*-
"""Run the first topology benchmark bunch.

This intentionally runs real optimizations and may take a while. The bunch uses
one known-good bin for each sequence so topology comparisons are not dominated
by bad launch-window selection.
"""

import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "WayfinderCore"))

from _Wayfinder import Wayfinder  # noqa: E402


DB_PATH = ROOT / "Tests" / "topology_benchmark_bunch1.sqlite"
TOPOLOGIES = ("unconnected", "ring", "fully_connected")
SEEDS = (0,)
OPT_LEVEL = "moderate"


CASES = [
    {
        "planet_pack": "JNSQ",
        "benchmark_name": "bench_jnsq_keemo_t1600_moderate",
        "swing_by_bodies": [["Kerbin"], ["Eve"], ["Eve"], ["Moho"]],
        "t0_min": 1460,
        "t0_bin": 365,
        "tof_min": 400,
        "tof_bin": 700,
        "n_t0_bins": 1,
        "n_tof_bins": 1,
        "opt_injection": "vinf",
        "injection_altitude": 100000,
    },
    {
        "planet_pack": "Vanilla",
        "benchmark_name": "bench_vanilla_kekkj_ref_moderate",
        "swing_by_bodies": [["Kerbin"], ["Eve"], ["Kerbin"], ["Kerbin"], ["Jool"]],
        "t0_min": 600,
        "t0_bin": 200,
        "tof_min": 3000,
        "tof_bin": 1000,
        "n_t0_bins": 1,
        "n_tof_bins": 1,
        "opt_injection": "vinf",
        "injection_altitude": 100000,
    },
]


def prepare_benchmark_jobs():
    created = []
    for case in CASES:
        plans = Wayfinder(planet_pack=case["planet_pack"])
        created.extend(
            plans.add_topology_benchmark_sqlite(
                db_path=DB_PATH,
                benchmark_name=case["benchmark_name"],
                swing_by_bodies=case["swing_by_bodies"],
                topologies=TOPOLOGIES,
                seeds=SEEDS,
                t0_min=case["t0_min"],
                t0_bin=case["t0_bin"],
                n_t0_bins=case["n_t0_bins"],
                tof_min=case["tof_min"],
                tof_bin=case["tof_bin"],
                n_tof_bins=case["n_tof_bins"],
                auto_tof=False,
                opt_level=OPT_LEVEL,
                opt_injection=case["opt_injection"],
                injection_altitude=case["injection_altitude"],
            )
        )
    return created


def run_pending_jobs():
    total = 0
    for case in CASES:
        plans = Wayfinder(planet_pack=case["planet_pack"])
        for topology in TOPOLOGIES:
            for seed in SEEDS:
                batch_name = f"{case['benchmark_name']}_{topology}_seed{seed}"
                total += plans.optimize_sqlite(DB_PATH, n=1, batch_name=batch_name)
    return total


def print_summary():
    frames = []
    for planet_pack in {case["planet_pack"] for case in CASES}:
        plans = Wayfinder(planet_pack=planet_pack)
        frame = plans.benchmark_results_sqlite(DB_PATH)
        if not frame.empty:
            frame.insert(0, "planet_pack", planet_pack)
            frames.append(frame)
    for frame in frames:
        columns = [
            "planet_pack",
            "batch_name",
            "sequence_short_name",
            "optimizer_topology",
            "optimizer_seed",
            "objective_dv",
            "runtime_seconds",
            "actual_n_island",
        ]
        print(frame[columns].to_string(index=False))


def main():
    created = prepare_benchmark_jobs()
    print(f"Prepared {len(created)} benchmark batches in {DB_PATH}")
    optimized = run_pending_jobs()
    print(f"Optimized {optimized} benchmark jobs")
    print_summary()


if __name__ == "__main__":
    main()
