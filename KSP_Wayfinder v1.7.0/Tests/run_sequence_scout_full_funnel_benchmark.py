# -*- coding: utf-8 -*-
"""Promote the best sequence-scout L0 result per bin to full funnels."""

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pygmo as pg


ROOT = Path(__file__).resolve().parents[1]
CORE_DIR = ROOT / "WayfinderCore"
TESTS_DIR = ROOT / "Tests"
sys.path.insert(0, str(CORE_DIR))
sys.path.insert(0, str(TESTS_DIR))

from _OptimizationService import OptimizationService  # noqa: E402
from _Wayfinder import Wayfinder, _make_mga_1dsm  # noqa: E402
from run_sequence_scout_l0_benchmark import MemoryTelemetryStore  # noqa: E402


EXACT_STRATEGY = "scout_archive_nm_64_mbh_between"


def _save(path, payload):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8",
    )
    temporary.replace(path)


def _select_l0_rows(payload, per_bin):
    selected = []
    bins = sorted({float(row["bin_start_days"]) for row in payload["rows"]})
    for bin_start in bins:
        rows = sorted(
            (
                row for row in payload["rows"]
                if float(row["bin_start_days"]) == bin_start
            ),
            key=lambda row: (
                float(row["l0_objective_dv"]), row["sequence"],
            ),
        )
        selected.extend(rows[:int(per_bin)])
    return selected


def _build_problem(plans, row):
    sequence = tuple(row["sequence"].split("-"))
    bounds = row["leg_tof_bounds_days"]
    edy_to_kdy = float(plans._Edy2Kdy)
    target = plans._fullname_dic[sequence[-1]]
    udp = _make_mga_1dsm(
        seq=[plans._fullname_dic[name] for name in sequence],
        t0=[
            float(row["bin_start_days"]) / edy_to_kdy,
            float(row["bin_end_days"]) / edy_to_kdy,
        ],
        tof=[
            [float(lower) / edy_to_kdy, float(upper) / edy_to_kdy]
            for lower, upper in bounds
        ],
        vinf=[0.8, 1.8],
        tof_encoding="direct",
        add_vinf_dep=True,
        add_vinf_arr=False,
        multi_objective=False,
        orbit_insertion=False,
        rp_target=plans.rp_target_ward(target, 100000.0),
        e_target=0.0,
        max_revs=0,
    )
    job = {
        "job_id": int(row["optimizer_seed"]),
        "ejection_altitude": 100000.0,
        "tof_encoding": "direct",
        "bodies_json": json.dumps(sequence),
        "leg_tof_bounds_json": json.dumps(bounds),
    }
    return udp, job


def _run_full(plans, row):
    udp, job = _build_problem(plans, row)
    seed = int(row["optimizer_seed"])
    telemetry = MemoryTelemetryStore()
    soi_radius_by_name = {
        name: plans.soi_radius(body)
        for name, body in plans._fullname_dic.items()
    }
    stage_plan = OptimizationService.funnel_stage_plan(
        16, 32, 20, 20, exact_strategy=EXACT_STRATEGY,
    )
    pg.set_global_rng_seed(seed)
    started = time.perf_counter()
    result = plans._run_sqlite_funnel_job(
        telemetry,
        job,
        udp,
        sqlite_run_id=0,
        sade_gen=20,
        requested_n_island=16,
        n_island=16,
        island_pop=32,
        n_evo_steps=20,
        topology="ring",
        versions={},
        adaptive_stop=None,
        soi_radius_by_name=soi_radius_by_name,
        exact_strategy=EXACT_STRATEGY,
        requested_optimizer_strategy="sequence_scout_full_funnel",
        effective_optimizer_seed=seed,
        stage_plan_override=stage_plan,
    )
    return result, telemetry, time.perf_counter() - started


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--l0-input", default="Tests/sequence_scout_l0_kj_0_1000.json",
    )
    parser.add_argument(
        "--output", default="Tests/sequence_scout_full_kj_0_1000.json",
    )
    parser.add_argument("--per-bin", type=int, default=2)
    parser.add_argument(
        "--max-jobs", type=int, default=None,
        help="Optional cap used for smoke tests; selection order is preserved.",
    )
    args = parser.parse_args()

    input_payload = json.loads(Path(args.l0_input).read_text(encoding="utf-8"))
    selected = _select_l0_rows(input_payload, args.per_bin)
    if args.max_jobs is not None:
        selected = selected[:max(0, int(args.max_jobs))]
    output = Path(args.output)
    if output.exists():
        payload = json.loads(output.read_text(encoding="utf-8"))
    else:
        payload = {
            "source": str(args.l0_input),
            "per_bin": int(args.per_bin),
            "selected_jobs": len(selected),
            "exact_strategy": EXACT_STRATEGY,
            "rows": [],
        }
    completed = {
        (row["sequence"], float(row["bin_start_days"]))
        for row in payload["rows"]
    }
    plans = Wayfinder(planet_pack="Vanilla")
    for index, source_row in enumerate(selected):
        identity = (
            source_row["sequence"], float(source_row["bin_start_days"]),
        )
        if identity in completed:
            continue
        print(
            "[{}/{}] {} T0 {:.0f}-{:.0f}, L0 {:.1f}".format(
                index + 1, len(selected), source_row["sequence_short_name"],
                source_row["bin_start_days"], source_row["bin_end_days"],
                source_row["l0_objective_dv"],
            ),
            flush=True,
        )
        result, telemetry, runtime = _run_full(plans, source_row)
        row = {
            "sequence": source_row["sequence"],
            "sequence_short_name": source_row["sequence_short_name"],
            "bin_start_days": float(source_row["bin_start_days"]),
            "bin_end_days": float(source_row["bin_end_days"]),
            "optimizer_seed": int(source_row["optimizer_seed"]),
            "screening_l0_dv": float(source_row["l0_objective_dv"]),
            "objective_dv": float(result["result_DV"]),
            "result_t0_days": float(result["result_t0"]),
            "result_tof_days": float(result["result_tof"]),
            "runtime_seconds": float(runtime),
            "gene": [float(value) for value in result["gene"]],
            "stage_summaries": result["stage_summaries"],
            "convergence": telemetry.snapshots,
        }
        payload["rows"].append(row)
        _save(output, payload)
        print(
            "  final DV={:.1f}, runtime={:.2f}s".format(
                row["objective_dv"], runtime,
            ),
            flush=True,
        )

    rows = payload["rows"]
    print(
        "Completed {} funnels; best {:.1f}; median {:.1f}; "
        "median runtime {:.2f}s".format(
            len(rows),
            min(float(row["objective_dv"]) for row in rows),
            float(np.median([row["objective_dv"] for row in rows])),
            float(np.median([row["runtime_seconds"] for row in rows])),
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
