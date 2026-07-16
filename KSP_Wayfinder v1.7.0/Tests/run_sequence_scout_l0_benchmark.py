# -*- coding: utf-8 -*-
"""Run the real Wayfinder L0 on promoted Tisserand/Lambert bin candidates."""

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pygmo as pg


ROOT = Path(__file__).resolve().parents[1]
CORE_DIR = ROOT / "WayfinderCore"
sys.path.insert(0, str(CORE_DIR))

from _OptimizationService import OptimizationService  # noqa: E402
from _Wayfinder import Wayfinder, _make_mga_1dsm  # noqa: E402


class MemoryTelemetryStore:
    """Minimal persistence boundary consumed by the funnel stage runner."""

    def __init__(self):
        self.snapshots = []
        self.stages = []

    def record_optimizer_snapshot(
        self, run_id, step, champions_f, champions_x, populations=None,
        telemetry=None,
    ):
        values = [float(value[0]) for value in champions_f]
        self.snapshots.append({
            "step": int(step),
            "best": min(values),
            "average_champion": float(np.mean(values)),
        })

    def record_optimizer_population(self, *args, **kwargs):
        return None

    def record_optimizer_stage(self, run_id, summary):
        self.stages.append(dict(summary))


def _sequence_name(sequence, abbreviations):
    return "".join(abbreviations.get(body, body) for body in sequence)


def _load_checkpoint(path):
    if path.exists():
        return json.loads(path.read_text(encoding="utf-8"))
    return {"rows": []}


def _save_checkpoint(path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8",
    )
    temporary.replace(path)


def _first_leg_inclusive_bounds(plans, sequence, lambert_tof_days):
    bounds = [
        list(pair) for pair in plans.estimate_direct_tof_bounds(
            sequence, profile="relaxed",
        )["direct_bounds_days"]
    ]
    hint = float(lambert_tof_days)
    bounds[0][0] = min(bounds[0][0], 0.8 * hint)
    bounds[0][1] = max(bounds[0][1], 1.2 * hint)
    return bounds


def _l0_stage_plan():
    plan = OptimizationService.funnel_stage_plan(
        16, 32, 20, 20, exact_strategy="scout_archive_nm_64",
    )
    stage = dict(plan[0])
    if stage["name"] != "scout_unconnected":
        raise RuntimeError("Expected scout_unconnected as the L0 stage")
    return [stage]


def _run_l0(plans, promoted, seed):
    sequence = tuple(promoted.candidate.sequence)
    bounds = _first_leg_inclusive_bounds(
        plans, sequence, promoted.solution.tof_days,
    )
    edy_to_kdy = float(plans._Edy2Kdy)
    target = plans._fullname_dic[sequence[-1]]
    udp = _make_mga_1dsm(
        seq=[plans._fullname_dic[name] for name in sequence],
        t0=[
            float(promoted.bin_start_days) / edy_to_kdy,
            float(promoted.bin_end_days) / edy_to_kdy,
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
        "job_id": int(seed),
        "ejection_altitude": 100000.0,
        "tof_encoding": "direct",
        "bodies_json": json.dumps(sequence),
        "leg_tof_bounds_json": json.dumps(bounds),
    }
    telemetry = MemoryTelemetryStore()
    soi_radius_by_name = {
        name: plans.soi_radius(body)
        for name, body in plans._fullname_dic.items()
    }
    pg.set_global_rng_seed(int(seed))
    started = time.perf_counter()
    result = plans._run_sqlite_funnel_job(
        telemetry,
        job,
        udp,
        sqlite_run_id=0,
        sade_gen=20,
        requested_n_island=64,
        n_island=64,
        island_pop=8,
        n_evo_steps=5,
        topology="unconnected",
        versions={},
        adaptive_stop=None,
        soi_radius_by_name=soi_radius_by_name,
        exact_strategy="scout_archive_nm_64",
        requested_optimizer_strategy="sequence_scout_l0_only",
        effective_optimizer_seed=int(seed),
        stage_plan_override=_l0_stage_plan(),
    )
    runtime = time.perf_counter() - started
    return result, telemetry, bounds, runtime


def _promoted_rows(plans, top_sequences=20):
    scout = plans.scout_sequences("Kerbin", "Jool")
    scan = plans.scan_scout_sequence_bins(
        scout,
        [0.0, 1000.0],
        bin_width_days=100.0,
        t0_step_days=10.0,
        tof_samples=40,
        maximum_ejection_ratio=1.05,
    )
    selected_sequences = {
        row.candidate.sequence
        for row in scan.ranked_unique(top_sequences)
    }
    promoted = [
        row for row in scan.candidates
        if row.candidate.sequence in selected_sequences
    ]
    return scan, sorted(
        promoted,
        key=lambda row: (
            row.bin_start_days, row.score, row.candidate.sequence,
        ),
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output", default="Tests/sequence_scout_l0_kj_0_1000.json",
    )
    parser.add_argument("--top-sequences", type=int, default=20)
    parser.add_argument("--max-jobs", type=int)
    parser.add_argument("--seed-base", type=int, default=42000)
    args = parser.parse_args()

    output = Path(args.output)
    plans = Wayfinder(planet_pack="Vanilla")
    scan_started = time.perf_counter()
    scan, promoted = _promoted_rows(plans, args.top_sequences)
    scan_runtime = time.perf_counter() - scan_started
    if args.max_jobs is not None:
        promoted = promoted[:max(0, int(args.max_jobs))]

    checkpoint = _load_checkpoint(output)
    completed = {
        (row["sequence"], float(row["bin_start_days"]))
        for row in checkpoint.get("rows", [])
    }
    checkpoint.update({
        "planet_pack": "Vanilla",
        "start": "Kerbin",
        "target": "Jool",
        "t0_bounds_days": [0.0, 1000.0],
        "bin_width_days": 100.0,
        "top_sequences": int(args.top_sequences),
        "promoted_jobs": len(promoted),
        "scan_runtime_seconds": scan_runtime,
        "direct_reference": scan.direct_reference.to_dict(),
        "l0_stage": _l0_stage_plan()[0],
    })
    checkpoint.setdefault("rows", [])
    _save_checkpoint(output, checkpoint)

    for index, promoted_row in enumerate(promoted):
        sequence = tuple(promoted_row.candidate.sequence)
        identity = ("-".join(sequence), float(promoted_row.bin_start_days))
        if identity in completed:
            continue
        seed = int(args.seed_base) + index
        label = _sequence_name(sequence, plans._Body_abrev_dic)
        print(
            "[{}/{}] {} T0 {:.0f}-{:.0f} seed {}".format(
                index + 1, len(promoted), label,
                promoted_row.bin_start_days,
                promoted_row.bin_end_days,
                seed,
            ),
            flush=True,
        )
        result, telemetry, bounds, runtime = _run_l0(
            plans, promoted_row, seed,
        )
        row = {
            "sequence": "-".join(sequence),
            "sequence_short_name": label,
            "bin_start_days": float(promoted_row.bin_start_days),
            "bin_end_days": float(promoted_row.bin_end_days),
            "scout_score": float(promoted_row.score),
            "scout_ejection_ratio": float(promoted_row.ejection_ratio),
            "lambert_t0_days": float(promoted_row.solution.t0_days),
            "lambert_tof_days": float(promoted_row.solution.tof_days),
            "leg_tof_bounds_days": bounds,
            "optimizer_seed": seed,
            "l0_objective_dv": float(result["result_DV"]),
            "l0_t0_days": float(result["result_t0"]),
            "l0_total_tof_days": float(result["result_tof"]),
            "l0_runtime_seconds": float(runtime),
            "l0_stage": result["stage_summaries"][0],
            "convergence": telemetry.snapshots,
        }
        checkpoint["rows"].append(row)
        _save_checkpoint(output, checkpoint)
        print(
            "  DV={:.1f} runtime={:.2f}s".format(
                row["l0_objective_dv"], runtime,
            ),
            flush=True,
        )

    rows = checkpoint["rows"]
    print(
        "Completed {} rows; median L0 DV {:.1f}; median runtime {:.2f}s".format(
            len(rows),
            float(np.median([row["l0_objective_dv"] for row in rows])),
            float(np.median([row["l0_runtime_seconds"] for row in rows])),
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
