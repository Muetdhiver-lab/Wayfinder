# -*- coding: utf-8 -*-
"""Plot L1/L2/L3 convergence from stored optimizer telemetry.

This is intentionally an analysis script: it reads the benchmark SQLite DBs
created by ``run_l0_pressure_gate_benchmark.py`` and reconstructs the selected
cascade result without launching new optimizer jobs.
"""

from __future__ import annotations

import argparse
import re
import sqlite3
import statistics
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Tests"))

from benchmark_analysis import CASE_BATCH_NAMES, case_from_sequence  # noqa: E402

STAGE_LABELS = {
    1: "L0 scout",
    2: "L1 wide",
    3: "MBH bridge",
    4: "L2 intermediate",
    5: "L3 exact",
}

PLOT_STAGE_INDICES = (2, 4, 5)


def _connect(db_path):
    con = sqlite3.connect(db_path)
    con.row_factory = sqlite3.Row
    return con


def _run_rows(db_path):
    query = """
        select
            b.name as batch,
            s.short_name as sequence,
            j.optimizer_seed as seed,
            r.id as run_id,
            r.runtime_seconds as runtime_seconds,
            res.objective_dv as objective_dv
        from batches b
        join batch_jobs bj on bj.batch_id = b.id
        join jobs j on j.id = bj.job_id
        join sequences s on s.id = j.sequence_id
        join runs r on r.job_id = j.id
        join results res on res.run_id = r.id
        where r.status = 'DONE'
    """
    with _connect(db_path) as con:
        return [dict(row) for row in con.execute(query)]


def _case_from_sequence(sequence):
    return case_from_sequence(sequence)


def _index_runs(baseline_db, pressure_db, strategy):
    baseline = {}
    same_seed = {}
    retry_seed = {}
    for row in _run_rows(baseline_db):
        case = _case_from_sequence(row["sequence"])
        baseline[(case, int(row["seed"]))] = dict(row, db_path=str(baseline_db))

    for row in _run_rows(pressure_db):
        case = _case_from_sequence(row["sequence"])
        batch_case = CASE_BATCH_NAMES.get(case)
        if batch_case is None:
            continue
        pressure_re = (
            r"^l0pressure_" + re.escape(batch_case) + "_"
            + re.escape(strategy) + r"_seed(\d+)_"
        )
        retry_re = (
            r"^l0retry_" + re.escape(batch_case) + "_"
            + re.escape(strategy) + r"_seed(\d+)_plus"
        )
        pressure_match = re.search(pressure_re, row["batch"])
        if pressure_match:
            original_seed = int(pressure_match.group(1))
            same_seed[(case, original_seed)] = dict(
                row, db_path=str(pressure_db), original_seed=original_seed,
            )
            continue
        retry_match = re.search(retry_re, row["batch"])
        if retry_match:
            original_seed = int(retry_match.group(1))
            retry_seed[(case, original_seed)] = dict(
                row, db_path=str(pressure_db), original_seed=original_seed,
            )
    return baseline, same_seed, retry_seed


def _select_cascade_runs(baseline, same_seed, retry_seed, thresholds):
    selected = {}
    branch_kind = {}
    for key, base in sorted(baseline.items()):
        case, seed = key
        current = base
        kind = "baseline"
        same = same_seed.get(key)
        if same is not None:
            kind = "same_seed"
            if float(same["objective_dv"]) < float(current["objective_dv"]):
                current = same
            threshold = thresholds.get(case)
            still_suspect = (
                threshold is not None
                and float(current["objective_dv"]) > float(threshold)
            )
            retry = retry_seed.get(key)
            if still_suspect and retry is not None:
                kind = "cascade_retry"
                if float(retry["objective_dv"]) < float(current["objective_dv"]):
                    current = retry
        selected[key] = current
        branch_kind[key] = kind
    return selected, branch_kind


def _run_stages(db_path, run_id):
    query = """
        select *
        from optimizer_stages
        where run_id = ?
        order by stage_index
    """
    with _connect(db_path) as con:
        return [dict(row) for row in con.execute(query, (int(run_id),))]


def _run_snapshots(db_path, run_id):
    query = """
        select step, best_fitness, average_fitness
        from optimizer_snapshots
        where run_id = ?
        order by step
    """
    with _connect(db_path) as con:
        return [dict(row) for row in con.execute(query, (int(run_id),))]


def _stage_boundaries(stages):
    boundaries = []
    start = 1
    for stage in stages:
        count = int(stage["actual_evo_steps"])
        end = start + count - 1
        boundaries.append((int(stage["stage_index"]), start, end, stage))
        start = end + 1
    return boundaries


def _snapshot_records(run, mode, branch_kind):
    stages = _run_stages(run["db_path"], run["run_id"])
    snapshots = _run_snapshots(run["db_path"], run["run_id"])
    boundaries = _stage_boundaries(stages)
    records = []
    for snapshot in snapshots:
        step = int(snapshot["step"])
        stage_index = None
        local_step = None
        stage = None
        for candidate_index, start, end, candidate_stage in boundaries:
            if start <= step <= end:
                stage_index = candidate_index
                local_step = step - start + 1
                stage = candidate_stage
                break
        if stage_index is None:
            continue
        records.append({
            "mode": mode,
            "branch_kind": branch_kind,
            "case": _case_from_sequence(run["sequence"]),
            "seed": int(run.get("original_seed", run["seed"])),
            "run_id": int(run["run_id"]),
            "stage_index": stage_index,
            "stage_name": stage["stage_name"],
            "stage_label": STAGE_LABELS.get(stage_index, stage["stage_name"]),
            "local_step": local_step,
            "global_step": step,
            "best_fitness": float(snapshot["best_fitness"]),
            "average_fitness": (
                float(snapshot["average_fitness"])
                if snapshot["average_fitness"] is not None else None
            ),
        })
    return records


def _stage_speed_records(run, mode, branch_kind):
    records = []
    stages = _run_stages(run["db_path"], run["run_id"])
    snapshots = _snapshot_records(run, mode, branch_kind)
    by_stage = defaultdict(list)
    for snapshot in snapshots:
        by_stage[int(snapshot["stage_index"])].append(snapshot)
    for stage in stages:
        stage_index = int(stage["stage_index"])
        points = by_stage.get(stage_index, [])
        if not points:
            continue
        start_best = float(points[0]["best_fitness"])
        end_best = float(stage["best_fitness"])
        runtime = float(stage["runtime_seconds"])
        improvement = start_best - end_best
        records.append({
            "mode": mode,
            "branch_kind": branch_kind,
            "case": _case_from_sequence(run["sequence"]),
            "seed": int(run.get("original_seed", run["seed"])),
            "stage_index": stage_index,
            "stage_label": STAGE_LABELS.get(stage_index, stage["stage_name"]),
            "runtime_seconds": runtime,
            "improvement": improvement,
            "improvement_per_second": (
                improvement / runtime if runtime > 0.0 else float("nan")
            ),
            "start_best": start_best,
            "end_best": end_best,
        })
    return records


def _median(values):
    values = [float(value) for value in values]
    return statistics.median(values) if values else float("nan")


def _quantile(values, fraction):
    values = sorted(float(value) for value in values)
    if not values:
        return float("nan")
    index = int(round(float(fraction) * (len(values) - 1)))
    return values[max(0, min(len(values) - 1, index))]


def _mode_style(mode):
    if mode == "baseline":
        return "#4C78A8"
    return "#F58518"


def plot_convergence(records, output_path):
    cases = sorted(set(record["case"] for record in records))
    fig, axes = plt.subplots(
        len(cases),
        len(PLOT_STAGE_INDICES),
        figsize=(4.0 * len(PLOT_STAGE_INDICES), 3.0 * len(cases)),
        squeeze=False,
        sharey=False,
    )
    for row_index, case in enumerate(cases):
        for col_index, stage_index in enumerate(PLOT_STAGE_INDICES):
            ax = axes[row_index][col_index]
            stage_label = STAGE_LABELS[stage_index]
            for mode in ("baseline", "cascade"):
                subset = [
                    record for record in records
                    if record["case"] == case
                    and record["stage_index"] == stage_index
                    and record["mode"] == mode
                ]
                by_step = defaultdict(list)
                for record in subset:
                    by_step[int(record["local_step"])].append(record["best_fitness"])
                if not by_step:
                    continue
                steps = sorted(by_step)
                medians = [_median(by_step[step]) for step in steps]
                lows = [_quantile(by_step[step], 0.25) for step in steps]
                highs = [_quantile(by_step[step], 0.75) for step in steps]
                color = _mode_style(mode)
                ax.plot(
                    steps,
                    medians,
                    label=mode,
                    color=color,
                    linewidth=2.0,
                )
                ax.fill_between(steps, lows, highs, color=color, alpha=0.16)
            ax.set_title(f"{case} — {stage_label}")
            ax.set_xlabel("local evo step")
            if col_index == 0:
                ax.set_ylabel("best objective DV (m/s)")
            ax.grid(True, alpha=0.25)
            ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=150)
    plt.close(fig)
    return output_path


def plot_speed(speed_records, output_path):
    cases = sorted(set(record["case"] for record in speed_records))
    fig, axes = plt.subplots(
        len(cases), 1, figsize=(10.5, 3.6 * len(cases)), squeeze=False,
    )
    for row_index, case in enumerate(cases):
        ax = axes[row_index][0]
        labels = [STAGE_LABELS[index] for index in PLOT_STAGE_INDICES]
        x = list(range(len(labels)))
        width = 0.36
        for offset, mode in [(-width / 2, "baseline"), (width / 2, "cascade")]:
            medians = []
            for stage_index in PLOT_STAGE_INDICES:
                values = [
                    record["improvement_per_second"]
                    for record in speed_records
                    if record["case"] == case
                    and record["stage_index"] == stage_index
                    and record["mode"] == mode
                ]
                medians.append(_median(values))
            ax.bar(
                [value + offset for value in x],
                medians,
                width=width,
                label=mode,
                color=_mode_style(mode),
                alpha=0.72,
            )
        ax.axhline(0.0, color="black", linewidth=0.8, alpha=0.5)
        ax.set_title(f"{case} — median stage improvement speed")
        ax.set_ylabel("best DV improvement per second")
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=20, ha="right")
        ax.grid(True, axis="y", alpha=0.25)
        ax.legend(loc="best")
    fig.tight_layout()
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=150)
    plt.close(fig)
    return output_path


def _print_summary(selected, branch_kind, speed_records):
    print("Selected cascade runs:")
    for (case, seed), run in sorted(selected.items()):
        print(
            f"  {case} seed={seed}: {branch_kind[(case, seed)]} "
            f"DV={float(run['objective_dv']):.3f} "
            f"runtime={float(run['runtime_seconds']):.2f}s"
        )
    print("\nMedian stage improvement speed (m/s/s):")
    for case in sorted(set(record["case"] for record in speed_records)):
        for mode in ("baseline", "cascade"):
            chunks = []
            for stage_index in PLOT_STAGE_INDICES:
                values = [
                    record["improvement_per_second"]
                    for record in speed_records
                    if record["case"] == case
                    and record["mode"] == mode
                    and record["stage_index"] == stage_index
                ]
                chunks.append(
                    f"{STAGE_LABELS[stage_index]}={_median(values):.1f}"
                )
            print(f"  {case} {mode}: " + ", ".join(chunks))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--baseline-db",
        default=str(ROOT / "Tests" / "l1_combined_10seeds_funnel_l0_recall_64_mbh_between.sqlite"),
    )
    parser.add_argument(
        "--pressure-db",
        default=str(ROOT / "Tests" / "l1_combined_10seeds_funnel_l0_recall_64_mbh_between_pressure_gate.sqlite"),
    )
    parser.add_argument(
        "--strategy",
        default="funnel_l0_recall_64_mbh_between",
    )
    parser.add_argument(
        "--quality-thresholds",
        default="KEEMo:5500,KEKKJ:1500",
    )
    parser.add_argument(
        "--convergence-plot",
        default=str(ROOT / "Tests" / "l123_convergence_cascade_10seeds.png"),
    )
    parser.add_argument(
        "--speed-plot",
        default=str(ROOT / "Tests" / "l123_stage_speed_cascade_10seeds.png"),
    )
    args = parser.parse_args()

    thresholds = {}
    for item in args.quality_thresholds.split(","):
        if not item.strip():
            continue
        key, value = item.split(":", 1)
        thresholds[key.strip()] = float(value)

    baseline_db = Path(args.baseline_db)
    pressure_db = Path(args.pressure_db)
    baseline, same_seed, retry_seed = _index_runs(
        baseline_db, pressure_db, args.strategy,
    )
    selected, branch_kind = _select_cascade_runs(
        baseline, same_seed, retry_seed, thresholds,
    )

    snapshot_records = []
    speed_records = []
    for key, run in sorted(baseline.items()):
        snapshot_records.extend(_snapshot_records(run, "baseline", "baseline"))
        speed_records.extend(_stage_speed_records(run, "baseline", "baseline"))
    for key, run in sorted(selected.items()):
        snapshot_records.extend(
            _snapshot_records(run, "cascade", branch_kind[key])
        )
        speed_records.extend(
            _stage_speed_records(run, "cascade", branch_kind[key])
        )

    convergence_plot = plot_convergence(snapshot_records, args.convergence_plot)
    speed_plot = plot_speed(speed_records, args.speed_plot)
    _print_summary(selected, branch_kind, speed_records)
    print(f"\nSaved convergence plot: {convergence_plot}")
    print(f"Saved speed plot: {speed_plot}")


if __name__ == "__main__":
    main()
