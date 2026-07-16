# -*- coding: utf-8 -*-
"""Retrospective threshold-free adaptive handoff benchmark.

The policy never sees final Delta-V when making a decision. It combines two
dimensionless signals available at the end of L2:

* relative MBH-to-L2 improvement;
* relative gap between the best and median exact-archive candidates.

Cut points are quantiles learned from the other runs (leave-one-out), rather
than sequence-specific objective thresholds. A triggered Hill-Valley branch is
kept only when its final exact result beats the current branch.
"""

from __future__ import annotations

import argparse
import json
import sqlite3
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "WayfinderCore"))

from _SQLiteStore import SQLiteJobStore  # noqa: E402


def _run_rows(db_path):
    query = """
        SELECT r.id AS run_id, s.short_name AS case_name,
               j.optimizer_seed AS seed, res.objective_dv AS dv,
               r.runtime_seconds AS runtime_seconds
        FROM runs r
        JOIN jobs j ON j.id = r.job_id
        JOIN sequences s ON s.id = j.sequence_id
        JOIN results res ON res.run_id = r.id
        ORDER BY s.short_name, j.optimizer_seed, r.id DESC
    """
    with sqlite3.connect(db_path) as connection:
        connection.row_factory = sqlite3.Row
        rows = connection.execute(query).fetchall()
    latest = {}
    for row in rows:
        key = (row["case_name"], int(row["seed"]))
        latest.setdefault(key, dict(row))
    return latest


def _archive_gap(store, run_id):
    points = [
        point
        for point in store.optimizer_population_points(
            run_id, source="exact_archive_stage3_seed",
        )
        if int(point["step"]) < 100
    ]
    if not points:
        raise ValueError("Run {} has no baseline exact archive".format(run_id))
    fitness = sorted(float(point["fitness"]) for point in points)
    median = float(np.median(fitness))
    return (median - fitness[0]) / max(abs(median), 1.0)


def _l2_signals(store, run_id):
    stages = [
        stage for stage in store.optimizer_stages(run_id)
        if int(stage["stage_index"]) < 100
    ]
    names = [stage["stage_name"].split("::")[-1] for stage in stages]
    mbh_index = names.index("mbh_refine")
    l2_index = names.index("intermediate")
    mbh_best = float(stages[mbh_index]["best_fitness"])
    l2_best = float(stages[l2_index]["best_fitness"])
    l2_improvement = (mbh_best - l2_best) / max(abs(mbh_best), 1.0)

    convergence = [
        point for point in store.optimizer_convergence(run_id)
        if int(point["step"]) < 100
    ]
    l2_steps = int(stages[l2_index]["actual_evo_steps"])
    l2_end = sum(int(stage["actual_evo_steps"]) for stage in stages[:l2_index + 1])
    l2_history = convergence[l2_end - l2_steps:l2_end]
    midpoint = max(0, len(l2_history) // 2 - 1)
    late_best_improvement = (
        float(l2_history[midpoint]["best_fitness"])
        - float(l2_history[-1]["best_fitness"])
    ) / max(abs(float(l2_history[midpoint]["best_fitness"])), 1.0)
    return {
        "l2_improvement": float(l2_improvement),
        "l2_late_improvement": float(late_best_improvement),
        "archive_gap": float(_archive_gap(store, run_id)),
    }


def _downstream_runtimes(equal_budget_json):
    payload = json.loads(Path(equal_budget_json).read_text(encoding="utf-8"))
    return {
        (row["case"], int(row["seed"])): float(row["total_runtime"])
        for row in payload["rows"]
        if row["branch"] == "hill_valley"
    }


def extract_rows(current_db, hill_valley_db, equal_budget_json=None):
    current = _run_rows(current_db)
    hill_valley = _run_rows(hill_valley_db)
    runtimes = (
        _downstream_runtimes(equal_budget_json)
        if equal_budget_json and Path(equal_budget_json).exists() else {}
    )
    if set(current) != set(hill_valley):
        raise ValueError("Current and Hill-Valley DBs do not contain the same runs")
    store = SQLiteJobStore(current_db)
    rows = []
    try:
        for key in sorted(current):
            base = current[key]
            alternative = hill_valley[key]
            signals = _l2_signals(store, base["run_id"])
            rows.append({
                "case": key[0],
                "seed": key[1],
                "current_dv": float(base["dv"]),
                "hill_valley_dv": float(alternative["dv"]),
                "current_runtime": float(base["runtime_seconds"]),
                "hill_valley_downstream_runtime": runtimes.get(
                    key, float(alternative["runtime_seconds"]),
                ),
                **signals,
            })
    finally:
        store.close()
    return rows


def quantile_decision(
    row, calibration_rows, progress_quantile=0.75,
    archive_low_quantile=0.30, archive_high_quantile=0.70,
):
    progress_cut = float(np.quantile(
        [item["l2_improvement"] for item in calibration_rows],
        progress_quantile,
    ))
    archive_values = [item["archive_gap"] for item in calibration_rows]
    archive_low = float(np.quantile(archive_values, archive_low_quantile))
    archive_high = float(np.quantile(archive_values, archive_high_quantile))
    trigger = (
        row["l2_improvement"] <= progress_cut
        and archive_low <= row["archive_gap"] <= archive_high
    )
    return trigger, {
        "progress_cut": progress_cut,
        "archive_low": archive_low,
        "archive_high": archive_high,
    }


def evaluate(rows, fixed_cuts=None):
    evaluated = []
    for index, row in enumerate(rows):
        if fixed_cuts is None:
            calibration = rows[:index] + rows[index + 1:]
            trigger, cuts = quantile_decision(row, calibration)
        else:
            cuts = dict(fixed_cuts)
            trigger = (
                row["l2_improvement"] <= cuts["progress_cut"]
                and cuts["archive_low"] <= row["archive_gap"]
                <= cuts["archive_high"]
            )
        adaptive_dv = (
            min(row["current_dv"], row["hill_valley_dv"])
            if trigger else row["current_dv"]
        )
        oracle_dv = min(row["current_dv"], row["hill_valley_dv"])
        adaptive_runtime = row["current_runtime"] + (
            row["hill_valley_downstream_runtime"] if trigger else 0.0
        )
        evaluated.append({
            **row,
            **cuts,
            "trigger_hill_valley": bool(trigger),
            "adaptive_dv": float(adaptive_dv),
            "oracle_dv": float(oracle_dv),
            "adaptive_runtime": float(adaptive_runtime),
            "hill_valley_relative_gain": (
                row["current_dv"] - row["hill_valley_dv"]
            ) / max(abs(row["current_dv"]), 1.0),
        })
    return evaluated


def summarize(rows):
    summary = {}
    for case in sorted({row["case"] for row in rows}):
        selected = [row for row in rows if row["case"] == case]
        current = [row["current_dv"] for row in selected]
        adaptive = [row["adaptive_dv"] for row in selected]
        oracle = [row["oracle_dv"] for row in selected]
        base_runtime = sum(row["current_runtime"] for row in selected)
        adaptive_runtime = sum(row["adaptive_runtime"] for row in selected)
        useful = [row for row in selected if row["hill_valley_relative_gain"] > 0.01]
        captured = [row for row in useful if row["trigger_hill_valley"]]
        summary[case] = {
            "runs": len(selected),
            "triggers": sum(row["trigger_hill_valley"] for row in selected),
            "useful_hill_valley_runs": len(useful),
            "useful_runs_captured": len(captured),
            "current_median_dv": float(np.median(current)),
            "adaptive_median_dv": float(np.median(adaptive)),
            "oracle_median_dv": float(np.median(oracle)),
            "adaptive_worst_dv": float(max(adaptive)),
            "oracle_regret_percent": float(
                100.0 * sum(a - o for a, o in zip(adaptive, oracle))
                / max(sum(oracle), 1.0)
            ),
            "projected_runtime_overhead_percent": float(
                100.0 * (adaptive_runtime - base_runtime)
                / max(base_runtime, 1.0)
            ),
        }
    all_base = sum(row["current_runtime"] for row in rows)
    all_adaptive = sum(row["adaptive_runtime"] for row in rows)
    summary["combined"] = {
        "runs": len(rows),
        "triggers": sum(row["trigger_hill_valley"] for row in rows),
        "useful_hill_valley_runs": sum(
            row["hill_valley_relative_gain"] > 0.01 for row in rows
        ),
        "useful_runs_captured": sum(
            row["hill_valley_relative_gain"] > 0.01
            and row["trigger_hill_valley"] for row in rows
        ),
        "projected_runtime_overhead_percent": float(
            100.0 * (all_adaptive - all_base) / max(all_base, 1.0)
        ),
    }
    return summary


def stress_equal_budget(rows, equal_budget_json):
    """Apply learned decisions to the independent identical-L0 replay."""
    payload = json.loads(Path(equal_budget_json).read_text(encoding="utf-8"))
    replay = {
        (row["case"], int(row["seed"]), row["branch"]): row
        for row in payload["rows"]
    }
    summary = {}
    for case in sorted({row["case"] for row in rows}):
        policy_rows = [row for row in rows if row["case"] == case]
        triples = [
            (
                replay[(case, row["seed"], "current")],
                replay[(case, row["seed"], "hill_valley")],
                row["trigger_hill_valley"],
            )
            for row in policy_rows
        ]
        current = [current["dv"] for current, _, _ in triples]
        adaptive = [
            min(current["dv"], hill_valley["dv"])
            if trigger else current["dv"]
            for current, hill_valley, trigger in triples
        ]
        oracle = [
            min(current["dv"], hill_valley["dv"])
            for current, hill_valley, _ in triples
        ]
        base_runtime = sum(current["total_runtime"] for current, _, _ in triples)
        extra_runtime = sum(
            hill_valley["total_runtime"]
            for _, hill_valley, trigger in triples if trigger
        )
        summary[case] = {
            "runs": len(triples),
            "triggers": sum(trigger for _, _, trigger in triples),
            "hill_valley_wins": sum(
                hill_valley["dv"] < current["dv"]
                for current, hill_valley, _ in triples
            ),
            "wins_captured": sum(
                trigger and hill_valley["dv"] < current["dv"]
                for current, hill_valley, trigger in triples
            ),
            "current_median_dv": float(np.median(current)),
            "adaptive_median_dv": float(np.median(adaptive)),
            "oracle_median_dv": float(np.median(oracle)),
            "projected_runtime_overhead_percent": float(
                100.0 * extra_runtime / max(base_runtime, 1.0)
            ),
        }
    return summary


def plot(rows, output_path, policy_label="leave-one-out"):
    fig, axes = plt.subplots(1, 3, figsize=(17, 5.2))
    ax = axes[0]
    full_progress_cut = float(np.median([row["progress_cut"] for row in rows]))
    archive_low = float(np.median([row["archive_low"] for row in rows]))
    archive_high = float(np.median([row["archive_high"] for row in rows]))
    for row in rows:
        useful = row["hill_valley_relative_gain"] > 0.01
        ax.scatter(
            row["l2_improvement"], row["archive_gap"],
            marker="*" if useful else "o",
            s=115 if useful else 55,
            color="#E45756" if row["trigger_hill_valley"] else "#4C78A8",
            edgecolor="black", linewidth=0.4,
        )
        ax.annotate(
            "{}{}".format(row["case"], row["seed"]),
            (row["l2_improvement"], row["archive_gap"]),
            xytext=(4, 3), textcoords="offset points", fontsize=7,
        )
    ax.axvline(full_progress_cut, color="black", linestyle="--", alpha=0.5)
    ax.axhspan(archive_low, archive_high, color="#F2CF5B", alpha=0.16)
    ax.set_xlabel("Relative MBH → L2 improvement")
    ax.set_ylabel("Exact archive best-to-median gap")
    ax.set_title("{} adaptive trigger\n(stars: useful HV >1%)".format(
        policy_label,
    ))
    ax.grid(True, alpha=0.2)

    ax = axes[1]
    cases = sorted({row["case"] for row in rows})
    x = np.arange(len(cases))
    width = 0.32
    values = {"current": [], "adaptive": []}
    for case in cases:
        selected = [row for row in rows if row["case"] == case]
        oracle_total = max(sum(row["oracle_dv"] for row in selected), 1.0)
        values["current"].append(
            100.0 * sum(
                row["current_dv"] - row["oracle_dv"] for row in selected
            ) / oracle_total
        )
        values["adaptive"].append(
            100.0 * sum(
                row["adaptive_dv"] - row["oracle_dv"] for row in selected
            ) / oracle_total
        )
    ax.bar(x - width / 2, values["current"], width, label="current", color="#4C78A8")
    ax.bar(x + width / 2, values["adaptive"], width, label="adaptive", color="#F58518")
    ax.set_xticks(x, cases)
    ax.set_ylabel("Aggregate regret versus oracle (%)")
    ax.set_title("Quality recovered from best-of-two oracle")
    ax.legend()
    ax.grid(True, axis="y", alpha=0.2)

    ax = axes[2]
    base = [sum(r["current_runtime"] for r in rows if r["case"] == c) for c in cases]
    adaptive = [sum(r["adaptive_runtime"] for r in rows if r["case"] == c) for c in cases]
    ax.bar(x - width / 2, base, width, label="current", color="#4C78A8")
    ax.bar(x + width / 2, adaptive, width, label="adaptive projected", color="#F58518")
    ax.set_xticks(x, cases)
    ax.set_ylabel("Aggregate runtime (s)")
    ax.set_title("Projected cost of conditional branches")
    ax.legend()
    ax.grid(True, axis="y", alpha=0.2)
    fig.tight_layout()
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=160)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--current-db",
        default=str(ROOT / "Tests" / "pressure_cascade_integrated_10seeds_funnel_l0_recall_64_mbh_between_pressure_cascade.sqlite"),
    )
    parser.add_argument(
        "--hill-valley-db",
        default=str(ROOT / "Tests" / "hill_valley_l0_integrated_10seeds_funnel_l0_recall_64_mbh_between_hill_valley_l0_pressure_cascade.sqlite"),
    )
    parser.add_argument(
        "--equal-budget-json",
        default=str(ROOT / "Tests" / "hill_valley_equal_budget_10seeds.json"),
    )
    parser.add_argument(
        "--output-json",
        default=str(ROOT / "Tests" / "adaptive_handoff_policy_loo_10seeds.json"),
    )
    parser.add_argument(
        "--output-plot",
        default=str(ROOT / "Tests" / "adaptive_handoff_policy_loo_10seeds.png"),
    )
    parser.add_argument(
        "--fixed-policy", action="store_true",
        help="Use the frozen seeds 0-9 calibration instead of leave-one-out.",
    )
    args = parser.parse_args()
    fixed_cuts = None
    if args.fixed_policy:
        fixed_cuts = {
            "progress_cut": 0.10628446377251824,
            "archive_low": 0.03394717357017691,
            "archive_high": 0.15573069674766785,
        }
    rows = evaluate(
        extract_rows(
            args.current_db, args.hill_valley_db, args.equal_budget_json,
        ),
        fixed_cuts=fixed_cuts,
    )
    payload = {
        "summary": summarize(rows),
        "policy": "fixed_seeds_0_9" if fixed_cuts else "leave_one_out",
        "fixed_cuts": fixed_cuts,
        "identical_l0_stress": (
            stress_equal_budget(rows, args.equal_budget_json)
            if args.equal_budget_json and Path(args.equal_budget_json).exists()
            and all(
                (row["case"], row["seed"], "current") in {
                    (item["case"], int(item["seed"]), item["branch"])
                    for item in json.loads(
                        Path(args.equal_budget_json).read_text(encoding="utf-8")
                    )["rows"]
                }
                for row in rows
            ) else None
        ),
        "rows": rows,
    }
    Path(args.output_json).write_text(
        json.dumps(payload, indent=2), encoding="utf-8",
    )
    plot(
        rows,
        args.output_plot,
        policy_label="Frozen seeds 0-9" if fixed_cuts else "Leave-one-out",
    )
    print(json.dumps(payload["summary"], indent=2))
    if payload["identical_l0_stress"] is not None:
        print("Identical-L0 stress:")
        print(json.dumps(payload["identical_l0_stress"], indent=2))
    print("JSON:", args.output_json)
    print("Plot:", args.output_plot)


if __name__ == "__main__":
    main()
