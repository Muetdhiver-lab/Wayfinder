# -*- coding: utf-8 -*-
"""Small shared helpers for Wayfinder benchmark analysis scripts.

The benchmark scripts in ``Tests/`` are intentionally lightweight command-line
harnesses. This module centralizes the boring SQLite reads, case labels, summary
statistics, and box plots so new experiments do not keep reimplementing them.
"""

from __future__ import annotations

import sqlite3
import statistics
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


SEQUENCE_LABELS = {
    "Kerbin-Eve-Kerbin-Kerbin-Jool": "KEKKJ",
    "Kerbin-Eve-Eve-Moho": "KEEMo",
    "Kerbin-Eve-Eve-Moho-Moho": "KEEMoMo",
}

CASE_BATCH_NAMES = {
    "KEKKJ": "vanilla_kekkj",
    "KEEMo": "jnsq_keemo",
    "KEEMoMo": "jnsq_keemomo",
}


def case_from_sequence(sequence):
    return SEQUENCE_LABELS.get(sequence, sequence)


def connect_rows(db_path, query, params=()):
    path = Path(db_path)
    if not path.exists():
        return []
    with sqlite3.connect(path) as con:
        con.row_factory = sqlite3.Row
        return [dict(row) for row in con.execute(query, tuple(params))]


def optimizer_result_rows(
    db_path,
    strategy_label,
    sequence_short_name=None,
    include_batch=True,
):
    batch_column = "b.name as batch," if include_batch else "NULL as batch,"
    where = []
    params = [strategy_label]
    if sequence_short_name is not None:
        where.append("s.short_name = ?")
        params.append(sequence_short_name)
    where_sql = "where " + " and ".join(where) if where else ""
    query = f"""
        select
            ? as strategy_label,
            {batch_column}
            s.short_name as sequence,
            j.planet_pack,
            j.optimizer_seed as seed,
            r.id as run_id,
            r.runtime_seconds as runtime_seconds,
            r.actual_n_evo_steps as actual_n_evo_steps,
            res.objective_dv as objective_dv,
            res.result_t0 as result_t0,
            res.result_tof as result_tof,
            res.ejection_vinf as ejection_vinf,
            res.arrival_vinf as arrival_vinf
        from results res
        join runs r on r.id = res.run_id
        join jobs j on j.id = r.job_id
        join sequences s on s.id = j.sequence_id
        join batch_jobs bj on bj.job_id = j.id
        join batches b on b.id = bj.batch_id
        {where_sql}
        order by s.short_name, j.optimizer_seed, res.objective_dv
    """
    return connect_rows(db_path, query, params)


def summarize_values(values, threshold=None):
    values = [float(value) for value in values]
    if not values:
        return {
            "n": 0,
            "median": None,
            "mean": None,
            "min": None,
            "max": None,
            "below_threshold": None,
        }
    return {
        "n": len(values),
        "median": statistics.median(values),
        "mean": statistics.mean(values),
        "min": min(values),
        "max": max(values),
        "below_threshold": (
            sum(value < float(threshold) for value in values)
            if threshold is not None else None
        ),
    }


def grouped_summary(rows, group_key="strategy_label", threshold_by_group=None):
    threshold_by_group = threshold_by_group or {}
    groups = {}
    for row in rows:
        groups.setdefault(row[group_key], []).append(row)
    summaries = {}
    for group, group_rows in sorted(groups.items()):
        dv_values = [float(row["objective_dv"]) for row in group_rows]
        runtime_values = [float(row["runtime_seconds"]) for row in group_rows]
        summaries[group] = {
            "objective_dv": summarize_values(
                dv_values, threshold=threshold_by_group.get(group)
            ),
            "runtime_seconds": summarize_values(runtime_values),
        }
    return summaries


def print_result_table(rows):
    if not rows:
        print("No results found.")
        return
    print(
        f"{'strategy':>16} {'case':>10} {'seed':>6} {'dv':>10} "
        f"{'runtime':>9} {'steps':>6} {'t0':>8} {'tof':>8} {'arr_vinf':>10}"
    )
    for row in rows:
        print(
            f"{row['strategy_label']:>16} "
            f"{case_from_sequence(row['sequence']):>10} "
            f"{int(row['seed']):>6} "
            f"{float(row['objective_dv']):>10.3f} "
            f"{float(row['runtime_seconds']):>9.2f} "
            f"{int(row['actual_n_evo_steps'] or 0):>6} "
            f"{float(row['result_t0']):>8.1f} "
            f"{float(row['result_tof']):>8.1f} "
            f"{float(row['arrival_vinf'] or 0.0):>10.1f}"
        )


def print_grouped_summary(rows, threshold_by_group=None):
    for group, summary in grouped_summary(
        rows, threshold_by_group=threshold_by_group,
    ).items():
        dv = summary["objective_dv"]
        runtime = summary["runtime_seconds"]
        success = ""
        if dv["below_threshold"] is not None:
            success = f" below_threshold={dv['below_threshold']}"
        print(
            f"{group}: n={dv['n']} median={dv['median']:.1f} "
            f"mean={dv['mean']:.1f} min={dv['min']:.1f} max={dv['max']:.1f}"
            f"{success} median_runtime={runtime['median']:.2f}s"
        )


def plot_quality_runtime_boxplots(
    rows,
    output_path,
    group_order=None,
    group_key="strategy_label",
    title_prefix="",
):
    if not rows:
        return None
    if group_order is None:
        group_order = []
        for row in rows:
            group = row[group_key]
            if group not in group_order:
                group_order.append(group)
    groups = [
        group for group in group_order
        if any(row[group_key] == group for row in rows)
    ]
    dv_data = [
        [float(row["objective_dv"]) for row in rows if row[group_key] == group]
        for group in groups
    ]
    runtime_data = [
        [float(row["runtime_seconds"]) for row in rows if row[group_key] == group]
        for group in groups
    ]
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.2))
    colors = ["#4C78A8", "#F58518", "#54A24B", "#B279A2", "#E45756"]
    for ax, data, ylabel, title in (
        (axes[0], dv_data, "Objective DV (m/s)", "quality"),
        (axes[1], runtime_data, "Runtime (s)", "runtime"),
    ):
        box = ax.boxplot(data, tick_labels=groups, patch_artist=True)
        for patch, color in zip(box["boxes"], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.35)
        for x_index, values in enumerate(data, start=1):
            ax.scatter(
                [x_index] * len(values),
                values,
                c="black",
                s=20,
                alpha=0.65,
                zorder=3,
            )
        prefix = (str(title_prefix).strip() + " ") if title_prefix else ""
        ax.set_title(prefix + title)
        ax.set_ylabel(ylabel)
        ax.grid(True, axis="y", alpha=0.25)
        ax.tick_params(axis="x", labelrotation=20)
    fig.tight_layout()
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=150)
    plt.close(fig)
    return output_path


def plot_quality_runtime_by_case(
    rows,
    output_path,
    strategy_order=None,
    case_order=None,
    threshold_by_case=None,
    title="",
):
    """Plot comparable DV/runtime distributions for each benchmark case."""
    if not rows:
        return None
    threshold_by_case = threshold_by_case or {}
    if strategy_order is None:
        strategy_order = list(dict.fromkeys(
            row["strategy_label"] for row in rows
        ))
    if case_order is None:
        case_order = list(dict.fromkeys(
            case_from_sequence(row["sequence"]) for row in rows
        ))
    fig, axes = plt.subplots(
        len(case_order), 2,
        figsize=(12.5, max(4.2, 3.8 * len(case_order))),
        squeeze=False,
    )
    colors = ["#4C78A8", "#F58518", "#54A24B", "#B279A2"]
    for case_index, case in enumerate(case_order):
        case_rows = [
            row for row in rows
            if case_from_sequence(row["sequence"]) == case
        ]
        for column, field, ylabel in (
            (0, "objective_dv", "Objective DV (m/s)"),
            (1, "runtime_seconds", "Runtime (s)"),
        ):
            ax = axes[case_index][column]
            data = [
                [
                    float(row[field]) for row in case_rows
                    if row["strategy_label"] == strategy
                ]
                for strategy in strategy_order
            ]
            box = ax.boxplot(
                data, tick_labels=strategy_order, patch_artist=True,
            )
            for patch, color in zip(box["boxes"], colors):
                patch.set_facecolor(color)
                patch.set_alpha(0.35)
            for x_index, values in enumerate(data, start=1):
                ax.scatter(
                    [x_index] * len(values), values,
                    c="black", s=18, alpha=0.65, zorder=3,
                )
            if column == 0 and case in threshold_by_case:
                ax.axhline(
                    float(threshold_by_case[case]), color="#E45756",
                    linestyle="--", linewidth=1.2, alpha=0.8,
                )
            ax.set_title("{} {}".format(case, "quality" if column == 0 else "runtime"))
            ax.set_ylabel(ylabel)
            ax.grid(True, axis="y", alpha=0.25)
            ax.tick_params(axis="x", labelrotation=18)
    if title:
        fig.suptitle(str(title))
    fig.tight_layout()
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=150)
    plt.close(fig)
    return output_path


def plot_sequence_scout_funnel(
    rows,
    output_path,
    direct_reference_dv=None,
    title="Tisserand/Lambert scout -> L0 -> full funnel",
):
    """Plot bin quality, L0-to-final gain, runtime, and best convergence.

    The helper consumes the JSON row format emitted by
    ``run_sequence_scout_full_funnel_benchmark.py``. Keeping this here makes
    future sequence-scout experiments comparable without duplicating plotting
    code in each harness.
    """
    if not rows:
        return None
    rows = sorted(
        rows,
        key=lambda row: (
            float(row["bin_start_days"]), float(row["objective_dv"]),
        ),
    )
    fig, axes = plt.subplots(2, 2, figsize=(13.5, 8.2))
    quality_ax, gain_ax, runtime_ax, convergence_ax = axes.ravel()
    colors = plt.get_cmap("tab20")

    for index, row in enumerate(rows):
        x = 0.5 * (
            float(row["bin_start_days"]) + float(row["bin_end_days"])
        )
        final = float(row["objective_dv"])
        label = row.get("sequence_short_name", row["sequence"])
        color = colors(index % 20)
        quality_ax.scatter(x, final, color=color, s=34, zorder=3)
        quality_ax.annotate(
            label, (x, final), xytext=(3, 4), textcoords="offset points",
            fontsize=7, rotation=18,
        )
        l0 = float(row["screening_l0_dv"])
        gain_ax.plot([l0, final], [index, index], color=color, linewidth=1.3)
        gain_ax.scatter(l0, index, color="#9ECAE1", s=22, zorder=3)
        gain_ax.scatter(final, index, color="#08519C", s=26, zorder=3)
        runtime_ax.scatter(
            x, float(row["runtime_seconds"]), color=color, s=34, zorder=3,
        )

    if direct_reference_dv is not None:
        for ax in (quality_ax, gain_ax):
            if ax is quality_ax:
                ax.axhline(
                    float(direct_reference_dv), color="#E45756",
                    linestyle="--", linewidth=1.2,
                    label="best direct ejection",
                )
            else:
                ax.axvline(
                    float(direct_reference_dv), color="#E45756",
                    linestyle="--", linewidth=1.2,
                )
        quality_ax.legend(loc="best", fontsize=8)

    quality_ax.set_title("Final quality by departure bin")
    quality_ax.set_xlabel("T0 bin midpoint (Kerbal days)")
    quality_ax.set_ylabel("Final objective DV (m/s)")
    quality_ax.grid(True, alpha=0.25)

    gain_ax.set_title("Promotion gain: screening L0 -> final")
    gain_ax.set_xlabel("Objective DV (m/s)")
    gain_ax.set_yticks(range(len(rows)))
    gain_ax.set_yticklabels([
        "{} @ {:.0f}".format(
            row.get("sequence_short_name", row["sequence"]),
            float(row["bin_start_days"]),
        )
        for row in rows
    ], fontsize=7)
    gain_ax.invert_yaxis()
    gain_ax.grid(True, axis="x", alpha=0.25)

    runtime_ax.set_title("Full-funnel runtime")
    runtime_ax.set_xlabel("T0 bin midpoint (Kerbal days)")
    runtime_ax.set_ylabel("Runtime (s)")
    runtime_ax.grid(True, alpha=0.25)

    best_row = min(rows, key=lambda row: float(row["objective_dv"]))
    snapshots = best_row.get("convergence", [])
    if snapshots:
        steps = [int(point["step"]) for point in snapshots]
        convergence_ax.plot(
            steps, [float(point["best"]) for point in snapshots],
            color="#1B9E77", linewidth=2.0, label="best island",
        )
        convergence_ax.plot(
            steps,
            [float(point["average_champion"]) for point in snapshots],
            color="#7570B3", linewidth=1.3, linestyle="--",
            label="average champion",
        )
        boundary = 0
        for summary in best_row.get("stage_summaries", [])[:-1]:
            boundary += int(summary.get("actual_evo_steps", 0))
            convergence_ax.axvline(
                boundary, color="#777777", alpha=0.35, linewidth=0.8,
            )
        convergence_ax.set_yscale("log")
        convergence_ax.legend(loc="best", fontsize=8)
    convergence_ax.set_title(
        "Best run convergence: {} ({:.1f} m/s)".format(
            best_row.get("sequence_short_name", best_row["sequence"]),
            float(best_row["objective_dv"]),
        )
    )
    convergence_ax.set_xlabel("Cumulative evolution step")
    convergence_ax.set_ylabel("Objective DV (m/s, log scale)")
    convergence_ax.grid(True, alpha=0.25)

    fig.suptitle(title)
    fig.tight_layout()
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=160)
    plt.close(fig)
    return output_path
