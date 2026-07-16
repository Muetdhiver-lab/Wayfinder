# -*- coding: utf-8 -*-
"""Benchmark the explicit L0 recall funnel policies.

The benchmark intentionally keeps one SQLite database per optimizer strategy.
SQLite jobs are keyed by the search box and seed, while ``optimizer_strategy``
is a run-time option; sharing a DB across strategies would therefore de-dup the
work item and silently invalidate the comparison.
"""

import argparse
import sqlite3
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]
CORE_DIR = ROOT / "WayfinderCore"
sys.path.insert(0, str(CORE_DIR))

from _Wayfinder import Wayfinder  # noqa: E402


CASES = {
    "vanilla_kekkj": {
        "planet_pack": "Vanilla",
        "swing_by_bodies": [["Kerbin"], ["Eve"], ["Kerbin"], ["Kerbin"], ["Jool"]],
        "t0_min": 600,
        "t0_bin": 200,
        "n_t0_bins": 1,
        "tof_profile": "relaxed",
        "arrival_mode": "flyby",
    },
    "jnsq_keemo": {
        "planet_pack": "JNSQ",
        "swing_by_bodies": [["Kerbin"], ["Eve"], ["Eve"], ["Moho"]],
        "t0_min": 1460,
        "t0_bin": 365,
        "n_t0_bins": 1,
        "tof_profile": "relaxed",
        "arrival_mode": "vinf",
    },
    "jnsq_keemomo": {
        "planet_pack": "JNSQ",
        "swing_by_bodies": [["Kerbin"], ["Eve"], ["Eve"], ["Moho"], ["Moho"]],
        "t0_min": 1460,
        "t0_bin": 365,
        "n_t0_bins": 1,
        "tof_profile": "relaxed",
        "arrival_mode": "vinf",
    },
}


def parse_csv(value):
    return [item.strip() for item in str(value).split(",") if item.strip()]


def db_for_strategy(output_dir, label, strategy):
    db_strategy = "funnel_l0_recall" if strategy == "funnel_l0_recall_64" else strategy
    safe_strategy = db_strategy.replace("/", "_").replace("\\", "_")
    return output_dir / f"{label}_{safe_strategy}.sqlite"


def prepare_and_run(
    db_path, case_name, case, strategy, opt_level, seed, auto_workers,
    batch_name=None, leg_tof_bounds_override=None,
):
    plans = Wayfinder(planet_pack=case["planet_pack"])
    if batch_name is None:
        batch_name = f"l0bench_{case_name}_{strategy}_seed{seed}"
    plans.add_direct_t0_batch_sqlite(
        swing_by_bodies=case["swing_by_bodies"],
        db_path=db_path,
        batch_name=batch_name,
        t0_min=case["t0_min"],
        t0_bin=case["t0_bin"],
        n_t0_bins=case["n_t0_bins"],
        tof_profile=case["tof_profile"],
        opt_level=opt_level,
        arrival_mode=case["arrival_mode"],
        optimizer_topology="ring",
        optimizer_seed=seed,
        purpose="l0_recall_benchmark",
        leg_tof_bounds_override=leg_tof_bounds_override,
    )
    optimized = plans.optimize_sqlite(
        db_path,
        n=1,
        batch_name=batch_name,
        optimizer_strategy=strategy,
        auto_workers=auto_workers,
    )
    return optimized


def rows_for_db(db_path, strategy):
    if not db_path.exists():
        return []
    query = """
        select
            ? as strategy_label,
            s.short_name as sequence,
            j.planet_pack,
            j.t0_min,
            j.t0_max,
            j.tof_profile,
            j.arrival_mode,
            j.optimizer_seed as job_seed,
            r.optimizer_strategy as run_strategy,
            r.effective_optimizer_seed,
            r.runtime_seconds,
            r.actual_n_evo_steps,
            res.objective_dv,
            res.result_t0,
            res.result_tof,
            res.ejection_vinf,
            res.arrival_vinf
        from results res
        join runs r on r.id = res.run_id
        join jobs j on j.id = r.job_id
        join sequences s on s.id = j.sequence_id
        order by s.short_name, j.optimizer_seed, res.objective_dv
    """
    with sqlite3.connect(db_path) as con:
        con.row_factory = sqlite3.Row
        return [dict(row) for row in con.execute(query, (strategy,))]


def print_summary(rows):
    if not rows:
        print("No results found.")
        return
    headers = [
        "case", "strategy", "seed", "dv", "runtime", "steps",
        "t0", "tof", "ej_vinf", "arr_vinf",
    ]
    print(" ".join(f"{header:>14}" for header in headers))
    for row in sorted(rows, key=lambda r: (r["sequence"], r["job_seed"], r["objective_dv"])):
        print(
            f"{row['sequence']:>14} "
            f"{row['strategy_label']:>14} "
            f"{int(row['job_seed']):>14} "
            f"{row['objective_dv']:>14.3f} "
            f"{row['runtime_seconds']:>14.2f} "
            f"{int(row['actual_n_evo_steps']):>14} "
            f"{row['result_t0']:>14.1f} "
            f"{row['result_tof']:>14.1f} "
            f"{row['ejection_vinf']:>14.1f} "
            f"{row['arrival_vinf']:>14.1f}"
        )


def plot_summary(rows, output_path):
    if not rows:
        return None
    cases = sorted(set(row["sequence"] for row in rows))
    strategies = []
    for row in rows:
        strategy = row["strategy_label"]
        if strategy not in strategies:
            strategies.append(strategy)

    fig, axes = plt.subplots(
        len(cases), 2, figsize=(12, max(4, 3.2 * len(cases))), squeeze=False
    )
    colors = plt.get_cmap("tab10").colors
    color_by_strategy = {
        strategy: colors[index % len(colors)]
        for index, strategy in enumerate(strategies)
    }

    for row_index, case in enumerate(cases):
        case_rows = [row for row in rows if row["sequence"] == case]
        dv_data = []
        runtime_data = []
        labels = []
        for strategy in strategies:
            values = [row for row in case_rows if row["strategy_label"] == strategy]
            if not values:
                continue
            labels.append(strategy.replace("funnel_", ""))
            dv_data.append([row["objective_dv"] for row in values])
            runtime_data.append([row["runtime_seconds"] for row in values])

        for col_index, (data, ylabel, title_suffix) in enumerate((
            (dv_data, "Objective DV (m/s)", "quality"),
            (runtime_data, "Runtime (s)", "runtime"),
        )):
            ax = axes[row_index][col_index]
            box = ax.boxplot(data, tick_labels=labels, patch_artist=True)
            for patch, label in zip(box["boxes"], labels):
                strategy = "funnel_" + label if label != "funnel" else "funnel"
                patch.set_facecolor(color_by_strategy.get(strategy, colors[0]))
                patch.set_alpha(0.35)
            for x_index, values in enumerate(data, start=1):
                ax.scatter(
                    [x_index] * len(values),
                    values,
                    color="black",
                    s=18,
                    alpha=0.65,
                    zorder=3,
                )
            ax.set_title(f"{case} — {title_suffix}")
            ax.set_ylabel(ylabel)
            ax.grid(True, axis="y", alpha=0.25)
            ax.tick_params(axis="x", labelrotation=20)

    fig.tight_layout()
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=150)
    plt.close(fig)
    return output_path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--label", default="l0_recall_benchmark")
    parser.add_argument("--output-dir", default=str(ROOT / "Tests"))
    parser.add_argument(
        "--cases", default="vanilla_kekkj,jnsq_keemo",
        help="Comma-separated case names.",
    )
    parser.add_argument(
        "--strategies", default="funnel,funnel_l0_recall",
        help="Comma-separated optimizer_strategy values.",
    )
    parser.add_argument("--seeds", default="0")
    parser.add_argument("--opt-level", default="benchmark_funnel")
    parser.add_argument("--no-auto-workers", action="store_true")
    parser.add_argument(
        "--summary-only", action="store_true",
        help="Only read existing benchmark DBs.",
    )
    parser.add_argument("--plot", default=None, help="Optional output PNG path.")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    case_names = parse_csv(args.cases)
    strategies = parse_csv(args.strategies)
    seeds = [int(value) for value in parse_csv(args.seeds)]

    for case_name in case_names:
        if case_name not in CASES:
            raise SystemExit(f"Unknown case {case_name!r}. Valid: {sorted(CASES)}")

    if not args.summary_only:
        for strategy in strategies:
            db_path = db_for_strategy(output_dir, args.label, strategy)
            for case_name in case_names:
                for seed in seeds:
                    print(
                        f"Running {case_name} seed={seed} strategy={strategy} "
                        f"db={db_path.name}",
                        flush=True,
                    )
                    optimized = prepare_and_run(
                        db_path,
                        case_name,
                        CASES[case_name],
                        strategy,
                        args.opt_level,
                        seed,
                        auto_workers=not args.no_auto_workers,
                    )
                    print(f"Optimized {optimized} job(s)", flush=True)

    all_rows = []
    for strategy in strategies:
        db_path = db_for_strategy(output_dir, args.label, strategy)
        all_rows.extend(rows_for_db(db_path, strategy))
    print_summary(all_rows)
    if args.plot:
        output_path = plot_summary(all_rows, args.plot)
        if output_path is not None:
            print(f"Saved plot: {output_path}")


if __name__ == "__main__":
    main()
