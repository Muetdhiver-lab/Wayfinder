# -*- coding: utf-8 -*-
"""Combine paired fixed-portfolio benchmark shards and plot the result."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
THRESHOLDS = {"KEEMo": 5500.0, "KEKKJ": 1500.0}


def load_rows(paths):
    rows = []
    for path in paths:
        rows.extend(json.loads(Path(path).read_text(encoding="utf-8"))["rows"])
    keys = [(row["case"], int(row["seed"]), row["branch"]) for row in rows]
    if len(keys) != len(set(keys)):
        raise ValueError("Duplicate case/seed/branch rows in portfolio inputs")
    return rows


def summarize(rows, portfolio_branch="portfolio_75_25"):
    branches = ("current", portfolio_branch)
    summary = {}
    for case in sorted({row["case"] for row in rows}):
        selected = [row for row in rows if row["case"] == case]
        by_key = {(row["seed"], row["branch"]): row for row in selected}
        seeds = sorted({row["seed"] for row in selected})
        if any((seed, branch) not in by_key for seed in seeds for branch in branches):
            raise ValueError("Unpaired portfolio rows for " + case)
        current = [by_key[(seed, "current")] for seed in seeds]
        portfolio = [by_key[(seed, portfolio_branch)] for seed in seeds]
        threshold = THRESHOLDS[case]
        summary[case] = {
            "runs": len(seeds),
            "current": {
                "median_dv": float(np.median([row["dv"] for row in current])),
                "best_dv": float(min(row["dv"] for row in current)),
                "worst_dv": float(max(row["dv"] for row in current)),
                "successes": sum(row["dv"] < threshold for row in current),
                "median_runtime": float(np.median([
                    row["total_runtime"] for row in current
                ])),
            },
            portfolio_branch: {
                "median_dv": float(np.median([row["dv"] for row in portfolio])),
                "best_dv": float(min(row["dv"] for row in portfolio)),
                "worst_dv": float(max(row["dv"] for row in portfolio)),
                "successes": sum(row["dv"] < threshold for row in portfolio),
                "median_runtime": float(np.median([
                    row["total_runtime"] for row in portfolio
                ])),
                "median_selection_runtime": float(np.median([
                    row["selection_runtime"] for row in portfolio
                ])),
            },
            "portfolio_wins": sum(
                portfolio_row["dv"] < current_row["dv"]
                for current_row, portfolio_row in zip(current, portfolio)
            ),
        }
    return summary


def plot(
    rows, output_path, portfolio_branch="portfolio_75_25",
    portfolio_label="portfolio 12+4",
):
    cases = ["KEEMo", "KEKKJ"]
    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    colors = {"current": "#4C78A8", portfolio_branch: "#54A24B"}
    for column, case in enumerate(cases):
        selected = [row for row in rows if row["case"] == case]
        lookup = {(row["seed"], row["branch"]): row for row in selected}
        seeds = sorted({row["seed"] for row in selected})
        ax = axes[0, column]
        for seed in seeds:
            ax.plot(
                [0, 1],
                [lookup[(seed, "current")]["dv"],
                 lookup[(seed, portfolio_branch)]["dv"]],
                color="#999999", alpha=0.55, marker="o", markersize=4,
            )
        ax.axhline(
            THRESHOLDS[case], color="#E45756", linestyle="--",
            alpha=0.75, label="success proxy",
        )
        ax.set_xticks([0, 1], ["current", portfolio_label])
        ax.set_ylabel("Objective DV (m/s)")
        ax.set_title("{} paired quality — 20 seeds".format(case))
        ax.grid(True, axis="y", alpha=0.22)
        ax.legend()

    summary = summarize(rows, portfolio_branch=portfolio_branch)
    ax = axes[1, 0]
    x = np.arange(len(cases))
    width = 0.34
    ax.bar(
        x - width / 2,
        [summary[case]["current"]["successes"] for case in cases],
        width, label="current", color=colors["current"],
    )
    ax.bar(
        x + width / 2,
        [summary[case][portfolio_branch]["successes"] for case in cases],
        width, label=portfolio_label, color=colors[portfolio_branch],
    )
    ax.set_xticks(x, cases)
    ax.set_ylim(0, 20)
    ax.set_ylabel("Runs below success proxy")
    ax.set_title("Consistency")
    ax.grid(True, axis="y", alpha=0.22)
    ax.legend()

    ax = axes[1, 1]
    ax.bar(
        x - width / 2,
        [summary[case]["current"]["median_runtime"] for case in cases],
        width, label="current", color=colors["current"],
    )
    ax.bar(
        x + width / 2,
        [summary[case][portfolio_branch]["median_runtime"] for case in cases],
        width, label=portfolio_label, color=colors[portfolio_branch],
    )
    ax.set_xticks(x, cases)
    ax.set_ylabel("Median selection + downstream runtime (s)")
    ax.set_title("Equal evolution budget runtime")
    ax.grid(True, axis="y", alpha=0.22)
    ax.legend()
    fig.tight_layout()
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=160)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "inputs", nargs="*",
        default=[
            str(ROOT / "Tests" / "portfolio_75_25_paired_seeds_0_9.json"),
            str(ROOT / "Tests" / "portfolio_75_25_paired_seeds_10_19.json"),
        ],
    )
    parser.add_argument(
        "--output-json",
        default=str(ROOT / "Tests" / "portfolio_75_25_paired_20seeds_summary.json"),
    )
    parser.add_argument(
        "--output-plot",
        default=str(ROOT / "Tests" / "portfolio_75_25_paired_20seeds.png"),
    )
    parser.add_argument("--portfolio-branch", default="portfolio_75_25")
    parser.add_argument("--portfolio-label", default="portfolio 12+4")
    args = parser.parse_args()
    rows = load_rows(args.inputs)
    payload = {
        "summary": summarize(rows, portfolio_branch=args.portfolio_branch),
        "portfolio_branch": args.portfolio_branch,
        "rows": rows,
    }
    Path(args.output_json).write_text(
        json.dumps(payload, indent=2), encoding="utf-8",
    )
    plot(
        rows, args.output_plot,
        portfolio_branch=args.portfolio_branch,
        portfolio_label=args.portfolio_label,
    )
    print(json.dumps(payload["summary"], indent=2))
    print("JSON:", args.output_json)
    print("Plot:", args.output_plot)


if __name__ == "__main__":
    main()
