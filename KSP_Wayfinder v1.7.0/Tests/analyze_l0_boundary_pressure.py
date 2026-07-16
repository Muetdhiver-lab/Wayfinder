# -*- coding: utf-8 -*-
"""Report leg-TOF boundary pressure in stored optimizer populations.

This diagnostic is intentionally read-only: it inspects SQLite optimizer
population archives and reports whether the best candidates are crowded near
the configured per-leg TOF bounds. It is useful before comparing optimizers:
if L0 already presses against a bound, the search box is likely suspect.
"""

from __future__ import annotations

import argparse
import json
import sqlite3
import statistics
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "WayfinderCore"))

from _Optimization import leg_tofs  # noqa: E402
from _Wayfinder import Wayfinder  # noqa: E402


def parse_csv_ints(value):
    return [int(part.strip()) for part in str(value).split(",") if part.strip()]


def percentile(sorted_values, fraction):
    if not sorted_values:
        return None
    index = int(round(float(fraction) * (len(sorted_values) - 1)))
    index = max(0, min(len(sorted_values) - 1, index))
    return sorted_values[index]


def run_rows(con, sequence):
    return list(con.execute(
        """
        select
            r.id as run_id,
            j.optimizer_seed as seed,
            j.planet_pack,
            j.tof_encoding,
            j.leg_tof_bounds_json
        from runs r
        join jobs j on j.id = r.job_id
        join sequences s on s.id = j.sequence_id
        where s.short_name = ?
        order by j.optimizer_seed, r.id
        """,
        (sequence,),
    ))


def population_rows(con, run_id, source):
    return list(con.execute(
        """
        select fitness, gene_json
        from optimizer_population_points
        where run_id = ? and source = ?
        order by fitness asc
        """,
        (run_id, source),
    ))


def pressure_summary(tof_vectors, bounds, near_fraction):
    summaries = []
    for leg_index, (low, high) in enumerate(bounds):
        span = float(high) - float(low)
        values = sorted(float(vector[leg_index]) for vector in tof_vectors)
        low_gaps = [(value - low) / span for value in values]
        high_gaps = [(high - value) / span for value in values]
        summaries.append({
            "leg": leg_index + 1,
            "low": float(low),
            "high": float(high),
            "min": min(values),
            "p05": percentile(values, 0.05),
            "median": statistics.median(values),
            "p95": percentile(values, 0.95),
            "max": max(values),
            "near_low": sum(gap <= near_fraction for gap in low_gaps),
            "near_high": sum(gap <= near_fraction for gap in high_gaps),
            "min_low_gap_days": min(value - low for value in values),
            "min_high_gap_days": min(high - value for value in values),
        })
    return summaries


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("db_path")
    parser.add_argument("--sequence", default="KEEMo")
    parser.add_argument("--source", default="stage_1_final")
    parser.add_argument("--top", default="16,32")
    parser.add_argument(
        "--near-fraction", type=float, default=0.03,
        help="Normalized distance-to-bound threshold counted as pressure.",
    )
    args = parser.parse_args()

    db_path = Path(args.db_path)
    top_counts = parse_csv_ints(args.top)

    con = sqlite3.connect(db_path)
    con.row_factory = sqlite3.Row
    try:
        runs = run_rows(con, args.sequence)
        if not runs:
            raise SystemExit(f"No runs found for sequence {args.sequence!r}")

        factor_by_pack = {}
        for run in runs:
            factor_by_pack.setdefault(
                run["planet_pack"],
                Wayfinder(planet_pack=run["planet_pack"])._Edy2Kdy,
            )

        print(f"DB: {db_path}")
        print(f"Sequence: {args.sequence}")
        print(f"Source: {args.source}")
        print(f"Near-bound threshold: {args.near_fraction:.3f} of leg span")

        for run in runs:
            bounds = json.loads(run["leg_tof_bounds_json"])
            factor = factor_by_pack[run["planet_pack"]]
            rows = population_rows(con, run["run_id"], args.source)
            if not rows:
                print(f"\nseed={run['seed']} no population rows")
                continue

            for top_count in top_counts:
                top_rows = rows[:top_count]
                tof_vectors = []
                for row in top_rows:
                    gene = json.loads(row["gene_json"])
                    tof_vectors.append([
                        value * factor
                        for value in leg_tofs(
                            gene, len(bounds), run["tof_encoding"],
                        )
                    ])
                print(
                    f"\nseed={run['seed']} top={top_count} "
                    f"bestfit={top_rows[0]['fitness']:.3f} "
                    f"top1={[round(value, 1) for value in tof_vectors[0]]}"
                )
                for summary in pressure_summary(
                    tof_vectors, bounds, args.near_fraction,
                ):
                    print(
                        "  leg{leg} [{low:.1f},{high:.1f}] "
                        "min/p05/med/p95/max="
                        "{min:.1f}/{p05:.1f}/{median:.1f}/{p95:.1f}/{max:.1f} "
                        "near_low={near_low}/{top} near_high={near_high}/{top} "
                        "gap_low={min_low_gap_days:.2f}d "
                        "gap_high={min_high_gap_days:.2f}d".format(
                            top=top_count, **summary,
                        )
                    )
    finally:
        con.close()


if __name__ == "__main__":
    main()
