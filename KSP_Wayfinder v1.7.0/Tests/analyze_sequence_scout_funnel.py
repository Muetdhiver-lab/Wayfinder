# -*- coding: utf-8 -*-
"""Summarize and plot a sequence-scout L0/full-funnel benchmark."""

import argparse
import json
import statistics
from pathlib import Path

from benchmark_analysis import plot_sequence_scout_funnel


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input", default="Tests/sequence_scout_full_kj_0_1000.json",
    )
    parser.add_argument(
        "--l0-input", default="Tests/sequence_scout_l0_kj_0_1000.json",
    )
    parser.add_argument(
        "--plot", default="Tests/sequence_scout_full_kj_0_1000.png",
    )
    args = parser.parse_args()

    payload = json.loads(Path(args.input).read_text(encoding="utf-8"))
    l0_payload = json.loads(Path(args.l0_input).read_text(encoding="utf-8"))
    rows = payload["rows"]
    direct = float(l0_payload["direct_reference"]["ejection_dv_mps"])
    plot_sequence_scout_funnel(rows, args.plot, direct_reference_dv=direct)

    best = min(rows, key=lambda row: float(row["objective_dv"]))
    runtimes = [float(row["runtime_seconds"]) for row in rows]
    l0_runtime = sum(
        float(row["l0_runtime_seconds"]) for row in l0_payload["rows"]
    )
    full_runtime = sum(runtimes)
    below_direct = [row for row in rows if float(row["objective_dv"]) < direct]
    print("funnels={}".format(len(rows)))
    print("best={} {:.1f} m/s T0 {:.0f}-{:.0f}".format(
        best["sequence_short_name"], float(best["objective_dv"]),
        float(best["bin_start_days"]), float(best["bin_end_days"]),
    ))
    print("median_dv={:.1f} m/s".format(statistics.median(
        float(row["objective_dv"]) for row in rows
    )))
    print("below_direct={}/{}".format(len(below_direct), len(rows)))
    print("median_funnel_runtime={:.2f} s".format(statistics.median(runtimes)))
    print("l0_total_runtime={:.2f} s".format(l0_runtime))
    print("full_total_runtime={:.2f} s".format(full_runtime))
    print("total_runtime={:.2f} s".format(l0_runtime + full_runtime))
    print("plot={}".format(args.plot))


if __name__ == "__main__":
    main()
