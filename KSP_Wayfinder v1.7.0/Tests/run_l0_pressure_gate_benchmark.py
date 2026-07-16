# -*- coding: utf-8 -*-
"""Benchmark a simple L0 boundary-pressure gate.

The gate is deliberately implemented as a benchmark harness rather than as a
core optimizer mutation. It reuses a normal baseline run, inspects the stored
L0 population, and reruns the same optimizer with one-sided relaxed leg-TOF
bounds only when the L0 elites visibly press against a boundary.

Runtime reporting distinguishes:

* baseline_runtime: the normal completed funnel runtime;
* l0_plus_rescue_runtime_projected: what a replace-by-rescue implementation
  would roughly pay (baseline L0 stage + rescue rerun when pressure triggers);
* gate_runtime_projected: the safer branch-and-keep-best runtime
  (baseline full funnel + relaxed rerun when pressure triggers);
* gate_runtime_harness: what this harness actually paid when it reused an
  already completed baseline funnel.
"""

from __future__ import annotations

import argparse
import json
import sqlite3
import statistics
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "WayfinderCore"))
sys.path.insert(0, str(ROOT / "Tests"))

from _Optimization import leg_tofs  # noqa: E402
from _Wayfinder import Wayfinder  # noqa: E402
from run_l0_recall_benchmark import CASES  # noqa: E402
from run_l0_recall_benchmark import db_for_strategy  # noqa: E402
from run_l0_recall_benchmark import parse_csv  # noqa: E402
from run_l0_recall_benchmark import prepare_and_run  # noqa: E402


def _batch_name(case_name, strategy, seed):
    return f"l0bench_{case_name}_{strategy}_seed{seed}"


def _pressure_batch_name(
    case_name, strategy, seed, top, near_fraction, relax_fraction,
    adjustment_mode,
):
    return (
        f"l0pressure_{case_name}_{strategy}_seed{seed}"
        f"_top{top}_nf{int(round(near_fraction * 1000))}"
        f"_rf{int(round(relax_fraction * 100))}"
        f"_{adjustment_mode}"
    )


def _legacy_pressure_batch_name(
    case_name, strategy, seed, top, near_fraction, relax_fraction,
):
    return (
        f"l0pressure_{case_name}_{strategy}_seed{seed}"
        f"_top{top}_nf{int(round(near_fraction * 1000))}"
        f"_rf{int(round(relax_fraction * 100))}"
    )


def _retry_pressure_batch_name(
    case_name, strategy, seed, top, near_fraction, relax_fraction,
    adjustment_mode, seed_offset,
):
    return (
        f"l0retry_{case_name}_{strategy}_seed{seed}"
        f"_plus{int(seed_offset)}"
        f"_top{top}_nf{int(round(near_fraction * 1000))}"
        f"_rf{int(round(relax_fraction * 100))}"
        f"_{adjustment_mode}"
    )


def _sequence_short_name(case):
    plans = Wayfinder(planet_pack=case["planet_pack"])
    return plans.generateShortSequences(case["swing_by_bodies"])[0]


def _result_for_batch(db_path, batch_name, sequence_short_name, seed):
    if not Path(db_path).exists():
        return None
    query = """
        select
            r.id as run_id,
            j.id as job_id,
            j.optimizer_seed as seed,
            s.short_name as sequence,
            j.planet_pack,
            j.tof_encoding,
            j.leg_tof_bounds_json,
            r.runtime_seconds,
            res.objective_dv,
            res.result_t0,
            res.result_tof,
            res.ejection_vinf,
            res.arrival_vinf
        from results res
        join runs r on r.id = res.run_id
        join jobs j on j.id = r.job_id
        join sequences s on s.id = j.sequence_id
        join batch_jobs bj on bj.job_id = j.id
        join batches b on b.id = bj.batch_id
        where b.name = ?
          and s.short_name = ?
          and j.optimizer_seed = ?
        order by r.id desc
        limit 1
    """
    with sqlite3.connect(db_path) as con:
        con.row_factory = sqlite3.Row
        row = con.execute(
            query, (batch_name, sequence_short_name, int(seed)),
        ).fetchone()
        return dict(row) if row else None


def _stage_runtime(db_path, run_id, stage_index=1):
    query = """
        select runtime_seconds
        from optimizer_stages
        where run_id = ? and stage_index = ?
    """
    with sqlite3.connect(db_path) as con:
        row = con.execute(query, (int(run_id), int(stage_index))).fetchone()
        return float(row[0]) if row else 0.0


def _stage_best(db_path, run_id, stage_index):
    query = """
        select best_fitness
        from optimizer_stages
        where run_id = ? and stage_index = ?
    """
    with sqlite3.connect(db_path) as con:
        row = con.execute(query, (int(run_id), int(stage_index))).fetchone()
        return float(row[0]) if row else None


def _population_tofs(db_path, run, source, top):
    bounds = json.loads(run["leg_tof_bounds_json"])
    factor = Wayfinder(planet_pack=run["planet_pack"])._Edy2Kdy
    query = """
        select fitness, gene_json
        from optimizer_population_points
        where run_id = ? and source = ?
        order by fitness asc
        limit ?
    """
    with sqlite3.connect(db_path) as con:
        con.row_factory = sqlite3.Row
        rows = list(con.execute(query, (int(run["run_id"]), source, int(top))))
    vectors = []
    for row in rows:
        gene = json.loads(row["gene_json"])
        vectors.append([
            float(value) * factor
            for value in leg_tofs(gene, len(bounds), run["tof_encoding"])
        ])
    return bounds, vectors


def _population_fitnesses(db_path, run_id, source, top):
    query = """
        select fitness
        from optimizer_population_points
        where run_id = ? and source = ?
        order by fitness asc
        limit ?
    """
    with sqlite3.connect(db_path) as con:
        return [
            float(row[0])
            for row in con.execute(query, (int(run_id), source, int(top)))
        ]


def _percentile(sorted_values, fraction):
    if not sorted_values:
        return None
    index = int(round(float(fraction) * (len(sorted_values) - 1)))
    return sorted_values[max(0, min(len(sorted_values) - 1, index))]


def pressure_actions(
    bounds,
    tof_vectors,
    near_fraction=0.03,
    min_pressure_count=2,
):
    """Return one-sided pressure actions for each leg.

    The p05/p95 rule catches a real flank even if only a few points land inside
    the strict near-bound bucket. This intentionally biases toward recall: a
    false positive costs one relaxed rerun, while a false negative can lose the
    useful basin.
    """
    actions = []
    if not tof_vectors:
        return actions
    for leg_index, (low, high) in enumerate(bounds):
        span = max(float(high) - float(low), 1e-12)
        values = sorted(float(vector[leg_index]) for vector in tof_vectors)
        near_low = sum((value - low) / span <= near_fraction for value in values)
        near_high = sum((high - value) / span <= near_fraction for value in values)
        p05 = _percentile(values, 0.05)
        p95 = _percentile(values, 0.95)
        low_pressed = (
            near_low >= min_pressure_count
            or (p05 is not None and (p05 - low) / span <= near_fraction)
        )
        high_pressed = (
            near_high >= min_pressure_count
            or (p95 is not None and (high - p95) / span <= near_fraction)
        )
        if low_pressed:
            actions.append({
                "leg": leg_index + 1,
                "side": "low",
                "near_count": near_low,
                "p05": p05,
                "p95": p95,
            })
        if high_pressed:
            actions.append({
                "leg": leg_index + 1,
                "side": "high",
                "near_count": near_high,
                "p05": p05,
                "p95": p95,
            })
    return actions


def adjust_bounds(
    bounds,
    actions,
    relax_fraction=0.20,
    min_relax_days=5.0,
    max_relax_days=60.0,
    min_leg_lower_days=1.0,
    adjustment_mode="widen",
):
    adjusted = [[float(low), float(high)] for low, high in bounds]
    deltas = {}
    modes = {}
    actions_by_leg = {}
    for action in actions:
        actions_by_leg.setdefault(int(action["leg"]), set()).add(action["side"])

    for leg, sides in sorted(actions_by_leg.items()):
        index = leg - 1
        low, high = adjusted[index]
        span = max(high - low, 1e-12)
        delta = max(
            float(min_relax_days),
            min(float(max_relax_days), span * float(relax_fraction)),
        )

        if adjustment_mode == "shift" and sides == {"low"}:
            actual_delta = min(delta, low - float(min_leg_lower_days))
            if actual_delta > 0:
                adjusted[index][0] = low - actual_delta
                adjusted[index][1] = high - actual_delta
            else:
                adjusted[index][0] = max(float(min_leg_lower_days), low - delta)
            deltas[(leg, "low")] = actual_delta if actual_delta > 0 else delta
            modes[(leg, "low")] = "shift"
        elif adjustment_mode == "shift" and sides == {"high"}:
            adjusted[index][0] = low + delta
            adjusted[index][1] = high + delta
            deltas[(leg, "high")] = delta
            modes[(leg, "high")] = "shift"
        else:
            # Either explicit widen mode, or ambiguous pressure on both flanks.
            if "low" in sides:
                adjusted[index][0] = max(float(min_leg_lower_days), low - delta)
                deltas[(leg, "low")] = delta
                modes[(leg, "low")] = "widen"
            if "high" in sides:
                adjusted[index][1] = high + delta
                deltas[(leg, "high")] = delta
                modes[(leg, "high")] = "widen"
    return adjusted, deltas, modes


def _actions_label(actions, deltas, modes):
    if not actions:
        return "-"
    chunks = []
    for action in actions:
        leg = int(action["leg"])
        side = "-" if action["side"] == "low" else "+"
        delta = deltas.get((leg, action["side"]), 0.0)
        mode = modes.get((leg, action["side"]), "widen")
        marker = "S" if mode == "shift" else "W"
        chunks.append(f"L{leg}{side}{delta:.1f}d{marker}")
    return ",".join(chunks)


def _parse_thresholds(value):
    thresholds = {}
    for item in parse_csv(value):
        if ":" not in item:
            raise ValueError(
                "Quality thresholds must use SEQUENCE:DV entries, got: "
                + str(item)
            )
        key, raw_threshold = item.split(":", 1)
        thresholds[key.strip()] = float(raw_threshold)
    return thresholds


def _should_branch(
    policy, sequence, baseline_dv, actions, thresholds,
    l1_improvement=None, max_l1_improvement=0.15,
    l1_best_to_l0_median=None, min_l1_best_to_l0_median=0.725,
):
    if not actions:
        return False
    if policy == "pressure":
        return True
    if policy == "quality_guard":
        threshold = thresholds.get(sequence)
        if threshold is None:
            return False
        return float(baseline_dv) > float(threshold)
    if policy == "l1_improvement":
        if l1_improvement is None:
            return False
        return float(l1_improvement) <= float(max_l1_improvement)
    if policy == "l1_combined":
        poor_improvement = (
            l1_improvement is not None
            and float(l1_improvement) <= float(max_l1_improvement)
        )
        high_relative_l1 = (
            l1_best_to_l0_median is not None
            and float(l1_best_to_l0_median)
            >= float(min_l1_best_to_l0_median)
        )
        return bool(poor_improvement or high_relative_l1)
    raise ValueError("Unknown branch policy: " + str(policy))


def _should_retry_after_same_seed(
    sequence,
    baseline_dv,
    current_best_dv,
    relaxed_dv,
    thresholds,
    min_improvement,
):
    """Return True when the first rescue still looks suspicious.

    If an explicit per-sequence quality threshold is available, use it. This is
    the cleanest benchmark interpretation: pay the second rescue only if the
    same-seed branch did not get the case back into the acceptable basin.

    Without a threshold, fall back to a relative-improvement guard so the mode
    still has deterministic behaviour in exploratory runs.
    """
    threshold = thresholds.get(sequence)
    if threshold is not None:
        return float(current_best_dv) > float(threshold)
    if relaxed_dv is None:
        return True
    if float(baseline_dv) <= 0.0:
        return False
    improvement = (float(baseline_dv) - float(current_best_dv)) / float(baseline_dv)
    return improvement <= float(min_improvement)


def _print_summary(rows):
    headers = [
        "case", "seed", "pressure", "branch", "l1_imp", "l1/l0med",
        "baseline_dv", "relaxed_dv", "retry_dv", "best_dv", "base_rt",
        "branch_rt", "actions",
    ]
    print(" ".join(f"{header:>14}" for header in headers))
    for row in rows:
        relaxed_dv = (
            f"{row['relaxed_dv']:14.3f}"
            if row["relaxed_dv"] is not None
            else f"{'-':>14}"
        )
        retry_dv = (
            f"{row['retry_dv']:14.3f}"
            if row["retry_dv"] is not None
            else f"{'-':>14}"
        )
        print(
            f"{row['case']:>14} "
            f"{int(row['seed']):>14} "
            f"{str(row['pressure_detected']):>14} "
            f"{str(row['triggered']):>14} "
            f"{row['l1_improvement']:>14.3f} "
            f"{row['l1_best_to_l0_median']:>14.3f} "
            f"{row['baseline_dv']:>14.3f} "
            f"{relaxed_dv} "
            f"{retry_dv} "
            f"{row['gate_dv']:>14.3f} "
            f"{row['baseline_runtime']:>14.2f} "
            f"{row['gate_runtime_projected']:>14.2f} "
            f"{row['actions']:>14}"
        )


def _plot(rows, output_path):
    if not rows:
        return None
    cases = sorted(set(row["case"] for row in rows))
    fig, axes = plt.subplots(
        len(cases), 2, figsize=(12, max(4, 3.2 * len(cases))), squeeze=False
    )
    for row_index, case in enumerate(cases):
        case_rows = [row for row in rows if row["case"] == case]
        quality = [
            [row["baseline_dv"] for row in case_rows],
            [row["gate_dv"] for row in case_rows],
        ]
        runtime = [
            [row["baseline_runtime"] for row in case_rows],
            [row["gate_runtime_projected"] for row in case_rows],
        ]
        for col_index, (data, title, ylabel) in enumerate((
            (quality, f"{case} quality", "Objective DV (m/s)"),
            (runtime, f"{case} projected branch runtime", "Runtime (s)"),
        )):
            ax = axes[row_index][col_index]
            box = ax.boxplot(
                data,
                tick_labels=["baseline", "L0 gate best-of-two"],
                patch_artist=True,
            )
            for patch, color in zip(box["boxes"], ["#4C78A8", "#F58518"]):
                patch.set_facecolor(color)
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
            ax.set_title(title)
            ax.set_ylabel(ylabel)
            ax.grid(True, axis="y", alpha=0.25)
    fig.tight_layout()
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=150)
    plt.close(fig)
    return output_path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--label", default="l0_pressure_gate")
    parser.add_argument("--output-dir", default=str(ROOT / "Tests"))
    parser.add_argument("--cases", default="vanilla_kekkj,jnsq_keemo")
    parser.add_argument(
        "--strategy", default="funnel_l0_recall_64_mbh_between",
        help="Baseline optimizer strategy to gate.",
    )
    parser.add_argument("--seeds", default="0,1,2,3,4")
    parser.add_argument("--opt-level", default="benchmark_funnel")
    parser.add_argument("--source", default="stage_1_final")
    parser.add_argument("--top", type=int, default=32)
    parser.add_argument("--near-fraction", type=float, default=0.03)
    parser.add_argument("--min-pressure-count", type=int, default=2)
    parser.add_argument("--relax-fraction", type=float, default=0.20)
    parser.add_argument("--min-relax-days", type=float, default=5.0)
    parser.add_argument("--max-relax-days", type=float, default=60.0)
    parser.add_argument(
        "--adjustment-mode",
        choices=("widen", "shift"),
        default="widen",
        help=(
            "widen expands the pressured flank. shift moves a leg window when "
            "only one flank is pressured, and widens when both flanks are."
        ),
    )
    parser.add_argument(
        "--branch-policy",
        choices=("pressure", "quality_guard", "l1_improvement", "l1_combined"),
        default="pressure",
        help=(
            "pressure branches on any L0 pressure. quality_guard requires "
            "pressure plus a per-sequence objective threshold. "
            "l1_improvement requires pressure plus poor L0->L1 improvement. "
            "l1_combined also branches when L1 best remains high relative "
            "to the L0 top population."
        ),
    )
    parser.add_argument(
        "--max-l1-improvement",
        type=float,
        default=0.15,
        help=(
            "Branch with --branch-policy l1_improvement when "
            "(L0_best-L1_best)/L0_best is at or below this value."
        ),
    )
    parser.add_argument(
        "--min-l1-best-to-l0-median",
        type=float,
        default=0.725,
        help=(
            "With --branch-policy l1_combined, branch when "
            "L1_best / median(L0 top population) is at or above this value."
        ),
    )
    parser.add_argument(
        "--quality-thresholds",
        default="",
        help="Comma-separated SEQUENCE:DV thresholds, e.g. KEEMo:5500,KEKKJ:1500.",
    )
    parser.add_argument(
        "--rescue-mode",
        choices=("same_seed", "retry_seed", "both", "cascade"),
        default="same_seed",
        help=(
            "same_seed runs the adjusted-bounds branch with the original seed. "
            "retry_seed reruns adjusted bounds with seed+offset. both keeps "
            "the better of both rescue branches. cascade runs same_seed first, "
            "then retry_seed only if the same-seed result remains suspicious."
        ),
    )
    parser.add_argument(
        "--cascade-min-improvement",
        type=float,
        default=0.10,
        help=(
            "Fallback suspicion guard for --rescue-mode cascade when no "
            "quality threshold exists: retry if the same-seed branch improves "
            "the current best by at most this fraction."
        ),
    )
    parser.add_argument(
        "--retry-seed-offset",
        type=int,
        default=100,
        help="Seed offset used by --rescue-mode retry_seed/both.",
    )
    parser.add_argument("--no-auto-workers", action="store_true")
    parser.add_argument(
        "--summary-only", action="store_true",
        help="Read existing DBs only; do not run missing optimizer jobs.",
    )
    parser.add_argument("--plot", default=None)
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    case_names = parse_csv(args.cases)
    seeds = [int(value) for value in parse_csv(args.seeds)]
    quality_thresholds = _parse_thresholds(args.quality_thresholds)

    for case_name in case_names:
        if case_name not in CASES:
            raise SystemExit(f"Unknown case {case_name!r}. Valid: {sorted(CASES)}")

    baseline_db = db_for_strategy(output_dir, args.label, args.strategy)
    pressure_db = output_dir / (
        baseline_db.stem + "_pressure_gate.sqlite"
    )

    if not args.summary_only:
        for case_name in case_names:
            for seed in seeds:
                print(
                    f"Running baseline {case_name} seed={seed} "
                    f"strategy={args.strategy} db={baseline_db.name}",
                    flush=True,
                )
                prepare_and_run(
                    baseline_db,
                    case_name,
                    CASES[case_name],
                    args.strategy,
                    args.opt_level,
                    seed,
                    auto_workers=not args.no_auto_workers,
                    batch_name=_batch_name(case_name, args.strategy, seed),
                )

    rows = []
    for case_name in case_names:
        case = CASES[case_name]
        sequence = _sequence_short_name(case)
        for seed in seeds:
            baseline_batch = _batch_name(case_name, args.strategy, seed)
            baseline = _result_for_batch(
                baseline_db, baseline_batch, sequence, seed,
            )
            if baseline is None:
                print(
                    f"Missing baseline result for {case_name} seed={seed}; skipped",
                    flush=True,
                )
                continue
            bounds, tof_vectors = _population_tofs(
                baseline_db, baseline, args.source, args.top,
            )
            actions = pressure_actions(
                bounds,
                tof_vectors,
                near_fraction=args.near_fraction,
                min_pressure_count=args.min_pressure_count,
            )
            relaxed_bounds, deltas, adjustment_modes = adjust_bounds(
                bounds,
                actions,
                relax_fraction=args.relax_fraction,
                min_relax_days=args.min_relax_days,
                max_relax_days=args.max_relax_days,
                adjustment_mode=args.adjustment_mode,
            )
            baseline_runtime = float(baseline["runtime_seconds"])
            baseline_l0_runtime = _stage_runtime(
                baseline_db, baseline["run_id"], stage_index=1,
            )
            l0_best = _stage_best(baseline_db, baseline["run_id"], stage_index=1)
            l1_best = _stage_best(baseline_db, baseline["run_id"], stage_index=2)
            l1_improvement = (
                (l0_best - l1_best) / l0_best
                if l0_best not in (None, 0.0) and l1_best is not None
                else None
            )
            l0_fitnesses = _population_fitnesses(
                baseline_db, baseline["run_id"], args.source, args.top,
            )
            l0_median = statistics.median(l0_fitnesses) if l0_fitnesses else None
            l1_best_to_l0_median = (
                l1_best / l0_median
                if l1_best is not None and l0_median not in (None, 0.0)
                else None
            )
            selected = baseline
            relaxed_runtime = 0.0
            relaxed_dv = None
            retry_runtime = 0.0
            retry_dv = None
            should_branch = _should_branch(
                args.branch_policy,
                sequence,
                float(baseline["objective_dv"]),
                actions,
                quality_thresholds,
                l1_improvement=l1_improvement,
                max_l1_improvement=args.max_l1_improvement,
                l1_best_to_l0_median=l1_best_to_l0_median,
                min_l1_best_to_l0_median=args.min_l1_best_to_l0_median,
            )
            if should_branch:
                if args.rescue_mode in ("same_seed", "both", "cascade"):
                    pressure_batch = _pressure_batch_name(
                        case_name,
                        args.strategy,
                        seed,
                        args.top,
                        args.near_fraction,
                        args.relax_fraction,
                        args.adjustment_mode,
                    )
                    if not args.summary_only:
                        print(
                            f"Pressure {case_name} seed={seed}: "
                            f"{_actions_label(actions, deltas, adjustment_modes)}; "
                            f"running adjusted same-seed db={pressure_db.name}",
                            flush=True,
                        )
                        prepare_and_run(
                            pressure_db,
                            case_name,
                            case,
                            args.strategy,
                            args.opt_level,
                            seed,
                            auto_workers=not args.no_auto_workers,
                            batch_name=pressure_batch,
                            leg_tof_bounds_override=relaxed_bounds,
                        )
                    relaxed = _result_for_batch(
                        pressure_db, pressure_batch, sequence, seed,
                    )
                    if relaxed is None and args.adjustment_mode == "widen":
                        legacy_pressure_batch = _legacy_pressure_batch_name(
                            case_name,
                            args.strategy,
                            seed,
                            args.top,
                            args.near_fraction,
                            args.relax_fraction,
                        )
                        relaxed = _result_for_batch(
                            pressure_db, legacy_pressure_batch, sequence, seed,
                        )
                    if relaxed is not None:
                        relaxed_runtime = float(relaxed["runtime_seconds"])
                        relaxed_dv = float(relaxed["objective_dv"])
                        if relaxed_dv < float(selected["objective_dv"]):
                            selected = relaxed
                    else:
                        print(
                            f"Missing relaxed result for {case_name} seed={seed}; "
                            "using current best in summary",
                            flush=True,
                        )

                run_retry = args.rescue_mode in ("retry_seed", "both")
                if args.rescue_mode == "cascade":
                    run_retry = _should_retry_after_same_seed(
                        sequence,
                        float(baseline["objective_dv"]),
                        float(selected["objective_dv"]),
                        relaxed_dv,
                        quality_thresholds,
                        args.cascade_min_improvement,
                    )

                if run_retry:
                    retry_seed = int(seed) + int(args.retry_seed_offset)
                    retry_batch = _retry_pressure_batch_name(
                        case_name,
                        args.strategy,
                        seed,
                        args.top,
                        args.near_fraction,
                        args.relax_fraction,
                        args.adjustment_mode,
                        args.retry_seed_offset,
                    )
                    if not args.summary_only:
                        cascade_note = (
                            " after same-seed still suspect"
                            if args.rescue_mode == "cascade" else ""
                        )
                        print(
                            f"Pressure {case_name} seed={seed}: "
                            f"{_actions_label(actions, deltas, adjustment_modes)}; "
                            f"running retry seed={retry_seed}{cascade_note} "
                            f"db={pressure_db.name}",
                            flush=True,
                        )
                        prepare_and_run(
                            pressure_db,
                            case_name,
                            case,
                            args.strategy,
                            args.opt_level,
                            retry_seed,
                            auto_workers=not args.no_auto_workers,
                            batch_name=retry_batch,
                            leg_tof_bounds_override=relaxed_bounds,
                        )
                    retry = _result_for_batch(
                        pressure_db, retry_batch, sequence, retry_seed,
                    )
                    if retry is not None:
                        retry_runtime = float(retry["runtime_seconds"])
                        retry_dv = float(retry["objective_dv"])
                        if retry_dv < float(selected["objective_dv"]):
                            selected = retry
                    else:
                        print(
                            f"Missing retry result for {case_name} seed={seed}; "
                            "using current best in summary",
                            flush=True,
                        )

            l0_plus_rescue_runtime_projected = (
                baseline_l0_runtime + relaxed_runtime + retry_runtime
                if should_branch and (relaxed_runtime > 0.0 or retry_runtime > 0.0)
                else baseline_runtime
            )
            gate_runtime_projected = (
                baseline_runtime + relaxed_runtime + retry_runtime
                if should_branch and (relaxed_runtime > 0.0 or retry_runtime > 0.0)
                else baseline_runtime
            )
            rows.append({
                "case": sequence,
                "seed": seed,
                "pressure_detected": bool(actions),
                "triggered": bool(should_branch),
                "l1_improvement": (
                    float(l1_improvement) if l1_improvement is not None else float("nan")
                ),
                "l1_best_to_l0_median": (
                    float(l1_best_to_l0_median)
                    if l1_best_to_l0_median is not None else float("nan")
                ),
                "baseline_dv": float(baseline["objective_dv"]),
                "relaxed_dv": relaxed_dv,
                "retry_dv": retry_dv,
                "gate_dv": float(selected["objective_dv"]),
                "baseline_runtime": baseline_runtime,
                "baseline_l0_runtime": baseline_l0_runtime,
                "l0_plus_rescue_runtime_projected": l0_plus_rescue_runtime_projected,
                "gate_runtime_projected": gate_runtime_projected,
                "actions": _actions_label(actions, deltas, adjustment_modes),
            })

    _print_summary(rows)
    if rows:
        pressure_detected = sum(1 for row in rows if row["pressure_detected"])
        triggered = sum(1 for row in rows if row["triggered"])
        print(
            f"Pressure detected {pressure_detected}/{len(rows)} runs; "
            f"branched {triggered}/{len(rows)} runs; "
            f"median baseline DV={statistics.median(row['baseline_dv'] for row in rows):.3f}; "
            f"median gate DV={statistics.median(row['gate_dv'] for row in rows):.3f}"
        )
        print(
            "Median projected runtime: "
            f"baseline={statistics.median(row['baseline_runtime'] for row in rows):.2f}s, "
            "L0+rescue="
            f"{statistics.median(row['l0_plus_rescue_runtime_projected'] for row in rows):.2f}s, "
            "full-gate="
            f"{statistics.median(row['gate_runtime_projected'] for row in rows):.2f}s"
        )
    if args.plot:
        output_path = _plot(rows, args.plot)
        if output_path is not None:
            print(f"Saved plot: {output_path}")


if __name__ == "__main__":
    main()
