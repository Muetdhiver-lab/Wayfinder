# -*- coding: utf-8 -*-
"""Compare L0 handoffs with one shared population and equal downstream budget.

The source database supplies the completed L0 populations. For every run, the
benchmark builds two sets of seeds from exactly the same L0:

* ``phase_farthest`` reproduces the current champion/diversity handoff;
* ``hill_valley_l0`` applies the strict Hill-Valley test to the full L0.

Both seed sets are then sent through the same L1/MBH/L2/L3 stage plan with the
same PyGMO seed. Selection overhead is timed separately from optimizer runtime.
"""

from __future__ import annotations

import argparse
import copy
import json
import sqlite3
import sys
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pygmo as pg


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "WayfinderCore"))

from _Optimization import WayfinderFitnessDecorator  # noqa: E402
from _OptimizationService import OptimizationService  # noqa: E402
from _SQLiteStore import SQLiteJobStore  # noqa: E402
from _Wayfinder import Wayfinder  # noqa: E402


class _MemoryTelemetryStore:
    """Minimal store used to retain convergence without mutating the source DB."""

    def __init__(self):
        self.snapshots = []
        self.stages = []

    def record_optimizer_snapshot(
        self, run_id, step, champion_fitness, champion_genes,
        populations=None, telemetry=None,
    ):
        values = [float(value[0]) for value in champion_fitness]
        population_values = [
            float(value[0])
            for population in (populations or [])
            for value in population.get_f()
        ]
        self.snapshots.append({
            "step": int(step),
            "best": min(values),
            "average": float(np.mean(population_values)),
        })

    def record_optimizer_population(self, *args, **kwargs):
        return None

    def record_optimizer_stage(self, run_id, summary):
        self.stages.append(dict(summary))


def _source_runs(db_path, cases, seeds):
    placeholders_cases = ",".join("?" for _ in cases)
    placeholders_seeds = ",".join("?" for _ in seeds)
    query = f"""
        SELECT r.id AS run_id, s.short_name AS case_name, j.planet_pack,
               j.optimizer_seed AS seed, res.objective_dv AS source_dv,
               r.runtime_seconds AS source_runtime
        FROM runs r
        JOIN jobs j ON j.id = r.job_id
        JOIN sequences s ON s.id = j.sequence_id
        JOIN results res ON res.run_id = r.id
        WHERE s.short_name IN ({placeholders_cases})
          AND j.optimizer_seed IN ({placeholders_seeds})
        ORDER BY s.short_name, j.optimizer_seed, r.id DESC
    """
    with sqlite3.connect(db_path) as connection:
        connection.row_factory = sqlite3.Row
        rows = connection.execute(query, [*cases, *seeds]).fetchall()
    latest = {}
    for row in rows:
        key = (row["case_name"], int(row["seed"]))
        latest.setdefault(key, dict(row))
    missing = [
        (case, seed) for case in cases for seed in seeds
        if (case, seed) not in latest
    ]
    if missing:
        raise ValueError("Missing source runs: {}".format(missing))
    return [latest[(case, seed)] for case in cases for seed in seeds]


def _selection_problem(wayfinder, context, udp):
    soi_radius_by_name = {
        name: wayfinder.soi_radius(body)
        for name, body in wayfinder._fullname_dic.items()
    }
    decorator = WayfinderFitnessDecorator(
        planet_pack=wayfinder.planet_pack,
        bodies_by_name=wayfinder._fullname_dic,
        soi_radius_by_name=soi_radius_by_name,
        ejection_altitude=context["ejection_altitude"],
        tof_encoding=context["tof_encoding"],
        ejection_model="approximate",
    )
    return (
        pg.problem(pg.decorator_problem(udp, fitness_decorator=decorator)),
        soi_radius_by_name,
    )


def _l0_inputs(store, run_id):
    points = store.optimizer_population_points(run_id, source="stage_1_final")
    if not points:
        raise ValueError("Run {} has no stage_1_final population".format(run_id))
    champions = {}
    for point in points:
        island = int(point["island_index"])
        if island not in champions or point["fitness"] < champions[island][0]:
            champions[island] = (float(point["fitness"]), point["gene"])
    champion_genes = [
        gene for _, gene in sorted(champions.values(), key=lambda item: item[0])
    ]
    return [point["gene"] for point in points], champion_genes


def _select_handoffs(
    wayfinder, context, problem, all_genes, champion_genes, count,
):
    started = time.perf_counter()
    current = wayfinder._select_phase_diverse_elites(
        context, problem, champion_genes, count,
        elite_fraction=0.75, selection_policy="phase_farthest",
    )
    current_seconds = time.perf_counter() - started
    started = time.perf_counter()
    hill_valley = wayfinder._select_phase_diverse_elites(
        context, problem, all_genes, count,
        elite_fraction=0.35, selection_policy="hill_valley_l0",
    )
    hill_valley_seconds = time.perf_counter() - started
    current_count = max(1, int(round(0.75 * count)))
    current_count = min(current_count, count)
    started = time.perf_counter()
    portfolio_current = wayfinder._select_phase_diverse_elites(
        context, problem, champion_genes, current_count,
        elite_fraction=0.75, selection_policy="phase_farthest",
    )
    portfolio_hill_valley = wayfinder._select_phase_diverse_elites(
        context, problem, all_genes, count,
        elite_fraction=0.35, selection_policy="hill_valley_l0",
    )
    selected_keys = {
        tuple(float(value) for value in gene) for gene in portfolio_current
    }
    portfolio = list(portfolio_current)
    for gene in portfolio_hill_valley:
        key = tuple(float(value) for value in gene)
        if key in selected_keys:
            continue
        portfolio.append(gene)
        selected_keys.add(key)
        if len(portfolio) == count:
            break
    if len(portfolio) != count:
        raise ValueError("Insufficient distinct Hill-Valley portfolio seeds")
    portfolio_seconds = time.perf_counter() - started
    started = time.perf_counter()
    portfolio_16_4 = list(current)
    selected_keys = {
        tuple(float(value) for value in gene) for gene in portfolio_16_4
    }
    for gene in hill_valley:
        key = tuple(float(value) for value in gene)
        if key in selected_keys:
            continue
        portfolio_16_4.append(gene)
        selected_keys.add(key)
        if len(portfolio_16_4) == count + 4:
            break
    if len(portfolio_16_4) != count + 4:
        raise ValueError("Insufficient distinct 16+4 portfolio seeds")
    portfolio_16_4_seconds = time.perf_counter() - started + (
        current_seconds + hill_valley_seconds
    )
    return {
        "current": (current, current_seconds),
        "hill_valley": (hill_valley, hill_valley_seconds),
        "portfolio_75_25": (portfolio, portfolio_seconds),
        "portfolio_16_4": (portfolio_16_4, portfolio_16_4_seconds),
    }


def _run_downstream(
    wayfinder, context, udp, soi_radius_by_name, stage_plan, seed_genes,
    seed, label,
):
    # Reproduce the public optimizer entry point and reset the global PyGMO
    # stream so branch order cannot influence either handoff.
    pg.set_global_rng_seed(int(seed))
    telemetry = _MemoryTelemetryStore()
    started = time.perf_counter()
    result = wayfinder._run_sqlite_funnel_job(
        telemetry,
        context,
        udp,
        sqlite_run_id=0,
        sade_gen=int(context["run_sade_gen"] or context["sade_gen"]),
        requested_n_island=int(
            context["requested_n_island"] or context["n_island"]
        ),
        n_island=int(context["actual_n_island"] or context["n_island"]),
        island_pop=int(context["run_island_pop"] or context["island_pop"]),
        n_evo_steps=int(context["run_n_evo_steps"] or context["n_evo_steps"]),
        topology=context["run_optimizer_topology"] or context["optimizer_topology"],
        versions={},
        adaptive_stop=None,
        soi_radius_by_name=soi_radius_by_name,
        exact_strategy="scout_archive_nm_64_g20_mbh_between",
        requested_optimizer_strategy="equal_budget_{}".format(label),
        effective_optimizer_seed=int(seed),
        initial_seed_genes=seed_genes,
        pressure_cascade=None,
        branch_label=label,
        stage_plan_override=stage_plan,
        stage_seed_offset=1,
    )
    runtime = time.perf_counter() - started
    return result, runtime, telemetry


def run_benchmark(
    source_db, cases, seeds, budget_multiplier=1.0,
    branches=("current", "hill_valley"),
):
    source_db = Path(source_db)
    source_runs = _source_runs(source_db, cases, seeds)
    store = SQLiteJobStore(source_db)
    rows = []
    try:
        for source in source_runs:
            context = store.run_context(source["run_id"])
            context["planet_pack"] = source["planet_pack"]
            wayfinder = Wayfinder(planet_pack=context["planet_pack"])
            udp = wayfinder._mga_problem_from_sqlite_context(context)
            problem, soi_radius_by_name = _selection_problem(
                wayfinder, context, udp,
            )
            full_plan = OptimizationService.funnel_stage_plan(
                int(context["actual_n_island"] or context["n_island"]),
                int(context["run_island_pop"] or context["island_pop"]),
                int(context["run_n_evo_steps"] or context["n_evo_steps"]),
                int(context["run_sade_gen"] or context["sade_gen"]),
                exact_strategy="scout_archive_nm_64_g20_mbh_between",
            )
            downstream_plan = full_plan[1:]
            # Equal-budget validation: production may stop L3 adaptively, but
            # here both handoffs must receive all ten planned refinement steps.
            for stage in downstream_plan:
                stage.pop("adaptive_stop", None)
                stage["evo_steps"] = max(
                    1, int(round(stage["evo_steps"] * float(budget_multiplier)))
                )
            handoff_count = int(downstream_plan[0]["n_island"])
            all_genes, champion_genes = _l0_inputs(store, source["run_id"])
            selections = _select_handoffs(
                wayfinder, context, problem, all_genes, champion_genes,
                handoff_count,
            )
            print(
                "{} seed {}: shared L0={} individuals, {} champions".format(
                    source["case_name"], source["seed"], len(all_genes),
                    len(champion_genes),
                ),
                flush=True,
            )
            for label in branches:
                if label not in selections:
                    raise ValueError("Unknown handoff branch: " + str(label))
                seed_genes, selection_runtime = selections[label]
                branch_plan = copy.deepcopy(downstream_plan)
                if label == "portfolio_16_4":
                    branch_plan[0].update({
                        "n_island": 20,
                        "initialization": "preselected",
                        "topology": "split_ring_16_4",
                        "split_ring": {
                            "current_islands": 16,
                            "alternative_islands": 4,
                            "bridge_weight": 0.25,
                        },
                    })
                result, downstream_runtime, telemetry = _run_downstream(
                    wayfinder, context, udp, soi_radius_by_name,
                    branch_plan, seed_genes, source["seed"], label,
                )
                row = {
                    "case": source["case_name"],
                    "seed": int(source["seed"]),
                    "source_run_id": int(source["run_id"]),
                    "source_dv": float(source["source_dv"]),
                    "branch": label,
                    "dv": float(result["result_DV"]),
                    "selection_runtime": float(selection_runtime),
                    "downstream_runtime": float(downstream_runtime),
                    "total_runtime": float(selection_runtime + downstream_runtime),
                    "actual_steps": int(result["actual_n_evo_steps"]),
                    "budget_multiplier": float(budget_multiplier),
                    "stage_summaries": result["stage_summaries"],
                    "convergence": telemetry.snapshots,
                }
                rows.append(row)
                print(
                    "  {:11s} DV={:8.1f} runtime={:6.2f}s "
                    "(selection {:5.2f}s)".format(
                        label, row["dv"], row["total_runtime"],
                        row["selection_runtime"],
                    ),
                    flush=True,
                )
    finally:
        store.close()
    return rows


def summarize(rows):
    summary = {}
    for case in sorted({row["case"] for row in rows}):
        case_rows = [row for row in rows if row["case"] == case]
        branch_summary = {}
        for branch in sorted({row["branch"] for row in case_rows}):
            selected = [row for row in case_rows if row["branch"] == branch]
            branch_summary[branch] = {
                "runs": len(selected),
                "median_dv": float(np.median([row["dv"] for row in selected])),
                "best_dv": float(min(row["dv"] for row in selected)),
                "worst_dv": float(max(row["dv"] for row in selected)),
                "median_runtime": float(np.median([
                    row["total_runtime"] for row in selected
                ])),
                "median_selection_runtime": float(np.median([
                    row["selection_runtime"] for row in selected
                ])),
            }
        summary[case] = branch_summary
    return summary


def plot(rows, output_path):
    cases = sorted({row["case"] for row in rows})
    fig, axes = plt.subplots(
        len(cases), 2, figsize=(12, 4.0 * len(cases)), squeeze=False,
    )
    for row_index, case in enumerate(cases):
        case_rows = [row for row in rows if row["case"] == case]
        seeds = sorted({row["seed"] for row in case_rows})
        lookup = {(row["seed"], row["branch"]): row for row in case_rows}
        branches = sorted({row["branch"] for row in case_rows})
        ax = axes[row_index, 0]
        if len(branches) > 1:
            for seed in seeds:
                ax.plot(
                    range(len(branches)),
                    [lookup[(seed, branch)]["dv"] for branch in branches],
                    color="#999999", alpha=0.65, marker="o",
                )
        else:
            ax.scatter(
                [0] * len(seeds),
                [lookup[(seed, branches[0])]["dv"] for seed in seeds],
                color="#4C78A8", edgecolor="black", linewidth=0.4,
            )
        ax.set_xticks(range(len(branches)), branches)
        ax.set_ylabel("Objective DV (m/s)")
        ax.set_title("{} — paired quality from identical L0".format(case))
        ax.grid(True, axis="y", alpha=0.25)

        ax = axes[row_index, 1]
        width = 0.8 / max(1, len(branches))
        x = np.arange(len(seeds))
        colors = ["#4C78A8", "#F58518", "#54A24B"]
        for index, branch in enumerate(branches):
            offset = (index - (len(branches) - 1) / 2.0) * width
            ax.bar(
                x + offset,
                [lookup[(seed, branch)]["total_runtime"] for seed in seeds],
                width, label=branch, color=colors[index % len(colors)],
            )
        ax.set_xticks(x, [str(seed) for seed in seeds])
        ax.set_xlabel("Seed")
        ax.set_ylabel("Selection + downstream runtime (s)")
        ax.set_title("{} — equal optimizer budget".format(case))
        ax.grid(True, axis="y", alpha=0.25)
        ax.legend()
    fig.tight_layout()
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=160)
    plt.close(fig)


def _parse_csv(value, cast=str):
    return [cast(item.strip()) for item in value.split(",") if item.strip()]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--source-db",
        default=str(
            ROOT / "Tests" /
            "pressure_cascade_integrated_10seeds_"
            "funnel_l0_recall_64_mbh_between_pressure_cascade.sqlite"
        ),
    )
    parser.add_argument("--cases", default="KEKKJ,KEEMo")
    parser.add_argument("--seeds", default="0,1,2,3,4,5,6,7,8,9")
    parser.add_argument("--budget-multiplier", type=float, default=1.0)
    parser.add_argument(
        "--branches", default="current,hill_valley",
        help=(
            "Comma-separated handoffs: current,hill_valley,portfolio_75_25,"
            "portfolio_16_4."
        ),
    )
    parser.add_argument(
        "--output-json",
        default=str(ROOT / "Tests" / "hill_valley_equal_budget_10seeds.json"),
    )
    parser.add_argument(
        "--output-plot",
        default=str(ROOT / "Tests" / "hill_valley_equal_budget_10seeds.png"),
    )
    args = parser.parse_args()
    rows = run_benchmark(
        args.source_db,
        _parse_csv(args.cases),
        _parse_csv(args.seeds, int),
        budget_multiplier=args.budget_multiplier,
        branches=tuple(_parse_csv(args.branches)),
    )
    payload = {"summary": summarize(rows), "rows": rows}
    output_json = Path(args.output_json)
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    plot(rows, args.output_plot)
    print(json.dumps(payload["summary"], indent=2))
    print("JSON:", output_json)
    print("Plot:", args.output_plot)


if __name__ == "__main__":
    main()
