# -*- coding: utf-8 -*-
"""Benchmark KEEMoMo seeded from prior KEEMo solutions.

The experiment treats KEEMo as a continuation prior for KEEMoMo:

    Kerbin-Eve-Eve-Moho -> Kerbin-Eve-Eve-Moho-Moho

Direct-encoded genes are compatible up to the terminal Eve->Moho leg. We copy
that prefix and append a small set of plausible Moho flyby / Moho->Moho leg
parameters, then inject those candidates into the initial funnel population.
"""

from __future__ import annotations

import argparse
import json
import sqlite3
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "WayfinderCore"))
sys.path.insert(0, str(ROOT / "Tests"))

from benchmark_analysis import optimizer_result_rows  # noqa: E402
from benchmark_analysis import plot_quality_runtime_boxplots  # noqa: E402
from benchmark_analysis import print_grouped_summary  # noqa: E402
from benchmark_analysis import print_result_table  # noqa: E402
from _Wayfinder import Wayfinder  # noqa: E402
from run_l0_recall_benchmark import CASES, db_for_strategy, parse_csv  # noqa: E402


DEFAULT_MOHO_EXTENSION_PARAMS = [
    # beta, rp/radius, eta, Moho->Moho TOF in Earth days.
    # These include simple neutral guesses and values observed in good
    # KEEMoMo basins from the exploratory benchmark.
    [0.0, 1.1, 0.5, 88.0],
    [0.0, 1.5, 0.7, 75.0],
    [0.7649, 1.0154, 0.7025, 75.8521],
    [0.4847, 1.6503, 0.7027, 75.2221],
    [2.6407, 1.0155, 0.6454, 53.6660],
]


def _top_keemo_genes(db_paths, limit):
    rows = []
    query = """
        select
            s.short_name as sequence,
            j.optimizer_seed as seed,
            r.id as run_id,
            res.objective_dv as objective_dv,
            g.gene_json as gene_json
        from results res
        join runs r on r.id = res.run_id
        join jobs j on j.id = r.job_id
        join sequences s on s.id = j.sequence_id
        join genes g on g.run_id = r.id
        where s.short_name = 'KEEMo'
        order by res.objective_dv asc
    """
    for db_path in db_paths:
        path = Path(db_path)
        if not path.exists():
            continue
        with sqlite3.connect(path) as con:
            con.row_factory = sqlite3.Row
            for row in con.execute(query):
                item = dict(row)
                item["source_db"] = str(path)
                item["gene"] = json.loads(item.pop("gene_json"))
                rows.append(item)
    rows.sort(key=lambda row: float(row["objective_dv"]))
    unique = []
    seen = set()
    for row in rows:
        key = tuple(round(float(value), 8) for value in row["gene"])
        if key in seen:
            continue
        seen.add(key)
        unique.append(row)
        if len(unique) >= int(limit):
            break
    return unique


def _target_job_context(plans, case, tmp_db):
    tmp_db = Path(tmp_db)
    if tmp_db.exists():
        tmp_db.unlink()
    plans.add_direct_t0_batch_sqlite(
        swing_by_bodies=case["swing_by_bodies"],
        db_path=tmp_db,
        batch_name="bounds_probe",
        t0_min=case["t0_min"],
        t0_bin=case["t0_bin"],
        n_t0_bins=case["n_t0_bins"],
        tof_profile=case["tof_profile"],
        opt_level="benchmark_funnel",
        arrival_mode=case["arrival_mode"],
        optimizer_seed=0,
    )
    with sqlite3.connect(tmp_db) as con:
        con.row_factory = sqlite3.Row
        row = con.execute(
            """
            select jobs.*, sequences.bodies_json
            from jobs
            join sequences on sequences.id = jobs.sequence_id
            limit 1
            """
        ).fetchone()
        return dict(row)


def _clip_gene(gene, lower_bounds, upper_bounds):
    return [
        min(max(float(value), float(low)), float(high))
        for value, low, high in zip(gene, lower_bounds, upper_bounds)
    ]


def build_keemomo_prior_genes(
    keemo_rows, plans, case, extension_params=None, tmp_db=None,
):
    extension_params = extension_params or DEFAULT_MOHO_EXTENSION_PARAMS
    tmp_db = tmp_db or (ROOT / "Tests" / "_tmp_keemomo_prior_bounds.sqlite")
    context = _target_job_context(plans, case, tmp_db)
    problem = plans._mga_problem_from_sqlite_context(context)
    lower_bounds, upper_bounds = __import__("pygmo").problem(problem).get_bounds()

    prior_genes = []
    for row in keemo_rows:
        source_gene = list(map(float, row["gene"]))
        if len(source_gene) != 14:
            continue
        for params in extension_params:
            candidate = source_gene + list(map(float, params))
            prior_genes.append(_clip_gene(candidate, lower_bounds, upper_bounds))
    return prior_genes


def _batch_name(case_name, strategy, seed):
    return f"priorseed_{case_name}_{strategy}_seed{seed}"


def run_prior_seed_benchmark(
    output_db, seeds, strategy, opt_level, prior_genes, auto_workers=True,
):
    case_name = "jnsq_keemomo"
    case = CASES[case_name]
    plans = Wayfinder(planet_pack=case["planet_pack"])
    output_db = Path(output_db)
    for seed in seeds:
        batch_name = _batch_name(case_name, strategy, seed)
        plans.add_direct_t0_batch_sqlite(
            swing_by_bodies=case["swing_by_bodies"],
            db_path=output_db,
            batch_name=batch_name,
            t0_min=case["t0_min"],
            t0_bin=case["t0_bin"],
            n_t0_bins=case["n_t0_bins"],
            tof_profile=case["tof_profile"],
            opt_level=opt_level,
            arrival_mode=case["arrival_mode"],
            optimizer_topology="ring",
            optimizer_seed=seed,
            purpose="keemomo_prior_seed_benchmark",
        )
        print(
            f"Running prior-seeded KEEMoMo seed={seed} "
            f"strategy={strategy} db={output_db.name}",
            flush=True,
        )
        plans.optimize_sqlite(
            output_db,
            n=1,
            batch_name=batch_name,
            optimizer_strategy=strategy,
            auto_workers=auto_workers,
            initial_seed_genes=prior_genes,
        )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--label", default="keemomo_prior_seed_10seeds")
    parser.add_argument("--output-dir", default=str(ROOT / "Tests"))
    parser.add_argument("--strategy", default="funnel_l0_recall_64_mbh_between")
    parser.add_argument("--opt-level", default="benchmark_funnel")
    parser.add_argument("--seeds", default="0,1,2,3,4,5,6,7,8,9")
    parser.add_argument("--prior-count", type=int, default=8)
    parser.add_argument("--no-auto-workers", action="store_true")
    parser.add_argument("--summary-only", action="store_true")
    parser.add_argument(
        "--keemo-db",
        action="append",
        default=[
            str(ROOT / "Tests" / "l1_combined_10seeds_funnel_l0_recall_64_mbh_between_pressure_gate.sqlite"),
            str(ROOT / "Tests" / "l1_combined_10seeds_funnel_l0_recall_64_mbh_between.sqlite"),
        ],
    )
    parser.add_argument(
        "--baseline-db",
        default=str(ROOT / "Tests" / "keemomo_l1_combined_10seeds_funnel_l0_recall_64_mbh_between.sqlite"),
    )
    parser.add_argument(
        "--plot",
        default=str(ROOT / "Tests" / "keemomo_prior_seed_vs_baseline_10seeds.png"),
    )
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_db = db_for_strategy(output_dir, args.label, args.strategy)
    seeds = [int(seed) for seed in parse_csv(args.seeds)]

    plans = Wayfinder(planet_pack="JNSQ")
    keemo_rows = _top_keemo_genes(args.keemo_db, args.prior_count)
    if not keemo_rows:
        raise SystemExit("No KEEMo prior genes found.")
    prior_genes = build_keemomo_prior_genes(
        keemo_rows,
        plans,
        CASES["jnsq_keemomo"],
    )
    print(
        f"Built {len(prior_genes)} KEEMoMo prior genes from "
        f"{len(keemo_rows)} KEEMo solutions.",
        flush=True,
    )

    if not args.summary_only:
        run_prior_seed_benchmark(
            output_db,
            seeds,
            args.strategy,
            args.opt_level,
            prior_genes,
            auto_workers=not args.no_auto_workers,
        )

    rows = []
    rows.extend(
        optimizer_result_rows(
            args.baseline_db, "baseline", sequence_short_name="KEEMoMo",
        )
    )
    rows.extend(
        optimizer_result_rows(
            output_db, "prior_seed", sequence_short_name="KEEMoMo",
        )
    )
    print_result_table(rows)
    print_grouped_summary(
        rows,
        threshold_by_group={"baseline": 5500.0, "prior_seed": 5500.0},
    )
    output_path = plot_quality_runtime_boxplots(
        rows,
        args.plot,
        group_order=["baseline", "prior_seed"],
        title_prefix="KEEMoMo",
    )
    if output_path is not None:
        print(f"Saved plot: {output_path}")


if __name__ == "__main__":
    main()
