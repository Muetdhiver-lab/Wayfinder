# -*- coding: utf-8 -*-
"""Validate the PyKEP -> KSP flight-plan layer on a short K-E-D chain.

The purpose of this script is deliberately modest: find or load a Vanilla
Kerbin-Eve-Duna flyby solution, build the KSP-oriented flight plan, and print
the pilot-facing targets needed for a manual in-game check.

K-E-D is a good translation validation case because:

* it has only one flyby B-plane to hit;
* it avoids a Kerbin flyby in the middle of the chain;
* Duna can be selected as the final target in KSP;
* the arcs are short enough that KSP patched-conic drift is less pathological
  than in KEKKJ-style chains.

This is not an autopilot or correction harness. kRPC injection, if desired, is
handled by ``Tests/inject_ksp_plan.py`` after this script identifies a run_id.
"""

import argparse
import json
import sqlite3
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CORE_DIR = ROOT / "WayfinderCore"
sys.path.insert(0, str(CORE_DIR))

from _KSPBridge import build_ksp_flight_plan  # noqa: E402
from _KSPBridge import format_ksp_flight_plan  # noqa: E402
from _Wayfinder import Wayfinder  # noqa: E402


DEFAULT_DB = ROOT / "Tests" / "vanilla_ked_translation_validation.sqlite"
DEFAULT_BATCH = "vanilla_ked_translation_validation"


def parse_sequence(value):
    """Parse a compact sequence such as ``KED`` or ``Kerbin,Eve,Duna``."""

    aliases = {
        "K": "Kerbin",
        "E": "Eve",
        "D": "Duna",
        "J": "Jool",
        "MO": "Moho",
        "M": "Moho",
        "EL": "Eeloo",
    }
    if "," in value or ">" in value:
        tokens = [part.strip() for part in value.replace(">", ",").split(",")]
    else:
        upper = value.upper()
        tokens = []
        i = 0
        while i < len(upper):
            two = upper[i : i + 2]
            if two in aliases:
                tokens.append(two)
                i += 2
            else:
                tokens.append(upper[i])
                i += 1
    bodies = []
    for token in tokens:
        if not token:
            continue
        body = aliases.get(token.upper(), token)
        bodies.append([body])
    if len(bodies) < 2:
        raise ValueError(f"sequence must contain at least two bodies, got {value!r}")
    return bodies


def sequence_label(sequence):
    return "".join(step[0][0] for step in sequence)


def best_run(db_path, batch_name):
    with sqlite3.connect(db_path) as con:
        con.row_factory = sqlite3.Row
        row = con.execute(
            """
            select
                r.id as run_id,
                res.objective_dv,
                res.result_t0,
                res.result_tof,
                res.ejection_vinf,
                res.arrival_vinf,
                r.runtime_seconds,
                j.*,
                s.bodies_json,
                s.short_name as sequence_short_name
            from results res
            join runs r on r.id=res.run_id
            join jobs j on j.id=r.job_id
            join sequences s on s.id=j.sequence_id
            join batch_jobs bj on bj.job_id=j.id
            join batches b on b.id=bj.batch_id
            where b.name=?
            order by res.objective_dv asc
            limit 1
            """,
            (batch_name,),
        ).fetchone()
        if row is None:
            return None
        gene = json.loads(
            con.execute(
                "select gene_json from genes where run_id=?", (row["run_id"],)
            ).fetchone()[0]
        )
        return dict(row), gene


def print_plan_summary(plan, context, result):
    print("=== KED translation validation candidate ===")
    print(f"run_id:                 {result['run_id']}")
    print(f"sequence:               {result['sequence_short_name']}")
    print(f"objective DV:           {result['objective_dv']:.3f} m/s")
    print(f"runtime:                {result['runtime_seconds']:.3f} s")
    print(f"T0 KUT:                 {result['result_t0']:.4f} days")
    print(f"TOF:                    {result['result_tof']:.4f} days")
    print(f"arrival mode:           {context['arrival_mode']}")
    print(f"operational node DV:    {plan.operational_total_dv_without_arrival:.3f} m/s")
    print()

    print("Nodes:")
    for index, node in enumerate(plan.nodes, start=1):
        print(
            "  {idx}. {label:<24} KUT {epoch:9.4f}  "
            "T+{met:9.4f}  P {p:+9.3f}  N {n:+9.3f}  R {r:+9.3f}  "
            "|dV| {dv:8.3f}  [{model}]".format(
                idx=index,
                label=node.label,
                epoch=node.epoch_days,
                met=node.met_days,
                p=node.prograde,
                n=node.normal,
                r=node.radial,
                dv=node.magnitude,
                model=node.model,
            )
        )
    print()

    print("Flyby targets:")
    for flyby in plan.flybys:
        print(
            "  {body} leg {leg}: Pe alt {pe:.3f} km, "
            "hyperbola inc {inc:.3f} deg, beta {beta:.3f} deg".format(
                body=flyby.body,
                leg=flyby.leg_index,
                pe=flyby.periapsis_altitude_km,
                inc=flyby.hyperbola_inclination_deg,
                beta=flyby.beta_angle_deg,
            )
        )
        print(
            "    SOI entry KUT {entry:.4f}, Pe KUT {pe_epoch:.4f}, "
            "SOI exit KUT {exit:.4f}".format(
                entry=flyby.soi_entry_epoch_days,
                pe_epoch=flyby.epoch_days,
                exit=flyby.soi_exit_epoch_days,
            )
        )
        print(
            "    plane normal {}, periapsis direction {}".format(
                [round(v, 5) for v in flyby.flyby_plane_normal],
                [round(v, 5) for v in flyby.periapsis_direction],
            )
        )

    if plan.arrival is not None:
        print()
        print("Arrival:")
        print(
            "  {body}: KUT {epoch:.4f}, T+{met:.4f}, v_inf {vinf:.3f} m/s, "
            "mode {mode}".format(
                body=plan.arrival.body,
                epoch=plan.arrival.epoch_days,
                met=plan.arrival.met_days,
                vinf=plan.arrival.arrival_vinf,
                mode=plan.arrival.arrival_mode,
            )
        )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--db", default=str(DEFAULT_DB))
    parser.add_argument("--batch-name", default=DEFAULT_BATCH)
    parser.add_argument("--optimize", action="store_true")
    parser.add_argument("--n", type=int, default=1)
    parser.add_argument("--t0-min", type=int, default=0)
    parser.add_argument("--t0-bin", type=int, default=400)
    parser.add_argument("--n-t0-bins", type=int, default=4)
    parser.add_argument("--sequence", default="KED")
    parser.add_argument("--tof-profile", default="relaxed")
    parser.add_argument("--opt-level", default="high")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument(
        "--strategy", choices=["selected", "split_soi"], default="split_soi"
    )
    parser.add_argument(
        "--arrival-mode",
        choices=["flyby", "vinf", "circular", "elliptical"],
        default="flyby",
    )
    parser.add_argument("--dump-json", action="store_true")
    args = parser.parse_args()

    db_path = Path(args.db)
    plans = Wayfinder(planet_pack="Vanilla")
    sequence = parse_sequence(args.sequence)
    label = sequence_label(sequence)

    if args.optimize:
        plans.add_direct_t0_batch_sqlite(
            swing_by_bodies=sequence,
            db_path=db_path,
            batch_name=args.batch_name,
            t0_min=args.t0_min,
            t0_bin=args.t0_bin,
            n_t0_bins=args.n_t0_bins,
            tof_profile=args.tof_profile,
            opt_level=args.opt_level,
            arrival_mode=args.arrival_mode,
            injection_altitude=100000,
            ejection_altitude=100000,
            purpose="translation_validation",
            optimizer_topology="ring",
            optimizer_seed=args.seed,
        )
        optimized = plans.optimize_sqlite(
            db_path,
            n=args.n,
            batch_name=args.batch_name,
        )
        print(f"Optimized {optimized} {label} validation job(s)")

    loaded = best_run(db_path, args.batch_name)
    if loaded is None:
        raise SystemExit(
            "No DONE KED result found. Re-run with --optimize or check --db/--batch-name."
        )
    context, gene = loaded
    udp = plans._mga_problem_from_sqlite_context(context)
    plan = build_ksp_flight_plan(
        udp,
        gene,
        planet_pack="Vanilla",
        parking_altitude=context["ejection_altitude"],
        strategy=args.strategy,
        arrival_mode=context.get("arrival_mode", args.arrival_mode),
    )

    print(
        format_ksp_flight_plan(
            plan,
            title=f"{label} translation validation flight plan",
            metadata=context,
        )
    )
    if args.dump_json:
        print()
        print(json.dumps(plan.as_dict(), indent=2))


if __name__ == "__main__":
    main()
