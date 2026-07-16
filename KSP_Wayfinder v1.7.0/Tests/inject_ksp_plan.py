# -*- coding: utf-8 -*-
"""Inject a Wayfinder validation run into KSP via kRPC.

Example:
    python Tests/inject_ksp_plan.py --run-id 33 --strategy selected

The kRPC server must be running in KSP and the connection must be accepted in
the kRPC window when autoAcceptConnections is disabled.

This script only creates the maneuver-node sequence. It does not attempt to
automatically correct flyby B-planes or execute the mission; use the printed
flight-plan Pe altitude and hyperbola-inclination targets for manual
leg-by-leg refinement in KSP.
"""

import argparse
import json
import sqlite3
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CORE_DIR = ROOT / "WayfinderCore"
sys.path.insert(0, str(CORE_DIR))

from _KSPBridge import KspKrpcClient  # noqa: E402
from _KSPBridge import build_ksp_flight_plan  # noqa: E402
from _KSPBridge import format_ksp_flight_plan  # noqa: E402
from _Wayfinder import Wayfinder  # noqa: E402


def load_run(db_path, run_id):
    con = sqlite3.connect(db_path)
    con.row_factory = sqlite3.Row
    try:
        cur = con.cursor()
        gene = json.loads(
            cur.execute(
                "select gene_json from genes where run_id=?", (run_id,)
            ).fetchone()[0]
        )
        context = dict(
            cur.execute(
                """
                select j.*, s.bodies_json, s.short_name as sequence_short_name
                from runs r
                join jobs j on j.id=r.job_id
                join sequences s on s.id=j.sequence_id
                where r.id=?
                """,
                (run_id,),
            ).fetchone()
        )
        return context, gene
    finally:
        con.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--db", default=str(ROOT / "wayfinder_validation.sqlite"))
    parser.add_argument("--run-id", type=int, default=33)
    parser.add_argument(
        "--strategy", choices=["selected", "split_soi"], default="selected"
    )
    parser.add_argument("--ideal-dsm", action="store_true")
    parser.add_argument(
        "--arrival-mode",
        choices=["flyby", "vinf", "circular", "elliptical"],
        default=None,
        help="Override the arrival mode stored with the run.",
    )
    parser.add_argument("--clear-existing", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--dump-json", action="store_true")
    parser.add_argument(
        "--min-node-dv",
        type=float,
        default=0.1,
        help="Skip maneuver nodes below this DV magnitude in m/s. Use 0 to inject all.",
    )
    parser.add_argument("--address", default="127.0.0.1")
    parser.add_argument("--rpc-port", type=int, default=50000)
    parser.add_argument(
        "--stream-port",
        type=int,
        default=None,
        help="Optional kRPC stream port. Omit for RPC-only mode.",
    )
    args = parser.parse_args()

    context, gene = load_run(Path(args.db), args.run_id)
    plans = Wayfinder(planet_pack=context["planet_pack"])
    udp = plans._mga_problem_from_sqlite_context(context)
    plan = build_ksp_flight_plan(
        udp,
        gene,
        planet_pack=context["planet_pack"],
        parking_altitude=context["ejection_altitude"],
        strategy=args.strategy,
        corrected_dsm=not args.ideal_dsm,
        arrival_mode=args.arrival_mode or context.get("arrival_mode", "flyby"),
    )

    metadata = {
        "run_id": args.run_id,
        "sequence_short_name": context.get("sequence_short_name"),
    }
    print(format_ksp_flight_plan(plan, title="Wayfinder KSP node injection plan", metadata=metadata))
    if args.dump_json:
        print()
        print(json.dumps(plan.as_dict(), indent=2))
    if args.dry_run:
        return

    if args.min_node_dv > 0.0:
        original_count = len(plan.nodes)
        plan.nodes = [node for node in plan.nodes if node.magnitude >= args.min_node_dv]
        skipped = original_count - len(plan.nodes)
        if skipped:
            print(
                f"Skipped {skipped} maneuver node(s) below "
                f"{args.min_node_dv:.3f} m/s before injection."
            )

    client = KspKrpcClient(
        name="Wayfinder inject",
        address=args.address,
        rpc_port=args.rpc_port,
        stream_port=args.stream_port,
    )
    client.connect()
    print("Connected:", client.active_vessel_summary())
    nodes = client.inject_plan(plan, clear_existing=args.clear_existing)
    print("Injected", len(nodes), "maneuver node(s).")


if __name__ == "__main__":
    main()
