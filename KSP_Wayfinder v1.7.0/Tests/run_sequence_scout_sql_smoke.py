# -*- coding: utf-8 -*-
"""End-to-end SQL smoke for scout -> L0 -> promoted continuation."""

import argparse
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "WayfinderCore"))

from _Wayfinder import Wayfinder  # noqa: E402


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--db", default="Tests/_tmp_sequence_scout_sql_smoke.sqlite",
    )
    parser.add_argument("--reset", action="store_true")
    args = parser.parse_args()
    db_path = Path(args.db)
    if args.reset:
        db_path.unlink(missing_ok=True)

    plans = Wayfinder(planet_pack="Vanilla")
    prepared = plans.prepare_sequence_scout_sqlite(
        db_path, "sql_scout_smoke", "Kerbin", "Jool",
    )
    print("PREPARED " + json.dumps(prepared, sort_keys=True), flush=True)
    l0_count = plans.optimize_sequence_scout_stage_sqlite(
        db_path, prepared["scout_run_id"], stage="l0", n=1,
        auto_workers=False,
    )
    print("L0_OPT {}".format(l0_count), flush=True)
    promoted = plans.promote_sequence_scout_sqlite(
        db_path, prepared["scout_run_id"], per_bin=1, allow_partial=True,
    )
    print("PROMOTED " + json.dumps(promoted, sort_keys=True), flush=True)
    funnel_count = plans.optimize_sequence_scout_stage_sqlite(
        db_path, prepared["scout_run_id"], stage="funnel", n=1,
        auto_workers=False,
    )
    print("FUNNEL_OPT {}".format(funnel_count), flush=True)
    print(
        "STATUS " + json.dumps(
            plans.sequence_scout_status_sqlite(
                db_path, prepared["scout_run_id"],
            ),
            sort_keys=True,
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
