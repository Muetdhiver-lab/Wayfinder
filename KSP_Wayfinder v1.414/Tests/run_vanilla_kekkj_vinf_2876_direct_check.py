# -*- coding: utf-8 -*-
"""Check whether the KEKKJ 2876 m/s alpha solution survives direct encoding."""

import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "WayfinderCore"))

from _SQLiteStore import SQLiteJobStore  # noqa: E402
from _Wayfinder import Wayfinder  # noqa: E402


DB_PATH = ROOT / "Tests" / "vanilla_kekkj_vinf_bin_check.sqlite"
RUN_ID = 2


def main():
    plans = Wayfinder(planet_pack="Vanilla")
    store = SQLiteJobStore(DB_PATH)
    try:
        context = store.run_context(RUN_ID)
    finally:
        store.close()

    print("source:")
    print(
        json.dumps(
            {
                "run_id": RUN_ID,
                "objective_dv": context["objective_dv"],
                "result_t0": context["result_t0"],
                "result_tof": context["result_tof"],
                "tof_encoding": context["tof_encoding"],
            },
            indent=2,
        )
    )

    result = plans.reoptimize_direct_near_run_sqlite(
        DB_PATH,
        run_id=RUN_ID,
        t0_wiggle_days=10,
        leg_tof_wiggle_days=10,
        sade_gen=80,
        pop_size=80,
    )
    print("direct check:")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
