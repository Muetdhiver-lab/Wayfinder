# -*- coding: utf-8 -*-
"""Tiny optimizer smoke test with direct SQLite persistence."""

import os
import sys
import uuid
from pathlib import Path

sys.path.append("../WayfinderCore")

from _SQLiteStore import SQLiteJobStore  # noqa: E402
from _Wayfinder import Wayfinder  # noqa: E402


DATASTORE_NAME = "Debug_Optimizer_SQLite_Smoke"
SWING_BY_BODIES = [["Kerbin"], ["Eve"], ["Moho"]]


def run_smoke_test():
    db_path = Path.cwd() / f"{DATASTORE_NAME}_{uuid.uuid4().hex}.sqlite"
    plans = Wayfinder(planet_pack="Vanilla")
    plans.add_batch(
        datastore_name=DATASTORE_NAME,
        overwrite=True,
        swing_by_bodies=SWING_BY_BODIES,
        t0_min=0,
        t0_bin=100,
        n_t0_bins=1,
        auto_tof=True,
        opt_level="debug",
        opt_injection="vinf",
        injection_altitude=100000,
    )

    plans.load_df(datastore_name=DATASTORE_NAME)
    plans.optimize(n=1, sqlite_db_path=db_path, sqlite_batch_name=DATASTORE_NAME)

    store = SQLiteJobStore(db_path)
    try:
        assert store.count_rows("jobs") == 1
        assert store.count_rows("results") == 1
        best = store.best_results(
            planet_pack="Vanilla",
            start_body="Kerbin",
            target_body="Moho",
            limit=1,
        )
        assert len(best) == 1
        assert best[0]["gene_json"]
        print(best[0])
    finally:
        store.close()
        db_path.unlink(missing_ok=True)


if __name__ == "__main__":
    run_smoke_test()
