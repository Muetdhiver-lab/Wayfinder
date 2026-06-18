# -*- coding: utf-8 -*-
"""Tiny SQL-only optimizer smoke test."""

import sys
import uuid
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "WayfinderCore"))

from _SQLiteStore import SQLiteJobStore  # noqa: E402
from _Wayfinder import Wayfinder  # noqa: E402


BATCH_NAME = "Debug_Optimizer_SQLite_Only_Smoke"
SWING_BY_BODIES = [["Kerbin"], ["Eve"], ["Moho"]]


def run_smoke_test():
    db_path = Path(__file__).resolve().parent / f"{BATCH_NAME}_{uuid.uuid4().hex}.sqlite"
    plans = Wayfinder(planet_pack="Vanilla")
    plans.add_batch_sqlite(
        db_path=db_path,
        batch_name=BATCH_NAME,
        swing_by_bodies=SWING_BY_BODIES,
        t0_min=0,
        t0_bin=100,
        n_t0_bins=1,
        auto_tof=True,
        opt_level="debug",
        opt_injection="vinf",
        injection_altitude=100000,
    )
    loaded = plans.optimize_sqlite(db_path, n=1, batch_name=BATCH_NAME)
    assert loaded == 1

    store = SQLiteJobStore(db_path)
    try:
        assert store.count_rows("jobs") == 1
        assert store.count_rows("results") == 1
        assert store.count_rows("optimizer_snapshots") == 1
        assert store.count_rows("optimizer_population_points") == 1
        best = store.best_results(
            planet_pack="Vanilla",
            start_body="Kerbin",
            target_body="Moho",
            limit=1,
        )
        assert len(best) == 1
        assert best[0]["gene_json"]
        porkchop_points = plans.porkchop_points_from_snapshots_sqlite(
            db_path,
            run_id=best[0]["run_id"],
        )
        assert len(porkchop_points) == 2
        plotted_points = plans.plot_porkchop_from_snapshots_sqlite(
            db_path,
            run_id=best[0]["run_id"],
        )
        assert len(plotted_points) == 2
        print(best[0])
    finally:
        store.close()
        db_path.unlink(missing_ok=True)


if __name__ == "__main__":
    run_smoke_test()
