# -*- coding: utf-8 -*-
"""Re-optimize Vanilla KEKKJ with a true flyby arrival objective.

This benchmark intentionally does not minimize terminal Jool v-infinity. It is
meant to recover the old low-DV KEKKJ family where the arrival is a free
encounter/flyby-style endpoint rather than a rendezvous-like v_inf match.
"""

import sys
import sqlite3
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "WayfinderCore"))

from _Wayfinder import Wayfinder  # noqa: E402


DB_PATH = ROOT / "Tests" / "vanilla_kekkj_flyby_reopt.sqlite"
BATCH_NAME = "vanilla_kekkj_flyby_reopt"


def main():
    plans = Wayfinder(planet_pack="Vanilla")
    plans.add_batch_sqlite(
        db_path=DB_PATH,
        batch_name=BATCH_NAME,
        swing_by_bodies=[["Kerbin"], ["Eve"], ["Kerbin"], ["Kerbin"], ["Jool"]],
        t0_min=600,
        t0_bin=200,
        n_t0_bins=1,
        tof_min=3000,
        tof_bin=1000,
        n_tof_bins=1,
        auto_tof=False,
        opt_level="ultra",
        opt_injection="flyby",
        injection_altitude=1000000,
        ejection_altitude=100000,
        purpose="validation",
        optimizer_topology="ring",
        optimizer_seed=0,
    )
    optimized = plans.optimize_sqlite(DB_PATH, n=1, batch_name=BATCH_NAME)
    print(f"Optimized {optimized} Vanilla KEKKJ flyby job(s)")
    with sqlite3.connect(DB_PATH) as con:
        con.row_factory = sqlite3.Row
        rows = con.execute(
            """
            select
                s.short_name as sequence_short_name,
                j.t0_min,
                j.t0_max,
                j.tof_min,
                j.tof_max,
                j.arrival_mode,
                j.add_vinf_arr,
                j.orbit_insertion,
                res.objective_dv,
                res.result_t0,
                res.result_tof,
                res.ejection_vinf,
                res.arrival_vinf,
                r.runtime_seconds
            from results res
            join runs r on r.id=res.run_id
            join jobs j on j.id=r.job_id
            join sequences s on s.id=j.sequence_id
            order by res.objective_dv asc
            """
        ).fetchall()
    for row in rows:
        print(dict(row))


if __name__ == "__main__":
    main()
