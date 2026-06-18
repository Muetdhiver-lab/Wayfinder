# -*- coding: utf-8 -*-
"""Generate corrected SQL-only JNSQ reference results.

This script intentionally runs the optimizer and may take a while. It writes a
SQLite datastore instead of overwriting the historical JNSQ examples, whose
stored times were produced through the old implicit Vanilla decode path.
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "WayfinderCore"))

from _Wayfinder import Wayfinder  # noqa: E402


DATASTORE_NAME = "JNSQ_Moho_reference_fixed"
DB_PATH = Path(__file__).resolve().parent / f"{DATASTORE_NAME}.sqlite"
SWING_BY_BODIES = [["Kerbin"], ["Eve"], ["Eve", "*"], ["Moho"]]


def generate_reference_jobs():
    plans = Wayfinder(planet_pack="JNSQ")
    plans.add_batch_sqlite(
        db_path=DB_PATH,
        batch_name=DATASTORE_NAME,
        swing_by_bodies=SWING_BY_BODIES,
        t0_min=0,
        t0_bin=200,
        n_t0_bins=2,
        tof_min=300,
        tof_bin=150,
        n_tof_bins=2,
        opt_level="moderate",
        opt_injection="vinf",
        injection_altitude=100000,
    )


def optimize_reference_jobs():
    plans = Wayfinder(planet_pack="JNSQ")
    plans.optimize_sqlite(DB_PATH, n=8, batch_name=DATASTORE_NAME)
    plans.find_best_known_plan_sqlite(
        db_path=DB_PATH,
        batch_name=DATASTORE_NAME,
        start_body="Kerbin",
        target_body="Moho",
        t0_range=[0, 1000],
    )


if __name__ == "__main__":
    generate_reference_jobs()
    optimize_reference_jobs()
