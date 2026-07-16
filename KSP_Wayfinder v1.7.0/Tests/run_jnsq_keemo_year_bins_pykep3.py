# -*- coding: utf-8 -*-
"""Run a SQL-only KEEMo JNSQ scan intended for the pykep3 environment."""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "WayfinderCore"))

from _Wayfinder import Wayfinder  # noqa: E402


DATASTORE_NAME = "JNSQ_KEEMo_5yr_high_vinf_pykep3"
DB_PATH = Path(__file__).resolve().parent / f"{DATASTORE_NAME}.sqlite"
SWING_BY_BODIES = [["Kerbin"], ["Eve"], ["Eve"], ["Moho"]]


def generate_jobs():
    plans = Wayfinder(planet_pack="JNSQ")
    plans.add_batch_sqlite(
        db_path=DB_PATH,
        batch_name=DATASTORE_NAME,
        swing_by_bodies=SWING_BY_BODIES,
        t0_min=0,
        t0_bin=365,
        n_t0_bins=5,
        auto_tof=True,
        opt_level="high",
        opt_injection="vinf",
        injection_altitude=100000,
    )


def optimize_jobs():
    plans = Wayfinder(planet_pack="JNSQ")
    plans.optimize_sqlite(DB_PATH, n=5, batch_name=DATASTORE_NAME)
    plans.find_best_known_plan_sqlite(
        db_path=DB_PATH,
        batch_name=DATASTORE_NAME,
        start_body="Kerbin",
        target_body="Moho",
        t0_range=[0, 365 * 5],
    )


if __name__ == "__main__":
    generate_jobs()
    optimize_jobs()
