# -*- coding: utf-8 -*-
"""Run a KEEMo JNSQ reference scan with one launch bin per JNSQ year."""

import sys

sys.path.append("../WayfinderCore")

from _Wayfinder import Wayfinder  # noqa: E402


DATASTORE_NAME = "JNSQ_KEEMo_5yr_high_vinf"
SWING_BY_BODIES = [["Kerbin"], ["Eve"], ["Eve"], ["Moho"]]


def generate_jobs():
    plans = Wayfinder(planet_pack="JNSQ")
    plans.add_batch(
        datastore_name=DATASTORE_NAME,
        overwrite=True,
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
    plans.load_df(datastore_name=DATASTORE_NAME)
    plans.optimize(n=5)
    plans.find_best_plan(SWING_BY_BODIES, t0_range=[0, 365 * 5])


if __name__ == "__main__":
    generate_jobs()
    optimize_jobs()
