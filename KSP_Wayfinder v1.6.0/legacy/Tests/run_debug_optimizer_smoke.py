# -*- coding: utf-8 -*-
"""Tiny optimizer smoke test.

This intentionally uses the debug optimization level. It only checks that a
Wayfinder job can pass through pygmo and the standalone trajectory helpers.
"""

import sys

sys.path.append("../WayfinderCore")

from _Wayfinder import Wayfinder  # noqa: E402


DATASTORE_NAME = "Debug_Optimizer_Smoke"
SWING_BY_BODIES = [["Kerbin"], ["Eve"], ["Moho"]]


def run_smoke_test():
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
    plans.optimize(n=1)

    job_index = plans._df.index[0]
    assert plans._df.at[job_index, "job_status"] == "DONE"
    assert len(plans._df.at[job_index, "gene"]) > 0
    print(plans._df.loc[job_index, ["job_status", "result_DV", "result_t0", "result_tof", "result_ej_vinf"]])


if __name__ == "__main__":
    run_smoke_test()
