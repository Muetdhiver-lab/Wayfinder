# -*- coding: utf-8 -*-
"""Run a JNSQ KEEMo optimization around T0 ~= 1600 and plot a porkchop."""

import sys
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "WayfinderCore"))

from _SQLiteStore import SQLiteJobStore  # noqa: E402
from _Wayfinder import Wayfinder  # noqa: E402


BATCH_NAME = "JNSQ_KEEMo_T1600_Porkchop_Deep"
DB_PATH = Path(__file__).resolve().parent / f"{BATCH_NAME}.sqlite"
PNG_PATH = Path(__file__).resolve().parent / f"{BATCH_NAME}.png"
OPTIMIZER_PNG_PATH = Path(__file__).resolve().parent / f"{BATCH_NAME}_optimizer.png"
SAMPLED_PNG_PATH = Path(__file__).resolve().parent / f"{BATCH_NAME}_sampled.png"


def run_test():
    if DB_PATH.exists():
        DB_PATH.unlink()
    if PNG_PATH.exists():
        PNG_PATH.unlink()
    if OPTIMIZER_PNG_PATH.exists():
        OPTIMIZER_PNG_PATH.unlink()
    if SAMPLED_PNG_PATH.exists():
        SAMPLED_PNG_PATH.unlink()

    plans = Wayfinder(planet_pack="JNSQ")
    plans._opt_levels_dic["porkchop_deep"] = [120, 32, 25, 160]
    plans.add_batch_sqlite(
        db_path=DB_PATH,
        batch_name=BATCH_NAME,
        swing_by_bodies=[["Kerbin"], ["Eve"], ["Eve"], ["Moho"]],
        t0_min=1460,
        t0_bin=365,
        n_t0_bins=1,
        auto_tof=True,
        opt_level="porkchop_deep",
        opt_injection="vinf",
    )
    started_at = datetime.now()
    loaded = plans.optimize_sqlite(DB_PATH, n=1, batch_name=BATCH_NAME)
    ended_at = datetime.now()
    assert loaded == 1

    store = SQLiteJobStore(DB_PATH)
    try:
        best = store.best_results("JNSQ", batch_name=BATCH_NAME, limit=1)[0]
        run_id = best["run_id"]
        n_snapshots = store.count_rows("optimizer_snapshots")
        n_population = store.count_rows("optimizer_population_points")
    finally:
        store.close()

    points = plans.plot_optimizer_porkchop_sqlite(
        DB_PATH,
        run_id=run_id,
        metric="result_DV",
    )
    plt.savefig(OPTIMIZER_PNG_PATH, dpi=160, bbox_inches="tight")
    plt.close()

    sampled_points = plans.sample_local_porkchop_sqlite(
        DB_PATH,
        run_id=run_id,
        sampler_name="best_gene_local_grid",
        t0_span=160,
        tof_span=260,
        n_t0=55,
        n_tof=55,
    )
    plans.plot_sampled_porkchop_sqlite(
        DB_PATH,
        run_id=run_id,
        sampler_name="best_gene_local_grid",
        metric="result_DV",
    )
    plt.savefig(SAMPLED_PNG_PATH, dpi=160, bbox_inches="tight")
    print({
        "db_path": str(DB_PATH),
        "optimizer_png_path": str(OPTIMIZER_PNG_PATH),
        "sampled_png_path": str(SAMPLED_PNG_PATH),
        "run_id": run_id,
        "n_optimizer_points": int(len(points)),
        "n_sampled_points": int(len(sampled_points)),
        "n_snapshots": int(n_snapshots),
        "n_population": int(n_population),
        "best_dv": float(best["objective_dv"]),
        "best_t0": float(best["result_t0"]),
        "best_tof": float(best["result_tof"]),
        "runtime_seconds": round((ended_at - started_at).total_seconds(), 1),
    })


if __name__ == "__main__":
    run_test()
