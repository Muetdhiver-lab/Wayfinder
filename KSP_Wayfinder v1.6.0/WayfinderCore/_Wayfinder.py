# Launch window scanner - coarse porkchop.
"""
Created on Wed May 29 08:02:50 2019

@author: v.fave
"""
import copy
import hashlib
import json
import logging
from pathlib import Path
import subprocess
import time
import pandas as pd
import pygmo as pg
import numpy as np
import pykep
from pykep.trajopt import mga_1dsm

from planet_packs import PACKS

from _Trajectory import decode_dV_tof
from _Trajectory import decode_trajectory
from _Trajectory import transx
from _Optimization import WayfinderFitnessDecorator
from _Optimization import alpha_gene_to_direct_gene
from _Optimization import alpha_leg_tofs
from _Optimization import direct_leg_tofs
from _Optimization import direct_tof_bounds_from_leg_tofs
from _OptimizationService import OptimizationService


WAYFINDER_VERSION = "1.6.0"
DEFAULT_OPTIMIZER_STRATEGY = "funnel"


from itertools import product
from matplotlib import pyplot as plt, cm
import matplotlib.colors as colors
from matplotlib.patches import Rectangle
import seaborn as sns


logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


def _porkchop_colormap(cmap="wayfinder_lwp", over_color="red", under_color="#08306b"):
    if cmap == "wayfinder_lwp":
        colormap = colors.LinearSegmentedColormap.from_list(
            "wayfinder_lwp",
            [
                "#08306b",  # deep blue
                "#2171b5",  # blue
                "#00bcd4",  # cyan
                "#2ca25f",  # green
                "#f5e642",  # yellow
                "#fdae00",  # orange
                "#d7301f",  # red
            ],
            N=256,
        )
    else:
        colormap = plt.get_cmap(cmap).copy()
    colormap.set_under(under_color)
    colormap.set_over(over_color)
    return colormap


def _rounded_porkchop_floor(value, step=50):
    return np.ceil(float(value) / float(step)) * float(step)


def _porkchop_levels(values, mode="linear_log", floor=None, ceiling=None, floor_round=50, count=28, linear_factor=2.0):
    finite_values = np.asarray(values, dtype=float)
    finite_values = finite_values[np.isfinite(finite_values)]
    if finite_values.size == 0:
        raise ValueError("No finite values available for porkchop levels")

    if floor is None:
        floor = _rounded_porkchop_floor(np.min(finite_values), floor_round)
    floor = float(floor)
    max_value = float(ceiling) if ceiling is not None else float(np.max(finite_values))
    if max_value <= floor:
        max_value = floor * 1.01

    if mode == "linear":
        return np.linspace(floor, max_value, count)

    if mode == "log":
        return np.geomspace(max(floor, 1e-9), max_value, count)

    if mode == "linear_log":
        transition = min(floor * float(linear_factor), max_value)
        linear_count = max(6, count // 2)
        linear_levels = np.linspace(floor, transition, linear_count)
        if transition >= max_value:
            return linear_levels
        log_count = max(4, count - linear_count + 1)
        log_levels = np.geomspace(transition, max_value, log_count)
        return np.unique(np.concatenate([linear_levels, log_levels[1:]]))

    raise ValueError("Unknown porkchop level mode: " + str(mode))


def _porkchop_boundary_colormap(
    levels,
    mode="linear_log",
    transition=None,
    over_color="red",
    under_color="#08306b",
    style="linear_to_yellow",
):
    levels = np.asarray(levels, dtype=float)
    interval_count = max(1, len(levels) - 1)

    if mode != "linear_log" or transition is None or transition <= levels[0] or transition >= levels[-1]:
        colormap = _porkchop_colormap("wayfinder_lwp", over_color=over_color, under_color=under_color)
        return colormap, colors.BoundaryNorm(levels, colormap.N)

    if style == "linear_to_red":
        linear_colormap = colors.LinearSegmentedColormap.from_list(
            "wayfinder_linear_to_red",
            ["#08306b", "#2171b5", "#00bcd4", "#2ca25f", "#f5e642", "#fdae00", "#d7301f"],
        )
        log_colormap = colors.LinearSegmentedColormap.from_list(
            "wayfinder_log_dark_red",
            ["#d7301f", "#7f0000", "#2b0000"],
        )
    else:
        linear_colormap = colors.LinearSegmentedColormap.from_list(
            "wayfinder_linear_low",
            ["#08306b", "#2171b5", "#00bcd4", "#2ca25f", "#f5e642"],
        )
        log_colormap = colors.LinearSegmentedColormap.from_list(
            "wayfinder_log_high",
            ["#f5e642", "#fdae00", "#d7301f"],
        )
    centers = (levels[:-1] + levels[1:]) / 2.0
    color_list = []
    for center in centers:
        if center <= transition:
            fraction = (center - levels[0]) / max(transition - levels[0], 1e-9)
            color_list.append(linear_colormap(np.clip(fraction, 0.0, 1.0)))
        else:
            fraction = (np.log(center) - np.log(transition)) / max(np.log(levels[-1]) - np.log(transition), 1e-9)
            color_list.append(log_colormap(np.clip(fraction, 0.0, 1.0)))
    colormap = colors.ListedColormap(color_list, name="wayfinder_linear_log")
    colormap.set_under(under_color)
    colormap.set_over(over_color)
    return colormap, colors.BoundaryNorm(levels, colormap.N)


def _make_mga_1dsm(**kwargs):
    """Build a pykep3 mga_1dsm UDP."""
    kwargs = dict(kwargs)
    kwargs.pop("max_revs", None)
    return mga_1dsm(**kwargs)
class Wayfinder:
    '''
    This class will allow to search/scan one or several gravity assists sequences
    for one or several T0 and/or ToF bins.
    
    Jobs and results are stored in SQLite. Optimisation reads job dictionaries
    from SQLite directly; pandas is only used for analysis and plotting views.
    
    The main functions are : add_mga_scans => allow to add a bunch of fly_by_sequences to scan.
    The fly_by_sequences are defined with a combinatorial method, including wildcard to allow skipping a step,
    e.g. ; [["Kerbin"],["Eve"],["Kerbin","Eve,"*"],["Moho"]] means from Kerbin to Moho, passing by Eve as a first
    step, then Either Kerbin, Eve or ignore the second step.
    
    Use add_batch_sqlite / optimize_sqlite for the main SQL-only workflow.
    
    Note: optimize can target a SQLite batch, or process TODO jobs in the selected datastore.
    
    Be aware that using the 'high' optimizsation level can get very time consuming, especially in large batches.
    
    
    TODO / IDEA list :
        Add a grid plot of sequence vs t0 with DV as color ?
        quipu style plot of transfers with color coded DV and planets passed as pearls, length between knots for the time ?
        /!\ Ejection altitude as a job parameter and patched in the fitness eval via the decorated problem.
        => Full EJ DV calc is too costly, for now simplified EJ_DV calc (underestimate the ejection inclination+ward on i at SOI to ~25°)
        => A trick would be to pre-calc a map from iSOI+Vinf => iEJ (estimate)
        => Will do for now, but will need to revisit.
        => A neat trick would be to toggle the full EJDV calc only when already deep in the optimisation. but how ? can't see a way.
        Add tests => done
        Maybe add a smart toggle in the decorated problem to trigger the use of the full ej calc when the solution gets "good enough"
        
        
    '''

    


    DEFAULT_OPTIMIZER_STRATEGY = DEFAULT_OPTIMIZER_STRATEGY

    def __init__(self,datastore_name = 'WayFinder_KSP',planet_pack = "Vanilla"):
        if planet_pack == "Vanilla" and datastore_name in PACKS:
            planet_pack = datastore_name
            datastore_name = None
        
        self.planet_pack = planet_pack
        if self.planet_pack not in PACKS:
            raise ValueError("Unknown planet pack: " + str(self.planet_pack))
        self._planet_pack_module = PACKS[self.planet_pack]
        self._Edy2Kdy = self._planet_pack_module.EDY_TO_KDY
        self._datastore_name = self._planet_pack_module.DEFAULT_DATASTORE_NAME
        self._Body_abrev_dic = self._planet_pack_module.BODY_ABBREVIATIONS
        self._fullname_dic = self._planet_pack_module.BODIES

        self._opt_levels_dic = {
                'debug'      : [1,1,7,1],             # minimum valid SADE population
                'low'        : [100,7,25,140],         # 1   UT
                'moderate'   : [100,16,25,140],       # 3   UT
                'high'       : [120,24,25,160],       # 4.5 UT   
                'wide'       : [100,48,50,140],       # 6   UT
                'ultra'      : [150,48,25,200],       # 11  UT
                'ultra+'     : [150,64,35,200],       # 19  UT         
                'deep'       : [150,24,25,200],       # 5.5 UT
                'benchmark'  : [20,8,32,5],           # ~25k nominal evaluations
                'benchmark_funnel': [20,16,32,5],     # 16 -> 8 -> 8 reference funnel
                'qualification_funnel': [50,16,32,20],# ~792k eval 16 -> 8 -> 8
                'qualification_equal': [50,16,32,25], # 25/25 stages, ~836k eval
                'benchmark_adaptive': [20,8,32,20],   # max ~103k, adaptive plateau stop
                }
        
        self._opt_insertion_dic = {
                # "add_vinf_arr", "mga_orbit_insertion", "mga_e_target", "mga_alt_target" 
                'circular'          : [True, True ,0.0],   #
                'elliptical'        : [True, True ,0.9],   #
                'vinf'              : [True ,False,0.0],   #
                'none'              : [False,False,0.0],   #
                'flyby'             : [False,False,0.0],   # arrival v_inf is reported, not penalized
                }

    @staticmethod
    def automatic_worker_count(reserve_cores=2):
        """Compatibility facade for :class:`OptimizationService`."""
        return OptimizationService.automatic_worker_count(reserve_cores)

    @classmethod
    def _resolve_island_count(cls, requested_islands, auto_workers=True, reserve_cores=2):
        return OptimizationService.resolve_island_count(
            requested_islands,
            auto_workers=auto_workers,
            reserve_cores=reserve_cores,
        )

    @staticmethod
    def _make_archipelago_topology(name, n_islands):
        return OptimizationService.make_archipelago_topology(name, n_islands)

    @staticmethod
    def _population_id_origins(populations):
        return OptimizationService.population_id_origins(populations)

    @staticmethod
    def _archipelago_telemetry(
        archipelago, populations, previous_id_origins, topology_name,
    ):
        return OptimizationService.archipelago_telemetry(
            archipelago,
            populations,
            previous_id_origins,
            topology_name,
        )

    @staticmethod
    def _balanced_annealing_schedule(problem, population_size, sade_gen):
        return OptimizationService.balanced_annealing_schedule(
            problem, population_size, sade_gen,
        )

    @staticmethod
    def _funnel_stage_plan(
        n_islands, island_pop, evo_steps, sade_gen, exact_strategy="legacy",
    ):
        return OptimizationService.funnel_stage_plan(
            n_islands, island_pop, evo_steps, sade_gen,
            exact_strategy=exact_strategy,
        )
    def _encounter_phase_embedding(self, context, problem, gene):
        return OptimizationService.encounter_phase_embedding(
            self._fullname_dic, context, problem, gene,
        )

    def _select_phase_diverse_elites(
        self, context, problem, genes, count, elite_fraction=0.35,
    ):
        return OptimizationService.select_phase_diverse_elites(
            problem,
            genes,
            count,
            embedding_for_gene=lambda gene: self._encounter_phase_embedding(
                context, problem, gene,
            ),
            elite_fraction=elite_fraction,
        )

    @staticmethod
    def _select_exact_diverse_seeds(problem, genes, count):
        return OptimizationService.select_exact_diverse_seeds(
            problem, genes, count,
        )

    def _update_exact_phase_archive(
        self, context, exact_problem, archive_genes, candidate_genes,
        max_size=32,
    ):
        return OptimizationService.update_exact_phase_archive(
            archive_genes,
            candidate_genes,
            max_size,
            select_elites=lambda genes, count: self._select_phase_diverse_elites(
                context,
                exact_problem,
                genes,
                count,
                elite_fraction=0.5,
            ),
        )

    def _planet_pack_hash(self):
        """Hash the orbital constants used by this Wayfinder instance."""
        def json_float(value):
            if value is None:
                return None
            return float(value)

        bodies = {}
        for name, body in sorted(self._fullname_dic.items()):
            bodies[name] = {
                "name": getattr(body, "name", name),
                "mu_self": json_float(getattr(body, "mu_self", None)),
                "mu_central_body": json_float(
                    getattr(body, "mu_central_body", None)
                ),
                "radius": json_float(getattr(body, "radius", None)),
                "safe_radius": json_float(getattr(body, "safe_radius", None)),
                "orbital_elements": [
                    json_float(value)
                    for value in getattr(body, "orbital_elements", [])
                ],
                "soi_radius": json_float(self.soi_radius(body)),
            }
        payload = {
            "planet_pack": self.planet_pack,
            "edy_to_kdy": self._Edy2Kdy,
            "bodies": bodies,
        }
        return hashlib.sha256(
            json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
        ).hexdigest()

    @staticmethod
    def _code_revision():
        """Return the current Git revision when the checkout is available."""
        try:
            root = Path(__file__).resolve().parents[2]
            completed = subprocess.run(
                ["git", "rev-parse", "HEAD"],
                cwd=str(root),
                capture_output=True,
                text=True,
                check=True,
            )
            return completed.stdout.strip()
        except Exception:
            return None
        
        
        
    def add_batch(self,swing_by_bodies,
                 t0_min=0  , t0_bin=100 ,  n_t0_bins=10,
                 tof_min=0 , tof_bin=1500, n_tof_bins=1,
                 opt_level = 'debug', datastore_name = 'default', 
                 auto_tof = False, 
                 overwrite = False,
                 opt_injection = 'circular',
                 injection_altitude = 100000, 
                 ejection_altitude = 100000,
                 lambert_max_revs = 0,
                 db_path = None,
                 batch_name = "default",
                 purpose="production",
                 optimizer_topology="ring",
                 optimizer_seed=None) :
        """Create a batch in SQLite.

        The old file-backed implementation has been removed. ``overwrite`` is
        kept in the signature for older scripts, but job identity is now handled
        by the SQLite parameter hash.
        """
        if db_path is None:
            raise ValueError("add_batch is SQL-only: provide db_path or call add_batch_sqlite(...).")
        if batch_name == "default":
            batch_name = datastore_name if datastore_name != "default" else self._datastore_name
        return self.add_batch_sqlite(
            swing_by_bodies=swing_by_bodies,
            db_path=db_path,
            batch_name=batch_name,
            t0_min=t0_min,
            t0_bin=t0_bin,
            n_t0_bins=n_t0_bins,
            tof_min=tof_min,
            tof_bin=tof_bin,
            n_tof_bins=n_tof_bins,
            opt_level=opt_level,
            auto_tof=auto_tof,
            opt_injection=opt_injection,
            injection_altitude=injection_altitude,
            ejection_altitude=ejection_altitude,
            lambert_max_revs=lambert_max_revs,
            purpose=purpose,
            optimizer_topology=optimizer_topology,
            optimizer_seed=optimizer_seed,
        )
    def add_batch_sqlite(self,swing_by_bodies,
                 db_path,
                 batch_name="default",
                 t0_min=0  , t0_bin=100 ,  n_t0_bins=10,
                 tof_min=0 , tof_bin=1500, n_tof_bins=1,
                 opt_level = 'debug',
                 auto_tof = False,
                 opt_injection = 'circular',
                 injection_altitude = 100000,
                 ejection_altitude = 100000,
                 lambert_max_revs = 0,
                 purpose="production",
                 optimizer_topology="ring",
                 optimizer_seed=None):
        """Create jobs directly in SQLite."""
        from _SQLiteStore import SQLiteJobStore

        if batch_name == "default":
            batch_name = self._datastore_name

        T0s  = range(t0_min,t0_min+n_t0_bins*t0_bin,t0_bin)
        ToFs = range(tof_min,tof_min+n_tof_bins*tof_bin,tof_bin)

        sequences = self.generateSequences(swing_by_bodies)
        shortSequences = self.generateShortSequences(swing_by_bodies)
        seqs2seqnames_dic = dict(zip(shortSequences, sequences))

        jobs = []
        for seq_shortname in shortSequences:
            seq_tof_min = tof_min
            seq_tof_bin = tof_bin
            seq_n_tof_bins = n_tof_bins
            if auto_tof:
                seq_tof_min, seq_tof_bin = self.auto_tof(seqs2seqnames_dic[seq_shortname])
                seq_n_tof_bins = 1

            target = self._fullname_dic[seqs2seqnames_dic[seq_shortname][-1]]
            for t0 in T0s:
                for tof in range(seq_tof_min, seq_tof_min + seq_n_tof_bins * seq_tof_bin, seq_tof_bin):
                    jobs.append({
                        "planet_pack": self.planet_pack,
                        "sequence": seqs2seqnames_dic[seq_shortname],
                        "sequence_short_name": seq_shortname,
                        "t0_min": float(t0),
                        "t0_max": float(t0 + t0_bin),
                        "tof_min": float(tof),
                        "tof_max": float(tof + seq_tof_bin),
                        "mga_t0": [t0 / self._Edy2Kdy, (t0 + t0_bin) / self._Edy2Kdy],
                        "mga_tof": [tof / self._Edy2Kdy, (tof + seq_tof_bin) / self._Edy2Kdy],
                        "vinf": [0.8, 1.8],
                        "tof_encoding": "alpha",
                        "opt_level": opt_level,
                        "opt_injection": opt_injection,
                        "ejection_altitude": float(ejection_altitude),
                        "injection_altitude": float(injection_altitude),
                        "rp_target": self.rp_target_ward(target, injection_altitude),
                        "e_target": self._opt_insertion_dic[opt_injection][2],
                        "add_vinf_dep": True,
                        "add_vinf_arr": self._opt_insertion_dic[opt_injection][0],
                        "orbit_insertion": self._opt_insertion_dic[opt_injection][1],
                        "multi_objective": False,
                        "lambert_max_revs": int(lambert_max_revs),
                        "sade_gen": self._opt_levels_dic[opt_level][0],
                        "n_island": self._opt_levels_dic[opt_level][1],
                        "island_pop": self._opt_levels_dic[opt_level][2],
                        "n_evo_steps": self._opt_levels_dic[opt_level][3],
                        "optimizer_topology": optimizer_topology,
                        "optimizer_seed": optimizer_seed,
                    })

        store = SQLiteJobStore(db_path)
        try:
            batch_id = store.upsert_batch(
                batch_name,
                self.planet_pack,
                template=swing_by_bodies,
                generation_options={
                    "t0_min": t0_min,
                    "t0_bin": t0_bin,
                    "n_t0_bins": n_t0_bins,
                    "tof_min": tof_min,
                    "tof_bin": tof_bin,
                    "n_tof_bins": n_tof_bins,
                    "auto_tof": auto_tof,
                    "opt_level": opt_level,
                    "opt_injection": opt_injection,
                    "injection_altitude": injection_altitude,
                    "ejection_altitude": ejection_altitude,
                    "lambert_max_revs": lambert_max_revs,
                    "optimizer_topology": optimizer_topology,
                    "optimizer_seed": optimizer_seed,
                },
                purpose=purpose,
            )
            for job in jobs:
                store.upsert_job(job, batch_id=batch_id, status="TODO")
        finally:
            store.close()
        return batch_id

    def add_direct_t0_batch_sqlite(
        self, swing_by_bodies, db_path, batch_name="default",
        t0_min=0, t0_bin=100, n_t0_bins=10, tof_profile="relaxed",
        opt_level="debug", arrival_mode="vinf", injection_altitude=100000,
        ejection_altitude=100000, purpose="production",
        optimizer_topology="ring", optimizer_seed=None,
    ):
        """Create direct MGA-1DSM jobs binned only along the launch epoch."""
        from _SQLiteStore import SQLiteJobStore

        if arrival_mode not in self._opt_insertion_dic:
            raise ValueError("Unknown arrival mode: " + str(arrival_mode))
        if batch_name == "default":
            batch_name = self._datastore_name
        sequences = self.generateSequences(swing_by_bodies)
        short_sequences = self.generateShortSequences(swing_by_bodies)
        add_vinf_arr, orbit_insertion, e_target = self._opt_insertion_dic[arrival_mode]
        jobs = []
        for short_name, sequence in zip(short_sequences, sequences):
            estimate = self.estimate_direct_tof_bounds(sequence, profile=tof_profile)
            target = self._fullname_dic[sequence[-1]]
            for bin_index in range(int(n_t0_bins)):
                lower_t0 = float(t0_min + bin_index * t0_bin)
                jobs.append({
                    "planet_pack": self.planet_pack,
                    "sequence": sequence,
                    "sequence_short_name": short_name,
                    "t0_min": lower_t0,
                    "t0_max": lower_t0 + float(t0_bin),
                    "tof_min": float(estimate["total_lower_days"]),
                    "tof_max": float(estimate["total_upper_days"]),
                    "vinf": [0.8, 1.8],
                    "tof_encoding": "direct",
                    "tof_profile": estimate["profile"],
                    "leg_tof_bounds": estimate["direct_bounds_days"],
                    "opt_level": opt_level,
                    "opt_injection": arrival_mode,
                    "arrival_mode": arrival_mode,
                    "ejection_altitude": float(ejection_altitude),
                    "injection_altitude": float(injection_altitude),
                    "rp_target": self.rp_target_ward(target, injection_altitude),
                    "e_target": float(e_target),
                    "add_vinf_dep": True,
                    "add_vinf_arr": bool(add_vinf_arr),
                    "orbit_insertion": bool(orbit_insertion),
                    "multi_objective": False,
                    "lambert_max_revs": 0,
                    "sade_gen": self._opt_levels_dic[opt_level][0],
                    "n_island": self._opt_levels_dic[opt_level][1],
                    "island_pop": self._opt_levels_dic[opt_level][2],
                    "n_evo_steps": self._opt_levels_dic[opt_level][3],
                    "optimizer_topology": optimizer_topology,
                    "optimizer_seed": optimizer_seed,
                })

        store = SQLiteJobStore(db_path)
        try:
            batch_id = store.upsert_batch(
                batch_name, self.planet_pack, template=swing_by_bodies,
                generation_options={
                    "search_axis": "t0", "t0_min": t0_min, "t0_bin": t0_bin,
                    "n_t0_bins": n_t0_bins, "tof_encoding": "direct",
                    "tof_profile": tof_profile, "arrival_mode": arrival_mode,
                    "opt_level": opt_level, "optimizer_topology": optimizer_topology,
                    "optimizer_seed": optimizer_seed,
                },
                purpose=purpose,
            )
            for job in jobs:
                store.upsert_job(job, batch_id=batch_id, status="TODO")
        finally:
            store.close()
        return batch_id

    def add_topology_benchmark_sqlite(
        self,
        swing_by_bodies,
        db_path,
        benchmark_name,
        topologies=("unconnected", "ring", "fully_connected"),
        seeds=range(5),
        **batch_kwargs,
    ):
        """Create benchmark batches for topology comparisons.

        Each topology/seed pair gets its own benchmark batch and distinct jobs,
        so results can be compared without contaminating normal route queries.
        """
        created_batches = []
        for topology in topologies:
            self._make_archipelago_topology(topology, 4)
            for seed in seeds:
                batch_name = f"{benchmark_name}_{topology}_seed{int(seed)}"
                batch_id = self.add_batch_sqlite(
                    swing_by_bodies=swing_by_bodies,
                    db_path=db_path,
                    batch_name=batch_name,
                    purpose="benchmark",
                    optimizer_topology=topology,
                    optimizer_seed=int(seed),
                    **batch_kwargs,
                )
                created_batches.append({
                    "batch_id": batch_id,
                    "batch_name": batch_name,
                    "topology": topology,
                    "seed": int(seed),
                })
        return created_batches

    def load_sqlite_jobs(self, db_path, limit=10, batch_name=None, status="TODO"):
        raise NotImplementedError(
            "load_sqlite_jobs was a temporary dataframe staging path. "
            "Use optimize_sqlite(), which reads jobs directly from SQLite."
        )

    def _optimize_sqlite_job(
        self,
        store,
        job,
        versions,
        auto_workers=True,
        reserve_cores=2,
        topology=None,
        adaptive_stop=None,
        optimizer_strategy=DEFAULT_OPTIMIZER_STRATEGY,
        claim_lease_seconds=86400,
    ):
        run_started_at = time.perf_counter()
        worker_id = job["worker_id"]
        claimed_at = job["claimed_at"]

        def renew_claim():
            if not store.renew_job_claim(
                job["job_id"], worker_id,
                lease_seconds=claim_lease_seconds,
                claimed_at=claimed_at,
            ):
                raise RuntimeError(
                    "Job claim was lost while optimization was running"
                )

        renew_claim()
        optimizer_topology = topology if topology is not None else job.get("optimizer_topology", "ring")
        optimizer_seed = job.get("optimizer_seed")
        sade_gen = int(job["sade_gen"])
        requested_n_island = int(job["n_island"])
        n_island = self._resolve_island_count(
            requested_n_island,
            auto_workers=auto_workers,
            reserve_cores=reserve_cores,
        )
        requested_island_pop = int(job["island_pop"])
        island_pop = max(7, requested_island_pop)
        if island_pop != requested_island_pop:
            logger.warning(
                "Increasing island population from %s to SADE minimum %s",
                requested_island_pop, island_pop,
            )
        n_evo_steps = int(job["n_evo_steps"])
        configured_seed = job.get("optimizer_seed")
        effective_optimizer_seed = (
            int(configured_seed)
            if configured_seed is not None
            else int(np.random.SeedSequence().generate_state(1)[0])
        )
        funnel_config = OptimizationService.funnel_run_config(
            optimizer_strategy,
            n_island,
            island_pop,
            n_evo_steps,
            sade_gen,
        )
        sqlite_run_id = store.start_run(
            job["job_id"],
            versions=versions,
            optimizer_metadata={
                "optimizer_topology": optimizer_topology,
                "optimizer_seed": optimizer_seed,
                "effective_optimizer_seed": effective_optimizer_seed,
                "requested_n_island": requested_n_island,
                "actual_n_island": n_island,
                "island_pop": island_pop,
                "sade_gen": sade_gen,
                "n_evo_steps": n_evo_steps,
                "adaptive_stop": adaptive_stop,
                "optimizer_strategy": optimizer_strategy,
                "funnel_config": funnel_config,
                "code_revision": self._code_revision(),
                "planet_pack_hash": self._planet_pack_hash(),
            },
        )
        try:
            pg.set_global_rng_seed(int(effective_optimizer_seed))
            udp = self._mga_problem_from_sqlite_context(job)
            result = self._run_sqlite_archipelago_job(
                store,
                job,
                udp,
                sqlite_run_id,
                sade_gen=sade_gen,
                requested_n_island=requested_n_island,
                n_island=n_island,
                island_pop=island_pop,
                n_evo_steps=n_evo_steps,
                topology=optimizer_topology,
                versions=versions,
                adaptive_stop=adaptive_stop,
                optimizer_strategy=optimizer_strategy,
                effective_optimizer_seed=effective_optimizer_seed,
                renew_claim=renew_claim,
            )
            runtime_seconds = time.perf_counter() - run_started_at
            store.complete_claimed_run(
                job["job_id"], sqlite_run_id, worker_id, claimed_at,
                result["result_DV"], result["result_t0"],
                result["result_tof"], result["result_ej_vinf"],
                result["gene"], arrival_vinf=result["result_arrival_vinf"],
                runtime_seconds=runtime_seconds,
                actual_n_evo_steps=result["actual_n_evo_steps"],
                stop_reason=result["stop_reason"],
            )
            result["runtime_seconds"] = runtime_seconds
            return result
        except Exception as exc:
            runtime_seconds = time.perf_counter() - run_started_at
            store.fail_claimed_run(
                job["job_id"], sqlite_run_id, worker_id, claimed_at, exc,
                runtime_seconds=runtime_seconds,
            )
            raise

    def _run_sqlite_archipelago_job(
        self,
        store,
        job,
        udp,
        sqlite_run_id,
        sade_gen,
        requested_n_island,
        n_island,
        island_pop,
        n_evo_steps,
        topology,
        versions,
        adaptive_stop=None,
        optimizer_strategy=DEFAULT_OPTIMIZER_STRATEGY,
        effective_optimizer_seed=None,
        renew_claim=None,
    ):
        soi_radius_by_name = {
            name: self.soi_radius(body) for name, body in self._fullname_dic.items()
        }
        funnel_strategies = OptimizationService.FUNNEL_STRATEGIES
        if optimizer_strategy in funnel_strategies:
            return self._run_sqlite_funnel_job(
                store, job, udp, sqlite_run_id,
                sade_gen=sade_gen,
                requested_n_island=requested_n_island,
                n_island=n_island,
                island_pop=island_pop,
                n_evo_steps=n_evo_steps,
                topology=topology,
                versions=versions,
                adaptive_stop=adaptive_stop,
                soi_radius_by_name=soi_radius_by_name,
                exact_strategy=funnel_strategies[optimizer_strategy],
                effective_optimizer_seed=effective_optimizer_seed,
                renew_claim=renew_claim,
            )
        if optimizer_strategy != "flat":
            raise ValueError("Unknown optimizer strategy: " + str(optimizer_strategy))
        fitness_decorator = WayfinderFitnessDecorator(
            planet_pack=self.planet_pack,
            bodies_by_name=self._fullname_dic,
            soi_radius_by_name=soi_radius_by_name,
            ejection_altitude=job["ejection_altitude"],
            tof_encoding=job["tof_encoding"],
        )
        dudp = pg.problem(pg.decorator_problem(udp, fitness_decorator=fitness_decorator))
        uda = pg.sade(gen=sade_gen, memory=True)
        archipelago_topology = self._make_archipelago_topology(topology, n_island)
        archi = pg.archipelago(
            algo=uda,
            prob=dudp,
            n=n_island,
            pop_size=island_pop,
            t=archipelago_topology,
        )
        if n_island != requested_n_island:
            logger.info(
                "Running SADE algo on %s %s islands (preset requested %s)",
                n_island,
                topology,
                requested_n_island,
            )
        else:
            logger.info("Running SADE algo on %s %s islands", n_island, topology)
        convergence_history = []
        plateau_windows = 0
        actual_n_evo_steps = 0
        stop_reason = "max_steps"
        for evo_step in range(1, n_evo_steps + 1):
            previous_populations = [island.get_population() for island in archi]
            previous_id_origins = self._population_id_origins(previous_populations)
            archi.evolve(1)
            archi.wait_check()
            if renew_claim is not None:
                renew_claim()
            populations = [island.get_population() for island in archi]
            telemetry = self._archipelago_telemetry(
                archi, populations, previous_id_origins, topology
            )
            population_fitness = [
                float(value[0])
                for population in populations
                for value in population.get_f()
            ]
            convergence_history.append({
                "best": min(population_fitness),
                "average": float(np.mean(population_fitness)),
            })
            store.record_optimizer_snapshot(
                sqlite_run_id,
                evo_step,
                archi.get_champions_f(),
                archi.get_champions_x(),
                populations=populations,
                telemetry=telemetry,
            )
            actual_n_evo_steps = evo_step
            if adaptive_stop:
                should_stop, plateau_windows = self._adaptive_stop_decision(
                    convergence_history,
                    plateau_windows,
                    adaptive_stop,
                )
                if should_stop:
                    stop_reason = "convergence_plateau"
                    break
        store.record_optimizer_population(
            sqlite_run_id,
            actual_n_evo_steps,
            [island.get_population() for island in archi],
            source="final",
        )
        sols = archi.get_champions_f()
        idx = sols.index(min(sols))
        gene = copy.copy(list(archi.get_champions_x()[idx]))
        decoded = decode_trajectory(udp, gene, planet_pack=self.planet_pack)
        result_dv = decoded["objective_dv"]
        result_t0 = decoded["t0"]
        result_tof = decoded["tof"]
        result_ej_vinf = decoded["ejection_vinf"]
        return {
            "job_id": int(job["job_id"]),
            "run_id": sqlite_run_id,
            "gene": gene,
            "result_DV": result_dv,
            "result_t0": result_t0,
            "result_tof": result_tof,
            "result_ej_vinf": result_ej_vinf,
            "result_arrival_vinf": decoded["arrival_vinf"],
            "topology": topology,
            "n_island": n_island,
            "requested_n_island": requested_n_island,
            "actual_n_evo_steps": actual_n_evo_steps,
            "stop_reason": stop_reason,
            "optimizer_strategy": "flat",
            "stage_summaries": [],
        }

    def _run_sqlite_funnel_job(
        self, store, job, udp, sqlite_run_id, sade_gen,
        requested_n_island, n_island, island_pop, n_evo_steps,
        topology, versions, adaptive_stop, soi_radius_by_name,
        exact_strategy="legacy", effective_optimizer_seed=None,
        renew_claim=None,
    ):
        """Run the corrected v1.5 wide/intermediate/exact-ejection funnel."""
        stage_plan = self._funnel_stage_plan(
            n_island, island_pop, n_evo_steps, sade_gen,
            exact_strategy=exact_strategy,
        )
        seed_genes = None
        exact_archive_genes = []
        global_step = 0
        total_steps = 0
        stage_summaries = []
        final_archipelago = None
        base_seed = (
            int(effective_optimizer_seed)
            if effective_optimizer_seed is not None
            else int(np.random.SeedSequence().generate_state(1)[0])
        )
        local_rng = np.random.default_rng(base_seed)
        exact_fitness_decorator = WayfinderFitnessDecorator(
            planet_pack=self.planet_pack,
            bodies_by_name=self._fullname_dic,
            soi_radius_by_name=soi_radius_by_name,
            ejection_altitude=job["ejection_altitude"],
            tof_encoding=job["tof_encoding"],
            ejection_model="vector_3d",
        )
        exact_problem = pg.problem(
            pg.decorator_problem(
                udp, fitness_decorator=exact_fitness_decorator
            )
        )

        def pygmo_seed(stage_number, island_number, offset=0):
            return int(
                (base_seed * 1000003 + stage_number * 10007
                 + island_number * 101 + offset) % (2 ** 32)
            )

        for stage_index, stage in enumerate(stage_plan, start=1):
            stage_started_at = time.perf_counter()
            stage_topology = stage.get("topology", topology)
            migration_rate = int(stage.get("migration_rate", 2))
            fitness_decorator = WayfinderFitnessDecorator(
                planet_pack=self.planet_pack,
                bodies_by_name=self._fullname_dic,
                soi_radius_by_name=soi_radius_by_name,
                ejection_altitude=job["ejection_altitude"],
                tof_encoding=job["tof_encoding"],
                ejection_model=stage["ejection_model"],
            )
            problem = pg.problem(
                pg.decorator_problem(udp, fitness_decorator=fitness_decorator)
            )
            archi = pg.archipelago(
                t=self._make_archipelago_topology(
                    stage_topology, stage["n_island"]
                ),
            )
            phase_elite_groups = None
            if seed_genes:
                if stage["initialization"] == "scout_diverse":
                    seed_genes = self._select_phase_diverse_elites(
                        job,
                        problem,
                        seed_genes,
                        stage["n_island"],
                        elite_fraction=0.75,
                    )
                elif stage["initialization"] == "phase_elites_mixed":
                    selected_elites = self._select_phase_diverse_elites(
                        job,
                        problem,
                        seed_genes,
                        stage["n_island"] * stage["island_pop"],
                        elite_fraction=stage.get("elite_fraction", 0.35),
                    )
                    phase_elite_groups = [
                        selected_elites[index::stage["n_island"]]
                        for index in range(stage["n_island"])
                    ]
                    if any(
                        len(group) != stage["island_pop"]
                        for group in phase_elite_groups
                    ):
                        raise ValueError(
                            "Insufficient phase elites for stage population"
                        )
                elif stage["ejection_model"] == "vector_3d":
                    seed_genes = self._select_exact_diverse_seeds(
                        problem, seed_genes, stage["n_island"]
                    )
                else:
                    seed_genes = sorted(
                        seed_genes,
                        key=lambda gene: float(problem.fitness(gene)[0]),
                    )[:stage["n_island"]]
            stage_algorithms = []
            for island_index in range(stage["n_island"]):
                if phase_elite_groups is not None:
                    population = pg.population(problem, size=0)
                    for elite_gene in phase_elite_groups[island_index]:
                        population.push_back(elite_gene)
                elif seed_genes:
                    use_global_population = (
                        stage["initialization"] == "local_global"
                        and stage["n_island"] > 1
                        and island_index == stage["n_island"] - 1
                    )
                    if stage["initialization"] in ("local", "local_global") and not use_global_population:
                        population = self._seeded_population(
                            problem,
                            seed_genes[island_index],
                            stage["island_pop"],
                            sigma_fraction=stage["sigma_fraction"],
                            rng=local_rng,
                        )
                    else:
                        population = pg.population(
                            problem, size=max(0, stage["island_pop"] - 1),
                            seed=pygmo_seed(stage_index, island_index, 17),
                        )
                        population.push_back(seed_genes[island_index])
                else:
                    population = pg.population(
                        problem, size=stage["island_pop"],
                        seed=pygmo_seed(stage_index, island_index, 17),
                    )

                if stage["algorithm"] == "sade_nlopt" and island_index % 2 == 1:
                    algorithm = pg.nlopt("neldermead")
                    algorithm.maxeval = max(
                        1, stage["sade_gen"] * stage["island_pop"] - 1
                    )
                    algorithm.xtol_rel = 0.0
                    algorithm.ftol_rel = 0.0
                    algorithm.selection = "best"
                    algorithm.replacement = "best"
                    algorithm_name = "nlopt_neldermead"
                elif stage["algorithm"] in ("sade", "sade_nlopt") or island_index % 2 == 0:
                    algorithm = pg.sade(
                        gen=stage["sade_gen"],
                        memory=True,
                        seed=pygmo_seed(stage_index, island_index, 31),
                    )
                    algorithm_name = "sade"
                else:
                    annealing_schedule = self._balanced_annealing_schedule(
                        problem, stage["island_pop"], stage["sade_gen"]
                    )
                    algorithm = pg.simulated_annealing(
                        Ts=stage["annealing"]["Ts"],
                        Tf=stage["annealing"]["Tf"],
                        n_T_adj=annealing_schedule["n_T_adj"],
                        n_range_adj=annealing_schedule["n_range_adj"],
                        bin_size=annealing_schedule["bin_size"],
                        seed=pygmo_seed(stage_index, island_index, 31),
                    )
                    algorithm_name = "simulated_annealing"
                if algorithm_name not in stage_algorithms:
                    stage_algorithms.append(algorithm_name)
                archi.push_back(pg.island(
                    algo=algorithm,
                    pop=population,
                    s_pol=pg.select_best(migration_rate),
                    r_pol=pg.fair_replace(migration_rate),
                ))

            logger.info(
                "Funnel stage %s (%s): %s islands, pop=%s, algorithm=%s, init=%s",
                stage_index, stage["name"], stage["n_island"], stage["island_pop"],
                stage["algorithm"], stage["initialization"],
            )
            history = []
            plateau_windows = 0
            stage_steps = 0
            stage_stop_reason = "max_steps"
            stage_adaptive_stop = adaptive_stop or stage.get("adaptive_stop")
            stage_migrants_published = 0
            stage_migration_islands_active = 0
            stage_migrations_accepted = 0
            stage_telemetry = None
            for _ in range(stage["evo_steps"]):
                previous_populations = [island.get_population() for island in archi]
                previous_id_origins = self._population_id_origins(
                    previous_populations
                )
                archi.evolve(1)
                archi.wait_check()
                if renew_claim is not None:
                    renew_claim()
                global_step += 1
                stage_steps += 1
                populations = [island.get_population() for island in archi]
                stage_telemetry = self._archipelago_telemetry(
                    archi, populations, previous_id_origins, stage_topology
                )
                stage_migrants_published += stage_telemetry["migrants_published"]
                stage_migration_islands_active += stage_telemetry[
                    "migration_islands_active"
                ]
                stage_migrations_accepted += stage_telemetry[
                    "migrations_accepted"
                ]
                population_fitness = [
                    float(value[0])
                    for population in populations
                    for value in population.get_f()
                ]
                history.append({
                    "best": min(population_fitness),
                    "average": float(np.mean(population_fitness)),
                })
                store.record_optimizer_snapshot(
                    sqlite_run_id,
                    global_step,
                    archi.get_champions_f(),
                    archi.get_champions_x(),
                    populations=populations,
                    telemetry=stage_telemetry,
                )
                archive_interval = int(stage.get("archive_interval", 0))
                if (
                    stage.get("archive_exact")
                    and archive_interval > 0
                    and stage_steps % archive_interval == 0
                ):
                    archive_candidates = []
                    for population in populations:
                        ranked = sorted(
                            zip(population.get_f(), population.get_x()),
                            key=lambda pair: float(pair[0][0]),
                        )
                        archive_candidates.extend(
                            list(gene) for _, gene in ranked[:2]
                        )
                    exact_archive_genes = self._update_exact_phase_archive(
                        job,
                        exact_problem,
                        exact_archive_genes,
                        archive_candidates,
                        max_size=stage.get("archive_size", 32),
                    )
                if stage_adaptive_stop:
                    should_stop, plateau_windows = self._adaptive_stop_decision(
                        history, plateau_windows, stage_adaptive_stop
                    )
                    if should_stop:
                        stage_stop_reason = "convergence_plateau"
                        break

            populations = [island.get_population() for island in archi]
            store.record_optimizer_population(
                sqlite_run_id,
                global_step,
                populations,
                source="stage_{}_final".format(stage_index),
            )
            champion_pairs = sorted(
                zip(archi.get_champions_f(), archi.get_champions_x()),
                key=lambda pair: float(pair[0][0]),
            )
            if stage_index < len(stage_plan):
                next_stage = stage_plan[stage_index]
                next_island_count = next_stage["n_island"]
                if (
                    next_stage["ejection_model"] == "vector_3d"
                    or next_stage["initialization"] == "phase_elites_mixed"
                ):
                    seed_genes = [
                        list(gene)
                        for population in populations
                        for gene in population.get_x()
                    ]
                    if next_stage.get("use_exact_archive"):
                        seed_genes = exact_archive_genes + seed_genes
                        if exact_archive_genes:
                            archive_population = pg.population(
                                exact_problem, size=0
                            )
                            for archive_gene in exact_archive_genes:
                                archive_population.push_back(archive_gene)
                            store.record_optimizer_population(
                                sqlite_run_id,
                                global_step,
                                [archive_population],
                                source="exact_archive_stage3_seed",
                            )
                elif next_stage["initialization"] == "scout_diverse":
                    seed_genes = [
                        list(gene) for _, gene in champion_pairs
                    ]
                else:
                    seed_genes = [
                        list(champion_pairs[index][1])
                        for index in range(next_island_count)
                    ]
            else:
                seed_genes = None
            stage_runtime = time.perf_counter() - stage_started_at
            stage_summary = {
                "stage_index": stage_index,
                "stage_name": stage["name"],
                "n_island": stage["n_island"],
                "island_pop": stage["island_pop"],
                "sade_gen": stage["sade_gen"],
                "planned_evo_steps": stage["evo_steps"],
                "actual_evo_steps": stage_steps,
                "ejection_model": stage["ejection_model"],
                "initialization": stage["initialization"],
                "topology": stage_topology,
                "migration_rate": migration_rate,
                "exact_archive_size": len(exact_archive_genes),
                "adaptive_stop": stage_adaptive_stop,
                "algorithms": stage_algorithms,
                "topology_vertices": stage_telemetry["topology_vertices"],
                "topology_edges": stage_telemetry["topology_edges"],
                "migrants_published": stage_migrants_published,
                "migration_islands_active": stage_migration_islands_active,
                "migrations_accepted": stage_migrations_accepted,
                "best_fitness": float(champion_pairs[0][0][0]),
                "average_fitness": history[-1]["average"],
                "runtime_seconds": stage_runtime,
                "stop_reason": stage_stop_reason,
            }
            store.record_optimizer_stage(sqlite_run_id, stage_summary)
            stage_summaries.append(stage_summary)
            total_steps += stage_steps
            final_archipelago = archi

        sols = final_archipelago.get_champions_f()
        idx = sols.index(min(sols))
        gene = copy.copy(list(final_archipelago.get_champions_x()[idx]))
        decoded = decode_trajectory(udp, gene, planet_pack=self.planet_pack)
        return {
            "job_id": int(job["job_id"]),
            "run_id": sqlite_run_id,
            "gene": gene,
            "result_DV": decoded["objective_dv"],
            "result_t0": decoded["t0"],
            "result_tof": decoded["tof"],
            "result_ej_vinf": decoded["ejection_vinf"],
            "result_arrival_vinf": decoded["arrival_vinf"],
            "topology": topology,
            "n_island": n_island,
            "requested_n_island": requested_n_island,
            "actual_n_evo_steps": total_steps,
            "stop_reason": "funnel_complete",
            "optimizer_strategy": (
                "funnel" if exact_strategy == "legacy"
                else (
                    (
                        "funnel_phase_elites_equal"
                        if exact_strategy == "phase_elites_nm_equal"
                        else (
                            {
                                "scout_archive_nm": "funnel_scout_archive",
                                "scout_archive_nm_32": "funnel_scout_archive_32",
                                "scout_archive_nm_64": "funnel_scout_archive_64",
                                "scout_archive_nm_128": "funnel_scout_archive_128",
                            }[exact_strategy]
                            if exact_strategy.startswith("scout_archive_nm")
                            else "funnel_phase_elites_nm"
                        )
                    )
                    if exact_strategy in (
                        "phase_elites_nm", "phase_elites_nm_equal",
                        "scout_archive_nm", "scout_archive_nm_32",
                        "scout_archive_nm_64", "scout_archive_nm_128",
                    )
                    else "funnel_{}_exact".format(exact_strategy)
                )
            ),
            "stage_summaries": stage_summaries,
        }

    @staticmethod
    def _adaptive_stop_decision(history, plateau_windows, options):
        """Return whether best and population-average progress have plateaued."""
        min_steps = int(options.get("min_steps", 5))
        window = int(options.get("window", 3))
        patience = int(options.get("patience", 2))
        best_tol = float(options.get("best_relative_tolerance", 0.005))
        average_tol = float(options.get("average_relative_tolerance", 0.02))
        require_average = bool(options.get("require_average_plateau", True))
        if len(history) < max(min_steps, window + 1):
            return False, 0

        old = history[-window - 1]
        new = history[-1]
        best_improvement = (old["best"] - new["best"]) / max(abs(old["best"]), 1.0)
        average_improvement = (
            old["average"] - new["average"]
        ) / max(abs(old["average"]), 1.0)
        plateau = best_improvement < best_tol and (
            not require_average or average_improvement < average_tol
        )
        plateau_windows = plateau_windows + 1 if plateau else 0
        return plateau_windows >= patience, plateau_windows

    def plot_optimizer_convergence_sqlite(
        self, db_path, run_id, figsize=(9, 5), log_y=False
    ):
        """Plot best and whole-population average fitness versus evolution step."""
        from _SQLiteStore import SQLiteJobStore

        store = SQLiteJobStore(db_path)
        try:
            rows = store.optimizer_convergence(run_id)
        finally:
            store.close()
        if not rows:
            raise ValueError("No optimizer snapshots for run_id: " + str(run_id))

        steps = [row["step"] for row in rows]
        best = [row["best_fitness"] for row in rows]
        average = [row["average_fitness"] for row in rows]
        fig, ax = plt.subplots(figsize=figsize)
        ax.plot(steps, best, label="Best DV", linewidth=2)
        if any(value is not None for value in average):
            ax.plot(steps, average, label="Average population DV", linewidth=2)
        ax.set_xlabel("Evolution step")
        ax.set_ylabel("Objective DV (m/s)")
        ax.set_title("Optimizer convergence — run " + str(run_id))
        ax.grid(True, alpha=0.25)
        ax.legend()
        if log_y:
            ax.set_yscale("log")
        fig.tight_layout()
        return fig, ax

    def optimize_sqlite(
        self,
        db_path,
        n=10,
        batch_name=None,
        auto_workers=True,
        reserve_cores=2,
        topology=None,
        adaptive_stop=False,
        adaptive_min_steps=5,
        adaptive_window=3,
        adaptive_patience=2,
        adaptive_best_relative_tolerance=0.005,
        adaptive_average_relative_tolerance=0.02,
        optimizer_strategy=DEFAULT_OPTIMIZER_STRATEGY,
        worker_id=None,
        claim_lease_seconds=86400,
    ):
        """Atomically claim and optimize SQLite jobs, then persist results."""
        from _SQLiteStore import SQLiteJobStore

        versions = {
            "pykep": getattr(pykep, "__version__", None),
            "pygmo": getattr(pg, "__version__", None),
            "wayfinder": WAYFINDER_VERSION,
        }
        store = SQLiteJobStore(db_path)
        try:
            count = 0
            max_jobs = max(0, int(n))
            adaptive_options = None
            if adaptive_stop:
                adaptive_options = {
                    "min_steps": adaptive_min_steps,
                    "window": adaptive_window,
                    "patience": adaptive_patience,
                    "best_relative_tolerance": adaptive_best_relative_tolerance,
                    "average_relative_tolerance": adaptive_average_relative_tolerance,
                }
            while count < max_jobs:
                jobs = store.claim_pending_jobs(
                    self.planet_pack,
                    limit=1,
                    batch_name=batch_name,
                    worker_id=worker_id,
                    lease_seconds=claim_lease_seconds,
                )
                if not jobs:
                    if count == 0:
                        logger.warning(
                            "No claimable TODO jobs found in SQLite datastore"
                        )
                    break
                job = jobs[0]
                self._optimize_sqlite_job(
                    store,
                    job,
                    versions,
                    auto_workers=auto_workers,
                    reserve_cores=reserve_cores,
                    topology=topology,
                    adaptive_stop=adaptive_options,
                    optimizer_strategy=optimizer_strategy,
                    claim_lease_seconds=claim_lease_seconds,
                )
                count += 1
                logger.info(
                    "Iteration %s complete: seq=%s lb_t0=%s lb_tof=%s",
                    count,
                    job["sequence_short_name"],
                    job["t0_min"],
                    job["tof_min"],
                )
            return count
        finally:
            store.close()

    def benchmark_results_sqlite(self, db_path, benchmark_name=None):
        """Return benchmark runs as a dataframe for comparison/plotting."""
        from _SQLiteStore import SQLiteJobStore

        store = SQLiteJobStore(db_path)
        try:
            return pd.DataFrame(
                store.benchmark_results(
                    planet_pack=self.planet_pack,
                    benchmark_name=benchmark_name,
                )
            )
        finally:
            store.close()


    def edit_batch(self,swing_by_bodies,action = 'reset',save_it = True):
        raise NotImplementedError(
            "edit_batch was tied to the legacy dataframe datastore. "
            "Use SQLite job queries/updates instead."
        )


    def build_one_mgaproblem(self,seq_name,t0,tof,):
        raise NotImplementedError(
            "build_one_mgaproblem was tied to the legacy dataframe staging path. "
            "Build MGA problems from SQLite job dictionaries instead."
        )
        
    
    def build_all_mgaproblems(self):
        raise NotImplementedError(
            "build_all_mgaproblems was tied to the legacy dataframe staging path. "
            "Use optimize_sqlite(), which builds problems directly from SQLite jobs."
        )


    def _mga_problem_from_sqlite_context(self, context):
        bodies = json.loads(context["bodies_json"])
        if context["tof_encoding"] == "direct":
            serialized_bounds = context.get("leg_tof_bounds_json")
            if not serialized_bounds:
                raise ValueError("Direct job is missing per-leg ToF bounds")
            leg_bounds = json.loads(serialized_bounds)
            if len(leg_bounds) != len(bodies) - 1:
                raise ValueError("Direct job ToF bounds do not match its number of legs")
            tof = [[float(lo) / self._Edy2Kdy, float(hi) / self._Edy2Kdy]
                   for lo, hi in leg_bounds]
        else:
            tof = [context["tof_min"] / self._Edy2Kdy,
                   context["tof_max"] / self._Edy2Kdy]
        return _make_mga_1dsm(
            seq=list(map(self._fullname_dic.get, bodies)),
            t0=[context["t0_min"] / self._Edy2Kdy, context["t0_max"] / self._Edy2Kdy],
            tof=tof,
            vinf=[context["vinf_min"], context["vinf_max"]],
            tof_encoding=context["tof_encoding"],
            add_vinf_dep=bool(context["add_vinf_dep"]),
            add_vinf_arr=bool(context["add_vinf_arr"]),
            multi_objective=bool(context["multi_objective"]),
            orbit_insertion=bool(context["orbit_insertion"]),
            rp_target=context["rp_target"],
            e_target=context["e_target"],
            max_revs=int(context["lambert_max_revs"]),
        )

    def _direct_local_problem_from_sqlite_context(
        self,
        context,
        alpha_gene,
        t0_wiggle_days=5,
        leg_tof_wiggle_days=5,
        t0_bounds_days=None,
    ):
        bodies = json.loads(context["bodies_json"])
        n_legs = len(bodies) - 1
        decoded_leg_tofs = alpha_leg_tofs(alpha_gene, n_legs)
        direct_gene = alpha_gene_to_direct_gene(alpha_gene, n_legs)
        direct_tof_bounds = direct_tof_bounds_from_leg_tofs(
            decoded_leg_tofs,
            float(leg_tof_wiggle_days) / self._Edy2Kdy,
        )
        if t0_bounds_days is None:
            center_t0 = float(alpha_gene[0])
            t0_bounds = [
                max(float(context["t0_min"]) / self._Edy2Kdy, center_t0 - float(t0_wiggle_days) / self._Edy2Kdy),
                min(float(context["t0_max"]) / self._Edy2Kdy, center_t0 + float(t0_wiggle_days) / self._Edy2Kdy),
            ]
        else:
            t0_bounds = [
                max(float(context["t0_min"]), float(t0_bounds_days[0])) / self._Edy2Kdy,
                min(float(context["t0_max"]), float(t0_bounds_days[1])) / self._Edy2Kdy,
            ]
            if t0_bounds[1] <= t0_bounds[0]:
                raise ValueError("Direct local T0 bounds are empty")
        udp = _make_mga_1dsm(
            seq=list(map(self._fullname_dic.get, bodies)),
            t0=t0_bounds,
            tof=direct_tof_bounds,
            vinf=[context["vinf_min"], context["vinf_max"]],
            tof_encoding="direct",
            add_vinf_dep=bool(context["add_vinf_dep"]),
            add_vinf_arr=bool(context["add_vinf_arr"]),
            multi_objective=bool(context["multi_objective"]),
            orbit_insertion=bool(context["orbit_insertion"]),
            rp_target=context["rp_target"],
            e_target=context["e_target"],
            max_revs=int(context["lambert_max_revs"]),
        )
        return udp, direct_gene, decoded_leg_tofs

    def reoptimize_direct_near_run_sqlite(
        self,
        db_path,
        run_id,
        t0_wiggle_days=5,
        leg_tof_wiggle_days=5,
        sade_gen=30,
        pop_size=48,
        seed_alpha_gene=None,
    ):
        """Run a small direct-mode optimization near an existing alpha-mode SQL run."""
        from _SQLiteStore import SQLiteJobStore

        store = SQLiteJobStore(db_path)
        try:
            context = store.run_context(run_id)
        finally:
            store.close()
        if context is None:
            raise ValueError("Unknown run_id: " + str(run_id))
        if context["tof_encoding"] != "alpha":
            raise ValueError("Direct local reopt currently expects an alpha source run")
        if seed_alpha_gene is None:
            if not context.get("gene_json"):
                raise ValueError("Run has no stored alpha gene: " + str(run_id))
            seed_alpha_gene = json.loads(context["gene_json"])

        udp, direct_gene, decoded_leg_tofs = self._direct_local_problem_from_sqlite_context(
            context,
            seed_alpha_gene,
            t0_wiggle_days=t0_wiggle_days,
            leg_tof_wiggle_days=leg_tof_wiggle_days,
        )
        soi_radius_by_name = {
            name: self.soi_radius(body) for name, body in self._fullname_dic.items()
        }
        fitness_decorator = WayfinderFitnessDecorator(
            planet_pack=self.planet_pack,
            bodies_by_name=self._fullname_dic,
            soi_radius_by_name=soi_radius_by_name,
            ejection_altitude=context["ejection_altitude"],
            tof_encoding="direct",
        )
        decorated_problem = pg.problem(pg.decorator_problem(udp, fitness_decorator=fitness_decorator))
        population = pg.population(decorated_problem, size=max(1, int(pop_size)))
        population.set_x(0, direct_gene)
        algorithm = pg.algorithm(pg.sade(gen=int(sade_gen)))
        population = algorithm.evolve(population)
        champion_gene = list(population.champion_x)
        champion_fitness = float(population.champion_f[0])
        decoded = decode_dV_tof(udp, champion_gene, planet_pack=self.planet_pack)
        seed_decoded = decode_dV_tof(udp, direct_gene, planet_pack=self.planet_pack)
        return {
            "source_run_id": int(run_id),
            "source_alpha_gene": seed_alpha_gene,
            "direct_seed_gene": direct_gene,
            "direct_champion_gene": champion_gene,
            "decoded_leg_tofs_days": [tof * self._Edy2Kdy for tof in decoded_leg_tofs],
            "seed_objective": seed_decoded[0],
            "seed_t0": seed_decoded[1],
            "seed_tof": seed_decoded[2],
            "seed_ejection_vinf": seed_decoded[3],
            "champion_fitness": champion_fitness,
            "champion_objective": decoded[0],
            "champion_t0": decoded[1],
            "champion_tof": decoded[2],
            "champion_ejection_vinf": decoded[3],
        }

    def _seeded_population(self, problem, seed_gene, size, sigma_fraction=0.12, rng=None):
        population = pg.population(problem, size=max(1, int(size)))
        lower_bounds, upper_bounds = problem.get_bounds()
        lower_bounds = np.asarray(lower_bounds, dtype=float)
        upper_bounds = np.asarray(upper_bounds, dtype=float)
        seed_gene = np.asarray(seed_gene, dtype=float)
        rng = np.random.default_rng() if rng is None else rng
        span = upper_bounds - lower_bounds
        population.set_x(0, seed_gene.tolist())
        for index in range(1, int(size)):
            perturbation = rng.normal(0.0, float(sigma_fraction), size=len(seed_gene)) * span
            candidate = np.clip(seed_gene + perturbation, lower_bounds, upper_bounds)
            population.set_x(index, candidate.tolist())
        return population

    def sample_direct_local_cloud_sqlite(
        self,
        db_path,
        run_id,
        sampler_name="direct_local_cloud",
        t0_wiggle_days=5,
        leg_tof_wiggle_days=5,
        sade_gen=40,
        n_island=8,
        island_pop=32,
        seed_alpha_gene=None,
        initialization="seeded",
        sigma_fraction=0.12,
        elite_fraction=0.25,
        max_elites_per_population=None,
        random_seed=None,
    ):
        """Sample a locally reoptimized direct-mode population cloud around an alpha solution."""
        from _SQLiteStore import SQLiteJobStore

        store = SQLiteJobStore(db_path)
        try:
            context = store.run_context(run_id)
        finally:
            store.close()
        if context is None:
            raise ValueError("Unknown run_id: " + str(run_id))
        if context["tof_encoding"] != "alpha":
            raise ValueError("Direct local cloud currently expects an alpha source run")
        if seed_alpha_gene is None:
            if not context.get("gene_json"):
                raise ValueError("Run has no stored alpha gene: " + str(run_id))
            seed_alpha_gene = json.loads(context["gene_json"])

        udp, direct_gene, decoded_leg_tofs = self._direct_local_problem_from_sqlite_context(
            context,
            seed_alpha_gene,
            t0_wiggle_days=t0_wiggle_days,
            leg_tof_wiggle_days=leg_tof_wiggle_days,
        )
        soi_radius_by_name = {
            name: self.soi_radius(body) for name, body in self._fullname_dic.items()
        }
        fitness_decorator = WayfinderFitnessDecorator(
            planet_pack=self.planet_pack,
            bodies_by_name=self._fullname_dic,
            soi_radius_by_name=soi_radius_by_name,
            ejection_altitude=context["ejection_altitude"],
            tof_encoding="direct",
        )
        decorated_problem = pg.problem(pg.decorator_problem(udp, fitness_decorator=fitness_decorator))
        algorithm = pg.algorithm(pg.sade(gen=int(sade_gen)))
        populations = []
        rng = np.random.default_rng(random_seed)
        for island_index in range(int(n_island)):
            if initialization == "seeded":
                population = self._seeded_population(
                    decorated_problem,
                    direct_gene,
                    int(island_pop),
                    sigma_fraction=sigma_fraction,
                    rng=rng,
                )
            elif initialization == "random":
                population = pg.population(decorated_problem, size=int(island_pop))
                population.set_x(0, direct_gene)
            else:
                raise ValueError("Unknown direct cloud initialization: " + str(initialization))
            populations.append(algorithm.evolve(population))

        rows = []
        for population in populations:
            genes = population.get_x()
            fitnesses = population.get_f()
            fitness_order = sorted(
                range(len(fitnesses)),
                key=lambda idx: fitnesses[idx][0] if isinstance(fitnesses[idx], (list, tuple, np.ndarray)) else fitnesses[idx],
            )
            elite_count = max(1, int(np.ceil(len(fitness_order) * float(elite_fraction))))
            if max_elites_per_population is not None:
                elite_count = min(elite_count, int(max_elites_per_population))
            elite_indexes = fitness_order[:elite_count]
            for individual_index in elite_indexes:
                gene = genes[individual_index]
                fitness_vector = fitnesses[individual_index]
                gene = list(gene)
                precise_t0 = float(gene[0]) * self._Edy2Kdy
                precise_tof = sum(direct_leg_tofs(gene, len(decoded_leg_tofs))) * self._Edy2Kdy
                try:
                    result_dv, result_t0, result_tof, ejection_vinf = decode_dV_tof(
                        udp,
                        gene,
                        planet_pack=self.planet_pack,
                    )
                    decode_error = None
                except (ValueError, ZeroDivisionError, FloatingPointError) as exc:
                    result_dv = None
                    ejection_vinf = None
                    decode_error = str(exc)
                fitness = fitness_vector[0] if isinstance(fitness_vector, (list, tuple, np.ndarray)) else fitness_vector
                rows.append({
                    "t0": precise_t0,
                    "tof": precise_tof,
                    "metric": result_dv,
                    "fitness": float(fitness),
                    "ejection_vinf": ejection_vinf,
                    "gene": gene,
                    "decode_error": decode_error,
                })

        deduped_rows = {}
        for row in rows:
            key = (round(float(row["t0"]), 9), round(float(row["tof"]), 9))
            previous = deduped_rows.get(key)
            if previous is None:
                deduped_rows[key] = row
                continue
            previous_metric = float("inf") if previous.get("metric") is None else float(previous["metric"])
            row_metric = float("inf") if row.get("metric") is None else float(row["metric"])
            if row_metric < previous_metric:
                deduped_rows[key] = row
        rows = list(deduped_rows.values())

        store = SQLiteJobStore(db_path)
        try:
            store.replace_porkchop_samples(
                run_id,
                sampler_name,
                rows,
                metadata={
                    "sampler_type": "direct_local_cloud",
                    "source_run_id": int(run_id),
                    "t0_wiggle_days": t0_wiggle_days,
                    "leg_tof_wiggle_days": leg_tof_wiggle_days,
                    "sade_gen": sade_gen,
                    "n_island": n_island,
                    "island_pop": island_pop,
                    "initialization": initialization,
                    "sigma_fraction": sigma_fraction,
                    "elite_fraction": elite_fraction,
                    "max_elites_per_population": max_elites_per_population,
                    "random_seed": random_seed,
                    "sample_count": len(rows),
                },
            )
            samples = store.porkchop_samples(run_id, sampler_name)
        finally:
            store.close()
        return pd.DataFrame(samples)

    def diverse_alpha_seeds_from_optimizer_sqlite(
        self,
        db_path,
        run_id,
        max_seeds=10,
        dv_margin=1000,
        t0_bin_days=1.0,
        tof_bin_days=1.0,
        include_final_population=True,
    ):
        """Select diverse good alpha genes from stored optimizer telemetry."""
        points = self.porkchop_points_from_snapshots_sqlite(
            db_path,
            run_id,
            include_final_population=include_final_population,
        )
        if points.empty:
            return points
        valid_points = points.dropna(subset=["result_DV", "result_t0", "result_tof"]).copy()
        if valid_points.empty:
            return valid_points

        best_dv = valid_points["result_DV"].min()
        valid_points = valid_points[valid_points["result_DV"] <= best_dv + float(dv_margin)].copy()
        valid_points["t0_bin"] = np.floor(valid_points["result_t0"] / float(t0_bin_days)) * float(t0_bin_days)
        valid_points["tof_bin"] = np.floor(valid_points["result_tof"] / float(tof_bin_days)) * float(tof_bin_days)
        idx = valid_points.groupby(["t0_bin", "tof_bin"])["result_DV"].idxmin()
        seeds = valid_points.loc[idx].sort_values("result_DV").head(int(max_seeds)).copy()
        return seeds.reset_index(drop=True)

    def expand_diverse_alpha_seeds_direct_sqlite(
        self,
        db_path,
        run_id,
        sampler_prefix="direct_seed_expand",
        max_seeds=10,
        dv_margin=1000,
        seed_t0_bin_days=1.0,
        seed_tof_bin_days=1.0,
        t0_wiggle_days=3,
        leg_tof_wiggle_days=3,
        sade_gen=80,
        n_island=12,
        island_pop=32,
        sigma_fraction=0.08,
        elite_fraction=0.25,
        max_elites_per_population=8,
        random_seed=None,
    ):
        """Expand diverse alpha optimizer seeds into local direct elite clouds."""
        seeds = self.diverse_alpha_seeds_from_optimizer_sqlite(
            db_path,
            run_id,
            max_seeds=max_seeds,
            dv_margin=dv_margin,
            t0_bin_days=seed_t0_bin_days,
            tof_bin_days=seed_tof_bin_days,
        )
        sampler_names = []
        frames = []
        for seed_index, seed in seeds.iterrows():
            sampler_name = f"{sampler_prefix}_{int(seed_index):03d}"
            seed_random = None if random_seed is None else int(random_seed) + int(seed_index)
            frame = self.sample_direct_local_cloud_sqlite(
                db_path,
                run_id,
                sampler_name=sampler_name,
                t0_wiggle_days=t0_wiggle_days,
                leg_tof_wiggle_days=leg_tof_wiggle_days,
                sade_gen=sade_gen,
                n_island=n_island,
                island_pop=island_pop,
                seed_alpha_gene=seed["gene"],
                initialization="seeded",
                sigma_fraction=sigma_fraction,
                elite_fraction=elite_fraction,
                max_elites_per_population=max_elites_per_population,
                random_seed=seed_random,
            )
            sampler_names.append(sampler_name)
            frames.append(frame)
        if frames:
            combined = pd.concat(frames, ignore_index=True)
        else:
            combined = pd.DataFrame()
        return {
            "seeds": seeds,
            "sampler_names": sampler_names,
            "points": combined,
        }

    def adaptive_direct_cell_refinement_sqlite(
        self,
        db_path,
        run_id,
        sampler_name="direct_adaptive_cells",
        center_t0=None,
        center_tof=None,
        span_t0_days=40,
        span_tof_days=40,
        initial_cell_days=4.0,
        min_cell_days=1.0,
        refine_threshold_dv=6000,
        alpha_seed_dv_margin=1500,
        alpha_seed_max=40,
        seed_search_t0_scale=2.0,
        seed_search_tof_scale=2.0,
        t0_wiggle_factor=0.75,
        leg_tof_wiggle_days=2.0,
        tof_window_penalty=2000.0,
        sade_gen=35,
        pop_size=32,
        sigma_fraction=0.08,
        random_seed=None,
        max_cells=200,
    ):
        """Adaptively refine a local T0/TOF grid using direct champion re-optimization."""
        from _SQLiteStore import SQLiteJobStore

        store = SQLiteJobStore(db_path)
        try:
            context = store.run_context(run_id)
        finally:
            store.close()
        if context is None:
            raise ValueError("Unknown run_id: " + str(run_id))
        if center_t0 is None:
            center_t0 = float(context["result_t0"])
        if center_tof is None:
            center_tof = float(context["result_tof"])

        alpha_seeds = self.diverse_alpha_seeds_from_optimizer_sqlite(
            db_path,
            run_id,
            max_seeds=alpha_seed_max,
            dv_margin=alpha_seed_dv_margin,
            t0_bin_days=max(0.5, float(initial_cell_days) / 2),
            tof_bin_days=max(0.5, float(initial_cell_days) / 2),
        )
        if alpha_seeds.empty:
            raise ValueError("No alpha seeds available for adaptive direct refinement")

        seed_rows = alpha_seeds.to_dict("records")
        rng = np.random.default_rng(random_seed)
        soi_radius_by_name = {
            name: self.soi_radius(body) for name, body in self._fullname_dic.items()
        }
        rows = []
        evaluated = 0

        def choose_seed(cell_center_t0, cell_center_tof, cell_days):
            t0_scale = max(float(cell_days) * float(seed_search_t0_scale), 1e-9)
            tof_scale = max(float(cell_days) * float(seed_search_tof_scale), 1e-9)
            return min(
                seed_rows,
                key=lambda seed: (
                    ((float(seed["result_t0"]) - cell_center_t0) / t0_scale) ** 2
                    + ((float(seed["result_tof"]) - cell_center_tof) / tof_scale) ** 2
                    + (float(seed["result_DV"]) - float(alpha_seeds["result_DV"].min())) / 10000.0
                ),
            )

        def evaluate_cell(t0_low, t0_high, tof_low, tof_high, cell_days):
            nonlocal evaluated
            if evaluated >= int(max_cells):
                return None
            evaluated += 1
            cell_center_t0 = (t0_low + t0_high) / 2
            cell_center_tof = (tof_low + tof_high) / 2
            seed = choose_seed(cell_center_t0, cell_center_tof, cell_days)
            seed_gene = list(seed["gene"])
            udp, direct_gene, decoded_leg_tofs = self._direct_local_problem_from_sqlite_context(
                context,
                seed_gene,
                t0_wiggle_days=max(float(cell_days) * float(t0_wiggle_factor), 0.25),
                leg_tof_wiggle_days=leg_tof_wiggle_days,
                t0_bounds_days=[t0_low, t0_high],
            )
            lower_bounds, upper_bounds = pg.problem(udp).get_bounds()
            direct_gene[0] = float(np.clip(float(cell_center_t0) / self._Edy2Kdy, lower_bounds[0], upper_bounds[0]))

            base_fitness_decorator = WayfinderFitnessDecorator(
                planet_pack=self.planet_pack,
                bodies_by_name=self._fullname_dic,
                soi_radius_by_name=soi_radius_by_name,
                ejection_altitude=context["ejection_altitude"],
                tof_encoding="direct",
            )
            n_legs = len(json.loads(context["bodies_json"])) - 1

            def fitness_decorator(orig_fitness_function):
                base_fitness = base_fitness_decorator(orig_fitness_function)

                def new_fitness_function(problem, dv):
                    fitness = base_fitness(problem, dv)
                    total_tof_days = sum(direct_leg_tofs(dv, n_legs)) * self._Edy2Kdy
                    if total_tof_days < tof_low:
                        fitness += float(tof_window_penalty) * (tof_low - total_tof_days)
                    elif total_tof_days > tof_high:
                        fitness += float(tof_window_penalty) * (total_tof_days - tof_high)
                    return fitness

                return new_fitness_function

            decorated_problem = pg.problem(pg.decorator_problem(udp, fitness_decorator=fitness_decorator))
            population = self._seeded_population(
                decorated_problem,
                direct_gene,
                pop_size,
                sigma_fraction=sigma_fraction,
                rng=rng,
            )
            algorithm = pg.algorithm(pg.sade(gen=int(sade_gen)))
            population = algorithm.evolve(population)
            champion_gene = list(population.champion_x)
            fitness = float(population.champion_f[0])
            try:
                result_dv, result_t0, result_tof, ejection_vinf = decode_dV_tof(
                    udp,
                    champion_gene,
                    planet_pack=self.planet_pack,
                )
                decode_error = None
            except (ValueError, ZeroDivisionError, FloatingPointError) as exc:
                result_dv = None
                result_t0 = cell_center_t0
                result_tof = cell_center_tof
                ejection_vinf = None
                decode_error = str(exc)
            row = {
                "t0": cell_center_t0,
                "tof": cell_center_tof,
                "metric": result_dv,
                "fitness": fitness,
                "ejection_vinf": ejection_vinf,
                "gene": champion_gene,
                "decode_error": decode_error,
                "cell_days": cell_days,
                "cell_t0_low": t0_low,
                "cell_t0_high": t0_high,
                "cell_tof_low": tof_low,
                "cell_tof_high": tof_high,
            }
            rows.append(row)
            return row

        def refine_cell(t0_low, t0_high, tof_low, tof_high, cell_days):
            row = evaluate_cell(t0_low, t0_high, tof_low, tof_high, cell_days)
            if row is None:
                return
            if (
                row.get("metric") is not None
                and float(row["metric"]) <= float(refine_threshold_dv)
                and cell_days / 2 >= float(min_cell_days)
                and evaluated < int(max_cells)
            ):
                t0_mid = (t0_low + t0_high) / 2
                tof_mid = (tof_low + tof_high) / 2
                for q_t0_low, q_t0_high, q_tof_low, q_tof_high in [
                    (t0_low, t0_mid, tof_low, tof_mid),
                    (t0_mid, t0_high, tof_low, tof_mid),
                    (t0_low, t0_mid, tof_mid, tof_high),
                    (t0_mid, t0_high, tof_mid, tof_high),
                ]:
                    refine_cell(q_t0_low, q_t0_high, q_tof_low, q_tof_high, cell_days / 2)

        t0_min = float(center_t0) - float(span_t0_days) / 2
        t0_max = float(center_t0) + float(span_t0_days) / 2
        tof_min = float(center_tof) - float(span_tof_days) / 2
        tof_max = float(center_tof) + float(span_tof_days) / 2
        t0_edges = np.arange(t0_min, t0_max, float(initial_cell_days))
        tof_edges = np.arange(tof_min, tof_max, float(initial_cell_days))
        for t0_low in t0_edges:
            for tof_low in tof_edges:
                refine_cell(
                    t0_low,
                    min(t0_low + float(initial_cell_days), t0_max),
                    tof_low,
                    min(tof_low + float(initial_cell_days), tof_max),
                    float(initial_cell_days),
                )
                if evaluated >= int(max_cells):
                    break
            if evaluated >= int(max_cells):
                break

        store = SQLiteJobStore(db_path)
        try:
            store.replace_porkchop_samples(
                run_id,
                sampler_name,
                rows,
                metadata={
                    "sampler_type": "adaptive_direct_cell_refinement",
                    "source_run_id": int(run_id),
                    "center_t0": center_t0,
                    "center_tof": center_tof,
                    "span_t0_days": span_t0_days,
                    "span_tof_days": span_tof_days,
                    "initial_cell_days": initial_cell_days,
                    "min_cell_days": min_cell_days,
                    "refine_threshold_dv": refine_threshold_dv,
                    "alpha_seed_dv_margin": alpha_seed_dv_margin,
                    "alpha_seed_max": alpha_seed_max,
                    "seed_search_t0_scale": seed_search_t0_scale,
                    "seed_search_tof_scale": seed_search_tof_scale,
                    "t0_wiggle_factor": t0_wiggle_factor,
                    "leg_tof_wiggle_days": leg_tof_wiggle_days,
                    "tof_window_penalty": tof_window_penalty,
                    "sade_gen": sade_gen,
                    "pop_size": pop_size,
                    "sigma_fraction": sigma_fraction,
                    "random_seed": random_seed,
                    "max_cells": max_cells,
                    "evaluated_cells": evaluated,
                    "sample_count": len(rows),
                },
            )
            samples = store.porkchop_samples(run_id, sampler_name)
        finally:
            store.close()
        result = pd.DataFrame(samples)
        if not result.empty:
            meta = pd.DataFrame(rows)[[
                "t0",
                "tof",
                "cell_days",
                "cell_t0_low",
                "cell_t0_high",
                "cell_tof_low",
                "cell_tof_high",
            ]]
            result = result.merge(
                meta,
                left_on=["result_t0", "result_tof"],
                right_on=["t0", "tof"],
                how="left",
            )
        return result

    def refine_local_porkchop_sqlite(
        self,
        db_path,
        run_id,
        sampler_name="direct_adaptive_cells",
        plot_styles=None,
        output_dir=None,
        plot_kwargs=None,
        **refinement_kwargs,
    ):
        """Run adaptive direct refinement and optionally render porkchop plots."""
        samples = self.adaptive_direct_cell_refinement_sqlite(
            db_path,
            run_id,
            sampler_name=sampler_name,
            **refinement_kwargs,
        )
        plots = {}
        if plot_styles is not None:
            plot_kwargs = {} if plot_kwargs is None else dict(plot_kwargs)
            for style in plot_styles:
                output_path = None
                if output_dir is not None:
                    output_path = str(Path(output_dir) / f"{sampler_name}_{style}.png")
                self.plot_adaptive_binned_sampled_porkchop_sqlite(
                    db_path,
                    run_id,
                    sampler_name=sampler_name,
                    style=style,
                    output_path=output_path,
                    **plot_kwargs,
                )
                plots[str(style)] = output_path
        return {
            "sampler_name": sampler_name,
            "samples": samples,
            "plots": plots,
        }


    def best_known_sqlite(
        self,
        db_path=None,
        batch_name=None,
        start_body=None,
        target_body=None,
        sequence_short_name=None,
        t0_range=None,
        contains_flyby=None,
        limit=10,
        include_benchmarks=False,
    ):
        """Return best known SQLite results across batches and binning schemes."""
        from _SQLiteStore import SQLiteJobStore

        if db_path is None:
            db_path = self._datastore_name + ".sqlite"
        store = SQLiteJobStore(db_path)
        try:
            return store.best_results(
                planet_pack=self.planet_pack,
                batch_name=batch_name,
                start_body=start_body,
                target_body=target_body,
                sequence_short_name=sequence_short_name,
                t0_range=t0_range,
                contains_flyby=contains_flyby,
                limit=limit,
                include_benchmarks=include_benchmarks,
            )
        finally:
            store.close()

    def find_best_known_plan_sqlite(
        self,
        db_path=None,
        batch_name=None,
        start_body=None,
        target_body=None,
        sequence_short_name=None,
        t0_range=None,
        contains_flyby=None,
        include_benchmarks=False,
    ):
        """Find and display the best known SQLite result as a flight plan."""
        rows = self.best_known_sqlite(
            db_path=db_path,
            batch_name=batch_name,
            start_body=start_body,
            target_body=target_body,
            sequence_short_name=sequence_short_name,
            t0_range=t0_range,
            contains_flyby=contains_flyby,
            limit=1,
            include_benchmarks=include_benchmarks,
        )
        if not rows:
            logger.warning("No matching DONE result found in SQLite datastore")
            return None

        best = rows[0]
        gene = json.loads(best["gene_json"])
        udp = self._mga_problem_from_sqlite_context(best)
        logger.info("Best dV is %.1f m/s", round(best["objective_dv"], 1))
        transx(udp, gene, planet_pack=self.planet_pack)
        return best

    def porkchop_points_from_snapshots_sqlite(self, db_path, run_id, include_final_population=True):
        """Decode stored optimizer genes as porkchop-ready points."""
        from _SQLiteStore import SQLiteJobStore

        store = SQLiteJobStore(db_path)
        try:
            context = store.run_context(run_id)
            if context is None:
                raise ValueError("Unknown run_id: " + str(run_id))
            udp = self._mga_problem_from_sqlite_context(context)
            points = [
                dict(point, source="snapshot_champion", individual_index=None)
                for point in store.optimizer_snapshot_points(run_id)
            ]
            if include_final_population:
                points.extend(store.optimizer_population_points(run_id, source="final"))
        finally:
            store.close()

        rows = []
        for point in points:
            try:
                result_dv, result_t0, result_tof, ejection_vinf = decode_dV_tof(
                    udp,
                    point["gene"],
                    planet_pack=self.planet_pack,
                )
                decode_error = None
            except (ValueError, ZeroDivisionError, FloatingPointError) as exc:
                result_dv = np.nan
                result_t0 = np.nan
                result_tof = np.nan
                ejection_vinf = np.nan
                decode_error = str(exc)
            rows.append({
                "run_id": point["run_id"],
                "step": point["step"],
                "island_index": point["island_index"],
                "individual_index": point["individual_index"],
                "source": point["source"],
                "fitness": point["fitness"],
                "best_fitness": point.get("best_fitness"),
                "result_DV": result_dv,
                "result_t0": result_t0,
                "result_tof": result_tof,
                "ejection_vinf": ejection_vinf,
                "decode_error": decode_error,
                "gene": point["gene"],
            })
        return pd.DataFrame(rows)

    optimizer_porkchop_points_sqlite = porkchop_points_from_snapshots_sqlite

    def plot_porkchop_from_snapshots_sqlite(
        self,
        db_path,
        run_id,
        include_final_population=True,
        metric="result_DV",
        figsize=(9, 7),
        cmap="wayfinder_lwp",
        annotate_best=True,
        show_samples=True,
    ):
        """Plot a local porkchop from stored optimizer snapshot/population genes."""
        points = self.porkchop_points_from_snapshots_sqlite(
            db_path,
            run_id,
            include_final_population=include_final_population,
        )
        if points.empty:
            logger.warning("No optimizer genes found for porkchop plot")
            return None
        if metric not in points.columns:
            raise ValueError("Unknown porkchop metric: " + str(metric))
        valid_points = points.dropna(subset=["result_t0", "result_tof", metric])
        if valid_points.empty:
            logger.warning("No valid optimizer genes found for porkchop plot")
            return points

        fig, ax = plt.subplots(figsize=figsize)
        unique_points = valid_points.drop_duplicates(subset=["result_t0", "result_tof"])
        colormap = _porkchop_colormap(cmap)
        surface = None
        if len(unique_points) >= 3:
            surface = ax.tricontourf(
                unique_points["result_t0"],
                unique_points["result_tof"],
                unique_points[metric],
                levels=24,
                cmap=colormap,
            )
            ax.tricontour(
                unique_points["result_t0"],
                unique_points["result_tof"],
                unique_points[metric],
                levels=12,
                colors="black",
                linewidths=0.35,
                alpha=0.35,
            )
        else:
            surface = ax.scatter(
                valid_points["result_t0"],
                valid_points["result_tof"],
                c=valid_points[metric],
                cmap=colormap,
                alpha=0.85,
                edgecolors="none",
                label="sampled solutions",
            )

        if show_samples:
            ax.scatter(
                valid_points["result_t0"],
                valid_points["result_tof"],
                marker="o",
                s=26,
                facecolors="none",
                edgecolors="white",
                linewidths=0.55,
                alpha=0.55,
                label="sampled solutions",
            )

        best_idx = valid_points[metric].idxmin()
        best = valid_points.loc[best_idx]
        ax.scatter(
            [best["result_t0"]],
            [best["result_tof"]],
            marker="*",
            s=160,
            color="black",
            label="best " + metric,
            zorder=5,
        )
        if annotate_best:
            ax.annotate(
                f"{best[metric]:.1f}",
                (best["result_t0"], best["result_tof"]),
                textcoords="offset points",
                xytext=(8, 8),
            )

        ax.set_xlabel("T0 (KSP days)")
        ax.set_ylabel("TOF (KSP days)")
        ax.set_title("Optimizer local porkchop")
        ax.legend(loc="best")
        if surface is not None:
            cbar = fig.colorbar(surface, ax=ax)
            cbar.set_label(metric)
        fig.tight_layout()
        return points

    plot_optimizer_porkchop_sqlite = plot_porkchop_from_snapshots_sqlite

    def plot_augmented_optimizer_porkchop_sqlite(
        self,
        db_path,
        run_id,
        refinement_sampler_names,
        include_final_population=True,
        metric="result_DV",
        figsize=(10, 8),
        background_cmap="wayfinder_lwp",
        refinement_cmap="wayfinder_lwp",
        background_levels=24,
        refinement_clip_percentile=(0, 95),
        show_optimizer_points=True,
        show_refinement_points=True,
        show_refinement_best_per_bin=True,
        refinement_bin_days=0.5,
    ):
        """Plot optimizer density augmented by direct local refinement samples."""
        optimizer_points = self.porkchop_points_from_snapshots_sqlite(
            db_path,
            run_id,
            include_final_population=include_final_population,
        )
        refinement_points = self.sampled_porkchop_points_sqlite(
            db_path,
            run_id,
            sampler_name=refinement_sampler_names,
        )
        if optimizer_points.empty:
            logger.warning("No optimizer genes found for augmented porkchop plot")
            return None
        optimizer_valid = optimizer_points.dropna(subset=["result_t0", "result_tof", metric])
        refinement_valid = refinement_points.dropna(subset=["result_t0", "result_tof", metric])
        if optimizer_valid.empty:
            logger.warning("No valid optimizer genes found for augmented porkchop plot")
            return optimizer_points, refinement_points

        fig, ax = plt.subplots(figsize=figsize)
        optimizer_unique = optimizer_valid.drop_duplicates(subset=["result_t0", "result_tof"])
        background_colormap = _porkchop_colormap(background_cmap)
        surface = None
        if len(optimizer_unique) >= 3:
            surface = ax.tricontourf(
                optimizer_unique["result_t0"],
                optimizer_unique["result_tof"],
                optimizer_unique[metric],
                levels=background_levels,
                cmap=background_colormap,
                alpha=0.82,
            )
            ax.tricontour(
                optimizer_unique["result_t0"],
                optimizer_unique["result_tof"],
                optimizer_unique[metric],
                levels=max(8, background_levels // 2),
                colors="black",
                linewidths=0.30,
                alpha=0.30,
            )
        if show_optimizer_points:
            ax.scatter(
                optimizer_valid["result_t0"],
                optimizer_valid["result_tof"],
                marker="o",
                s=23,
                facecolors="none",
                edgecolors="white",
                linewidths=0.55,
                alpha=0.52,
                label="optimizer samples",
            )

        refinement_scatter = None
        if not refinement_valid.empty:
            if refinement_clip_percentile is None:
                vmin = refinement_valid[metric].min()
                vmax = refinement_valid[metric].max()
            else:
                vmin = refinement_valid[metric].quantile(float(refinement_clip_percentile[0]) / 100.0)
                vmax = refinement_valid[metric].quantile(float(refinement_clip_percentile[1]) / 100.0)
            refinement_colormap = _porkchop_colormap(refinement_cmap, over_color="red")
            refinement_norm = colors.Normalize(vmin=vmin, vmax=vmax, clip=False)
            if show_refinement_points:
                refinement_scatter = ax.scatter(
                    refinement_valid["result_t0"],
                    refinement_valid["result_tof"],
                    c=refinement_valid[metric],
                    cmap=refinement_colormap,
                    norm=refinement_norm,
                    marker="o",
                    s=30,
                    alpha=0.80,
                    edgecolors="black",
                    linewidths=0.25,
                    label="direct refinements",
                    zorder=4,
                )
            if show_refinement_best_per_bin:
                binned = self.binned_sampled_porkchop_points_sqlite(
                    db_path,
                    run_id,
                    sampler_name=refinement_sampler_names,
                    metric=metric,
                    t0_bin_days=refinement_bin_days,
                    tof_bin_days=refinement_bin_days,
                )
                if not binned.empty:
                    ax.scatter(
                        binned["bin_t0_center"],
                        binned["bin_tof_center"],
                        c=binned[metric],
                        cmap=refinement_colormap,
                        norm=refinement_norm,
                        marker="D",
                        s=48,
                        alpha=0.95,
                        edgecolors="black",
                        linewidths=0.45,
                        label=f"best / {refinement_bin_days:g}d bin",
                        zorder=5,
                    )

        combined = optimizer_valid
        if not refinement_valid.empty:
            combined = pd.concat([optimizer_valid, refinement_valid], ignore_index=True)
        best_idx = combined[metric].idxmin()
        best = combined.loc[best_idx]
        ax.scatter(
            [best["result_t0"]],
            [best["result_tof"]],
            marker="*",
            s=180,
            color="black",
            label="best " + metric,
            zorder=6,
        )
        ax.annotate(
            f"{best[metric]:.1f}",
            (best["result_t0"], best["result_tof"]),
            textcoords="offset points",
            xytext=(8, 8),
        )

        ax.set_xlabel("T0 (KSP days)")
        ax.set_ylabel("TOF (KSP days)")
        ax.set_title("Optimizer porkchop augmented with direct refinements")
        ax.legend(loc="best")
        if surface is not None:
            cbar = fig.colorbar(surface, ax=ax)
            cbar.set_label("optimizer " + metric)
        if refinement_scatter is not None:
            cbar_refined = fig.colorbar(refinement_scatter, ax=ax, fraction=0.046, pad=0.08)
            cbar_refined.set_label("direct " + metric)
        fig.tight_layout()
        return optimizer_points, refinement_points

    def sample_local_porkchop_sqlite(
        self,
        db_path,
        run_id,
        sampler_name="local_grid",
        t0_span=120,
        tof_span=220,
        n_t0=45,
        n_tof=45,
        seed_gene=None,
    ):
        """Sample a uniform local T0/TOF porkchop around a known solution."""
        from _SQLiteStore import SQLiteJobStore

        store = SQLiteJobStore(db_path)
        try:
            context = store.run_context(run_id)
            if context is None:
                raise ValueError("Unknown run_id: " + str(run_id))
            if seed_gene is None:
                if not context.get("gene_json"):
                    raise ValueError("Run has no stored seed gene: " + str(run_id))
                seed_gene = json.loads(context["gene_json"])
            udp = self._mga_problem_from_sqlite_context(context)
            problem = pg.problem(udp)
            lower_bounds, upper_bounds = problem.get_bounds()

            center_t0 = float(context["result_t0"])
            center_tof = float(context["result_tof"])
            t0_values = np.linspace(
                max(float(context["t0_min"]), center_t0 - t0_span / 2),
                min(float(context["t0_max"]), center_t0 + t0_span / 2),
                int(n_t0),
            )
            tof_values = np.linspace(
                max(float(context["tof_min"]), center_tof - tof_span / 2),
                min(float(context["tof_max"]), center_tof + tof_span / 2),
                int(n_tof),
            )

            rows = []
            for t0 in t0_values:
                for tof in tof_values:
                    gene = list(seed_gene)
                    gene[0] = float(t0) / self._Edy2Kdy
                    gene[-1] = float(tof) / self._Edy2Kdy
                    gene = np.clip(np.array(gene, dtype=float), lower_bounds, upper_bounds).tolist()
                    try:
                        fitness = float(problem.fitness(gene)[0])
                        result_dv, result_t0, result_tof, ejection_vinf = decode_dV_tof(
                            udp,
                            gene,
                            planet_pack=self.planet_pack,
                        )
                        decode_error = None
                    except (ValueError, ZeroDivisionError, FloatingPointError) as exc:
                        fitness = None
                        result_dv = None
                        result_t0 = float(t0)
                        result_tof = float(tof)
                        ejection_vinf = None
                        decode_error = str(exc)
                    rows.append({
                        "t0": result_t0,
                        "tof": result_tof,
                        "metric": result_dv,
                        "fitness": fitness,
                        "ejection_vinf": ejection_vinf,
                        "gene": gene,
                        "decode_error": decode_error,
                    })

            store.replace_porkchop_samples(run_id, sampler_name, rows)
            samples = store.porkchop_samples(run_id, sampler_name)
        finally:
            store.close()
        return pd.DataFrame(samples)

    def sampled_porkchop_points_sqlite(self, db_path, run_id, sampler_name="local_grid"):
        from _SQLiteStore import SQLiteJobStore

        store = SQLiteJobStore(db_path)
        try:
            return pd.DataFrame(store.porkchop_samples(run_id, sampler_name))
        finally:
            store.close()

    def binned_sampled_porkchop_points_sqlite(
        self,
        db_path,
        run_id,
        sampler_name="local_grid",
        metric="result_DV",
        t0_bin_days=0.5,
        tof_bin_days=0.5,
    ):
        points = self.sampled_porkchop_points_sqlite(db_path, run_id, sampler_name=sampler_name)
        if points.empty:
            return points
        valid_points = points.dropna(subset=["result_t0", "result_tof", metric]).copy()
        if valid_points.empty:
            return valid_points

        valid_points["t0_bin"] = np.floor(valid_points["result_t0"] / float(t0_bin_days)) * float(t0_bin_days)
        valid_points["tof_bin"] = np.floor(valid_points["result_tof"] / float(tof_bin_days)) * float(tof_bin_days)
        idx = valid_points.groupby(["t0_bin", "tof_bin"])[metric].idxmin()
        binned = valid_points.loc[idx].copy()
        binned["bin_t0_center"] = binned["t0_bin"] + float(t0_bin_days) / 2
        binned["bin_tof_center"] = binned["tof_bin"] + float(tof_bin_days) / 2
        binned["bin_t0_days"] = float(t0_bin_days)
        binned["bin_tof_days"] = float(tof_bin_days)
        return binned.sort_values(["bin_t0_center", "bin_tof_center"])

    def adaptive_binned_sampled_porkchop_points_sqlite(
        self,
        db_path,
        run_id,
        sampler_name="local_grid",
        metric="result_DV",
        coarse_bin_days=2.0,
        min_bin_days=0.5,
        min_points_to_split=8,
        t0_range=None,
        tof_range=None,
    ):
        points = self.sampled_porkchop_points_sqlite(db_path, run_id, sampler_name=sampler_name)
        if points.empty:
            return points
        valid_points = points.dropna(subset=["result_t0", "result_tof", metric]).copy()
        if valid_points.empty:
            return valid_points
        if t0_range is not None:
            valid_points = valid_points[
                (valid_points["result_t0"] >= float(t0_range[0]))
                & (valid_points["result_t0"] <= float(t0_range[1]))
            ].copy()
        if tof_range is not None:
            valid_points = valid_points[
                (valid_points["result_tof"] >= float(tof_range[0]))
                & (valid_points["result_tof"] <= float(tof_range[1]))
            ].copy()
        if valid_points.empty:
            return valid_points

        coarse_bin_days = float(coarse_bin_days)
        min_bin_days = float(min_bin_days)
        t0_min = np.floor(valid_points["result_t0"].min() / coarse_bin_days) * coarse_bin_days
        t0_max = np.ceil(valid_points["result_t0"].max() / coarse_bin_days) * coarse_bin_days
        tof_min = np.floor(valid_points["result_tof"].min() / coarse_bin_days) * coarse_bin_days
        tof_max = np.ceil(valid_points["result_tof"].max() / coarse_bin_days) * coarse_bin_days

        leaves = []

        def split_cell(cell_points, t0_low, t0_high, tof_low, tof_high):
            t0_width = t0_high - t0_low
            tof_width = tof_high - tof_low
            can_split = (
                len(cell_points) >= int(min_points_to_split)
                and t0_width / 2 >= min_bin_days
                and tof_width / 2 >= min_bin_days
            )
            if not can_split:
                if not cell_points.empty:
                    best = cell_points.loc[cell_points[metric].idxmin()].copy()
                    best["bin_t0"] = t0_low
                    best["bin_tof"] = tof_low
                    best["bin_t0_center"] = (t0_low + t0_high) / 2
                    best["bin_tof_center"] = (tof_low + tof_high) / 2
                    best["bin_t0_days"] = t0_width
                    best["bin_tof_days"] = tof_width
                    best["bin_point_count"] = len(cell_points)
                    leaves.append(best)
                return

            t0_mid = (t0_low + t0_high) / 2
            tof_mid = (tof_low + tof_high) / 2
            quadrants = [
                (t0_low, t0_mid, tof_low, tof_mid),
                (t0_mid, t0_high, tof_low, tof_mid),
                (t0_low, t0_mid, tof_mid, tof_high),
                (t0_mid, t0_high, tof_mid, tof_high),
            ]
            for q_t0_low, q_t0_high, q_tof_low, q_tof_high in quadrants:
                q_points = cell_points[
                    (cell_points["result_t0"] >= q_t0_low)
                    & (cell_points["result_t0"] < q_t0_high)
                    & (cell_points["result_tof"] >= q_tof_low)
                    & (cell_points["result_tof"] < q_tof_high)
                ]
                split_cell(q_points, q_t0_low, q_t0_high, q_tof_low, q_tof_high)

        t0_edges = np.arange(t0_min, t0_max + coarse_bin_days, coarse_bin_days)
        tof_edges = np.arange(tof_min, tof_max + coarse_bin_days, coarse_bin_days)
        for t0_low, t0_high in zip(t0_edges[:-1], t0_edges[1:]):
            for tof_low, tof_high in zip(tof_edges[:-1], tof_edges[1:]):
                cell_points = valid_points[
                    (valid_points["result_t0"] >= t0_low)
                    & (valid_points["result_t0"] < t0_high)
                    & (valid_points["result_tof"] >= tof_low)
                    & (valid_points["result_tof"] < tof_high)
                ]
                split_cell(cell_points, t0_low, t0_high, tof_low, tof_high)

        if not leaves:
            return valid_points.iloc[0:0].copy()
        return pd.DataFrame(leaves).sort_values(["bin_t0_center", "bin_tof_center"])

    def plot_sampled_porkchop_sqlite(
        self,
        db_path,
        run_id,
        sampler_name="local_grid",
        metric="result_DV",
        figsize=(9, 7),
        cmap="wayfinder_lwp",
        show_samples=False,
        clip_percentile=None,
        clipped_color="red",
        level_mode="linear_log",
        level_floor=None,
        level_floor_round=50,
        level_count=28,
        linear_factor=2.0,
        linear_log_color_style="linear_to_red",
        show_level_split=True,
        split_color="magenta",
    ):
        points = self.sampled_porkchop_points_sqlite(db_path, run_id, sampler_name=sampler_name)
        if points.empty:
            logger.warning("No sampled porkchop points found")
            return None
        valid_points = points.dropna(subset=["result_t0", "result_tof", metric])
        if valid_points.empty:
            logger.warning("No valid sampled porkchop points found")
            return points

        fig, ax = plt.subplots(figsize=figsize)
        unique_points = valid_points.drop_duplicates(subset=["result_t0", "result_tof"])
        levels = _porkchop_levels(
            valid_points[metric],
            mode=level_mode,
            floor=level_floor,
            ceiling=valid_points[metric].quantile(float(clip_percentile[1]) / 100.0)
            if clip_percentile is not None
            else None,
            floor_round=level_floor_round,
            count=level_count,
            linear_factor=linear_factor,
        )
        transition = min(levels[0] * float(linear_factor), levels[-1])
        if level_mode == "linear_log" and cmap == "wayfinder_lwp":
            colormap, norm = _porkchop_boundary_colormap(
                levels,
                mode=level_mode,
                transition=transition,
                over_color=clipped_color,
                style=linear_log_color_style,
            )
        else:
            colormap = _porkchop_colormap(cmap, over_color=clipped_color)
            norm = colors.BoundaryNorm(levels, colormap.N)
        surface = ax.tricontourf(
            unique_points["result_t0"],
            unique_points["result_tof"],
            unique_points[metric],
            levels=levels,
            cmap=colormap,
            norm=norm,
            extend="both",
        )
        ax.tricontour(
            unique_points["result_t0"],
            unique_points["result_tof"],
            unique_points[metric],
            levels=levels,
            colors="black",
            linewidths=0.35,
            alpha=0.35,
        )
        if show_samples:
            ax.scatter(
                valid_points["result_t0"],
                valid_points["result_tof"],
                marker="o",
                s=12,
                facecolors="none",
                edgecolors="white",
                linewidths=0.35,
                alpha=0.35,
                label="sampled grid",
            )
        best_idx = valid_points[metric].idxmin()
        best = valid_points.loc[best_idx]
        ax.scatter(
            [best["result_t0"]],
            [best["result_tof"]],
            marker="*",
            s=160,
            color="black",
            label="best sampled " + metric,
            zorder=5,
        )
        ax.annotate(
            f"{best[metric]:.1f}",
            (best["result_t0"], best["result_tof"]),
            textcoords="offset points",
            xytext=(8, 8),
        )
        ax.set_xlabel("T0 (KSP days)")
        ax.set_ylabel("TOF (KSP days)")
        ax.set_title("Sampled local porkchop")
        ax.legend(loc="best")
        cbar = fig.colorbar(surface, ax=ax)
        cbar.set_label(metric)
        if show_level_split and level_mode == "linear_log" and levels[0] < transition < levels[-1]:
            cbar.ax.axhline(transition, color=split_color, linewidth=1.6)
            cbar.ax.text(
                1.08,
                transition,
                "linear/log",
                color=split_color,
                fontsize=8,
                va="center",
                ha="left",
                transform=cbar.ax.get_yaxis_transform(),
            )
        fig.tight_layout()
        return points

    def plot_sampled_porkchop_scatter_sqlite(
        self,
        db_path,
        run_id,
        sampler_name="local_grid",
        metric="result_DV",
        figsize=(9, 7),
        cmap="wayfinder_lwp",
        clip_percentile=None,
        point_size=24,
        alpha=0.8,
    ):
        points = self.sampled_porkchop_points_sqlite(db_path, run_id, sampler_name=sampler_name)
        if points.empty:
            logger.warning("No sampled porkchop points found")
            return None
        valid_points = points.dropna(subset=["result_t0", "result_tof", metric])
        if valid_points.empty:
            logger.warning("No valid sampled porkchop points found")
            return points

        if clip_percentile is None:
            vmin = valid_points[metric].min()
            vmax = valid_points[metric].max()
        else:
            vmin = valid_points[metric].quantile(float(clip_percentile[0]) / 100.0)
            vmax = valid_points[metric].quantile(float(clip_percentile[1]) / 100.0)
        colormap = _porkchop_colormap(cmap, over_color="red")
        norm = colors.Normalize(vmin=vmin, vmax=vmax, clip=False)

        fig, ax = plt.subplots(figsize=figsize)
        scatter = ax.scatter(
            valid_points["result_t0"],
            valid_points["result_tof"],
            c=valid_points[metric],
            cmap=colormap,
            norm=norm,
            s=point_size,
            alpha=alpha,
            edgecolors="black",
            linewidths=0.2,
        )
        best_idx = valid_points[metric].idxmin()
        best = valid_points.loc[best_idx]
        ax.scatter(
            [best["result_t0"]],
            [best["result_tof"]],
            marker="*",
            s=180,
            color="black",
            label="best " + metric,
            zorder=5,
        )
        ax.annotate(
            f"{best[metric]:.1f}",
            (best["result_t0"], best["result_tof"]),
            textcoords="offset points",
            xytext=(8, 8),
        )
        ax.set_xlabel("T0 (KSP days)")
        ax.set_ylabel("TOF (KSP days)")
        ax.set_title("Sampled local porkchop scatter")
        ax.legend(loc="best")
        cbar = fig.colorbar(scatter, ax=ax)
        cbar.set_label(metric)
        fig.tight_layout()
        return points

    def plot_binned_sampled_porkchop_sqlite(
        self,
        db_path,
        run_id,
        sampler_name="local_grid",
        metric="result_DV",
        t0_bin_days=0.5,
        tof_bin_days=0.5,
        figsize=(9, 7),
        cmap="wayfinder_lwp",
        clip_percentile=None,
        point_size=64,
        show_raw=True,
    ):
        raw_points = self.sampled_porkchop_points_sqlite(db_path, run_id, sampler_name=sampler_name)
        binned_points = self.binned_sampled_porkchop_points_sqlite(
            db_path,
            run_id,
            sampler_name=sampler_name,
            metric=metric,
            t0_bin_days=t0_bin_days,
            tof_bin_days=tof_bin_days,
        )
        if binned_points.empty:
            logger.warning("No binned sampled porkchop points found")
            return binned_points

        if clip_percentile is None:
            vmin = binned_points[metric].min()
            vmax = binned_points[metric].max()
        else:
            vmin = binned_points[metric].quantile(float(clip_percentile[0]) / 100.0)
            vmax = binned_points[metric].quantile(float(clip_percentile[1]) / 100.0)
        colormap = _porkchop_colormap(cmap, over_color="red")
        norm = colors.Normalize(vmin=vmin, vmax=vmax, clip=False)

        fig, ax = plt.subplots(figsize=figsize)
        if show_raw and not raw_points.empty:
            raw_valid = raw_points.dropna(subset=["result_t0", "result_tof", metric])
            ax.scatter(
                raw_valid["result_t0"],
                raw_valid["result_tof"],
                s=12,
                color="0.70",
                alpha=0.25,
                edgecolors="none",
                label="raw samples",
            )

        scatter = ax.scatter(
            binned_points["bin_t0_center"],
            binned_points["bin_tof_center"],
            c=binned_points[metric],
            cmap=colormap,
            norm=norm,
            s=point_size,
            alpha=0.95,
            edgecolors="black",
            linewidths=0.35,
            label="best per bin",
        )
        best_idx = binned_points[metric].idxmin()
        best = binned_points.loc[best_idx]
        ax.scatter(
            [best["bin_t0_center"]],
            [best["bin_tof_center"]],
            marker="*",
            s=180,
            color="black",
            label="best " + metric,
            zorder=5,
        )
        ax.annotate(
            f"{best[metric]:.1f}",
            (best["bin_t0_center"], best["bin_tof_center"]),
            textcoords="offset points",
            xytext=(8, 8),
        )
        ax.set_xlabel("T0 (KSP days)")
        ax.set_ylabel("TOF (KSP days)")
        ax.set_title(f"Binned local porkchop ({t0_bin_days:g} x {tof_bin_days:g} days)")
        ax.legend(loc="best")
        cbar = fig.colorbar(scatter, ax=ax)
        cbar.set_label(metric)
        fig.tight_layout()
        return binned_points

    def plot_adaptive_binned_sampled_porkchop_sqlite(
        self,
        db_path,
        run_id,
        sampler_name="local_grid",
        metric="result_DV",
        coarse_bin_days=2.0,
        min_bin_days=0.5,
        min_points_to_split=8,
        figsize=(9, 7),
        cmap="wayfinder_lwp",
        clip_percentile=None,
        vmin=None,
        vmax=None,
        show_raw=True,
        t0_range=None,
        tof_range=None,
        title=None,
        output_path=None,
        style="cells",
        continuous_levels=32,
        color_levels=None,
        low_detail_step=100,
        low_detail_factor=1.25,
        coarse_detail_step=500,
        raw_point_size=12,
        raw_point_alpha=0.55,
        raw_point_facecolor="none",
        raw_point_edgecolor="white",
    ):
        raw_points = self.sampled_porkchop_points_sqlite(db_path, run_id, sampler_name=sampler_name)
        binned_points = self.adaptive_binned_sampled_porkchop_points_sqlite(
            db_path,
            run_id,
            sampler_name=sampler_name,
            metric=metric,
            coarse_bin_days=coarse_bin_days,
            min_bin_days=min_bin_days,
            min_points_to_split=min_points_to_split,
            t0_range=t0_range,
            tof_range=tof_range,
        )
        if binned_points.empty:
            logger.warning("No adaptive binned sampled porkchop points found")
            return binned_points

        if vmin is None:
            if clip_percentile is None:
                vmin = binned_points[metric].min()
            else:
                vmin = binned_points[metric].quantile(float(clip_percentile[0]) / 100.0)
        if isinstance(vmax, str):
            if vmax != "double_floor":
                raise ValueError("Unknown adaptive porkchop vmax mode: " + vmax)
            vmax = _rounded_porkchop_floor(binned_points[metric].min(), 50) * 2.0
        if vmax is None:
            if clip_percentile is None:
                vmax = binned_points[metric].max()
            else:
                vmax = binned_points[metric].quantile(float(clip_percentile[1]) / 100.0)
        colormap = _porkchop_colormap(cmap, over_color="red")
        levels = None
        if isinstance(color_levels, str):
            if color_levels != "low_detail":
                raise ValueError("Unknown adaptive porkchop color_levels mode: " + color_levels)
            low_ceiling = min(float(vmax), _rounded_porkchop_floor(float(vmin) * float(low_detail_factor), 50))
            fine_levels = np.arange(float(vmin), low_ceiling + float(low_detail_step), float(low_detail_step))
            coarse_levels = np.arange(
                low_ceiling + float(coarse_detail_step),
                float(vmax) + float(coarse_detail_step),
                float(coarse_detail_step),
            )
            levels = np.unique(np.concatenate([fine_levels, coarse_levels, [float(vmax)]]))
        elif color_levels is not None:
            levels = np.asarray(color_levels, dtype=float)
        if levels is not None:
            levels = levels[(levels >= float(vmin)) & (levels <= float(vmax))]
            levels = np.unique(np.concatenate([[float(vmin)], levels, [float(vmax)]]))
            if len(levels) < 2:
                raise ValueError("Adaptive porkchop color_levels must contain at least 2 levels")
            norm = colors.BoundaryNorm(levels, colormap.N)
        else:
            norm = colors.Normalize(vmin=vmin, vmax=vmax, clip=False)

        fig, ax = plt.subplots(figsize=figsize)
        raw_valid = pd.DataFrame()
        if show_raw and not raw_points.empty:
            raw_valid = raw_points.dropna(subset=["result_t0", "result_tof", metric])
            if t0_range is not None:
                raw_valid = raw_valid[
                    (raw_valid["result_t0"] >= float(t0_range[0]))
                    & (raw_valid["result_t0"] <= float(t0_range[1]))
                ]
            if tof_range is not None:
                raw_valid = raw_valid[
                    (raw_valid["result_tof"] >= float(tof_range[0]))
                    & (raw_valid["result_tof"] <= float(tof_range[1]))
                ]

        style = str(style).lower()
        surface = None
        if style == "cells":
            for _, row in binned_points.iterrows():
                rectangle = Rectangle(
                    (row["bin_t0"], row["bin_tof"]),
                    row["bin_t0_days"],
                    row["bin_tof_days"],
                    facecolor=colormap(norm(row[metric])),
                    edgecolor="black",
                    linewidth=0.45,
                    alpha=0.88,
                )
                ax.add_patch(rectangle)
        elif style == "continuous":
            unique_points = binned_points.drop_duplicates(subset=["bin_t0_center", "bin_tof_center"])
            if len(unique_points) < 3:
                raise ValueError("Continuous adaptive porkchop plot needs at least 3 unique points")
            if levels is None:
                levels = np.linspace(float(vmin), float(vmax), int(continuous_levels))
            surface = ax.tricontourf(
                unique_points["bin_t0_center"],
                unique_points["bin_tof_center"],
                unique_points[metric],
                levels=levels,
                cmap=colormap,
                norm=norm,
                extend="both",
            )
            ax.tricontour(
                unique_points["bin_t0_center"],
                unique_points["bin_tof_center"],
                unique_points[metric],
                levels=levels,
                colors="black",
                linewidths=0.35,
                alpha=0.35,
            )
        else:
            raise ValueError("Unknown adaptive porkchop style: " + str(style))

        if show_raw and not raw_valid.empty:
            ax.scatter(
                raw_valid["result_t0"],
                raw_valid["result_tof"],
                s=raw_point_size,
                facecolors=raw_point_facecolor,
                edgecolors=raw_point_edgecolor,
                linewidths=0.75,
                alpha=raw_point_alpha,
                label="raw samples",
                zorder=4,
            )

        best_idx = binned_points[metric].idxmin()
        best = binned_points.loc[best_idx]
        ax.scatter(
            [best["bin_t0_center"]],
            [best["bin_tof_center"]],
            marker="*",
            s=180,
            color="black",
            label="best " + metric,
            zorder=5,
        )
        ax.annotate(
            f"{best[metric]:.1f}",
            (best["bin_t0_center"], best["bin_tof_center"]),
            textcoords="offset points",
            xytext=(8, 8),
        )
        ax.set_xlim(
            binned_points["bin_t0"].min(),
            (binned_points["bin_t0"] + binned_points["bin_t0_days"]).max(),
        )
        ax.set_ylim(
            binned_points["bin_tof"].min(),
            (binned_points["bin_tof"] + binned_points["bin_tof_days"]).max(),
        )
        ax.set_xlabel("T0 (KSP days)")
        ax.set_ylabel("TOF (KSP days)")
        if title is None:
            title = f"Adaptive binned local porkchop ({coarse_bin_days:g} -> {min_bin_days:g} days, {style})"
        ax.set_title(title)
        ax.legend(loc="best")
        scalar_mappable = cm.ScalarMappable(norm=norm, cmap=colormap)
        if surface is not None:
            scalar_mappable = surface
        cbar = fig.colorbar(scalar_mappable, ax=ax)
        cbar.set_label(metric)
        fig.tight_layout()
        if output_path is not None:
            fig.savefig(output_path, dpi=150)
        return binned_points

    def plot_DVvsT0_sqlite(
        self,
        db_path,
        batch_name=None,
        start_body=None,
        target_body=None,
        sequence_short_name=None,
        t0_range=None,
        contains_flyby=None,
        include_benchmarks=False,
    ):
        """Plot DV vs launch date from SQLite results across batches."""
        from _SQLiteStore import SQLiteJobStore

        store = SQLiteJobStore(db_path)
        try:
            rows = store.result_rows(
                planet_pack=self.planet_pack,
                batch_name=batch_name,
                start_body=start_body,
                target_body=target_body,
                sequence_short_name=sequence_short_name,
                t0_range=t0_range,
                contains_flyby=contains_flyby,
                include_benchmarks=include_benchmarks,
            )
        finally:
            store.close()

        if not rows:
            logger.warning("No DONE result found in SQLite datastore for plot")
            return None

        df = pd.DataFrame(rows)
        idx = df.groupby(["sequence_short_name", "t0_min"])["objective_dv"].idxmin()
        filtered_df = df.loc[idx].sort_values(["sequence_short_name", "result_t0"])
        plt.figure(figsize=(10, 8))
        sns.lineplot(
            x='result_t0',
            y='objective_dv',
            hue='sequence_short_name',
            marker='o',
            data=filtered_df,
        )
        plt.ylabel("result_DV")
        plt.legend(title='Sequence', bbox_to_anchor=(1.05, 1), loc='upper left')
        return filtered_df

    def plot_by_sequences_sqlite(
        self,
        db_path,
        batch_name=None,
        start_body=None,
        target_body=None,
        sequence_short_name=None,
        t0_range=None,
        contains_flyby=None,
        include_benchmarks=False,
    ):
        """Plot the best known DV per sequence from SQLite results."""
        from _SQLiteStore import SQLiteJobStore

        store = SQLiteJobStore(db_path)
        try:
            rows = store.result_rows(
                planet_pack=self.planet_pack,
                batch_name=batch_name,
                start_body=start_body,
                target_body=target_body,
                sequence_short_name=sequence_short_name,
                t0_range=t0_range,
                contains_flyby=contains_flyby,
                include_benchmarks=include_benchmarks,
            )
        finally:
            store.close()

        if not rows:
            logger.warning("No DONE result found in SQLite datastore for plot")
            return None

        df = pd.DataFrame(rows)
        idx = df.groupby("sequence_short_name")["objective_dv"].idxmin()
        best_by_sequence = df.loc[idx].sort_values("objective_dv")
        fig, ax = plt.subplots(figsize=(9,6))
        palette = colors.Normalize(
            vmin=best_by_sequence["result_tof"].min(),
            vmax=best_by_sequence["result_tof"].max(),
        )
        plt.bar(
            range(len(best_by_sequence)),
            best_by_sequence["objective_dv"],
            color=plt.get_cmap('coolwarm')(palette(best_by_sequence["result_tof"])),
            edgecolor='black',
            linewidth=1,
            align='center',
        )
        sm = plt.cm.ScalarMappable(cmap='coolwarm', norm=palette)
        sm._A = []
        ax.set_ylabel('Total DV')
        cbar = fig.colorbar(sm, ax=ax)
        if self.planet_pack == "JNSQ" :
            cbar.set_label('ToF in JNSQ Days', rotation=90)
        elif self.planet_pack == "Vanilla" :
            cbar.set_label('ToF in KSP Days', rotation=90)
        plt.xticks(
            range(len(best_by_sequence)),
            list(best_by_sequence["sequence_short_name"]),
            rotation='vertical',
        )
        fig.tight_layout()
        return best_by_sequence
    def debugPrint(self):
        raise NotImplementedError("debugPrint was tied to the legacy dataframe staging path.")

    def rp_target_ward(self,target,injection_altitude):
        '''This function is a ward to avoid injection below the safe altitude or in atmosphere'''
        '''It will display a warning, and set a higher injection altitude of 1000 km'''
        if target.radius+injection_altitude < target.safe_radius:            
            logger.warning(
                "Injection altitude below safe limits, injection altitude for target %s set to 1000 km",
                target.name,
            )
            return target.radius+1000000
        else :
            return target.radius+injection_altitude

    def generateSequences(self,swing_by_bodies=[["Kerbin"],["Eve","Duna"],["Kerbin"]]):
        '''
        Generates flyby sequences possibilities from start to end passing by the different combinations of bodies in the swing_by_bodies list
        The overall sequence length will be at minimum 3 and at most equal to the length parameter
        '''
        
        seqs = list(i for i in product(*swing_by_bodies) if ((len(set(i[1:-1])) >= len(i[1:-1])-1 or len(i)-i.count("*")) < 5 and len(set(i[0:3])) > 0))
        
            
        sequences = []
        for i in seqs:
            temp = []
            for j in i:
                if isinstance(j, list):
                    for k in j:
                        temp.append(k)
                else:
                    temp.append(j)
            sequences.append(temp)

        '''Cleanup the *'s              '''
        
        for sbb in range(len(swing_by_bodies)-2):
            for i in range(len(sequences)):
                try:
                    sequences[i].remove('*')
                except ValueError:
                    pass  # do nothing!

        '''Remove duplicates            '''
        
        new_k = []
        for elem in sequences:
            if elem not in new_k:
                new_k.append(elem)
        sequences = new_k
        
        return sequences
        
    def generateShortSequences(self,swing_by_bodies=[["Kerbin"],["Eve","Duna"],["Kerbin"]]):
        '''
        Generate the list in shorthand notation for the flyby sequences.      
        '''
        shortSequences = []
        sequences = self.generateSequences(swing_by_bodies)
        for seq in sequences :
            seq_name = ''
            for body in seq:
                seq_name += self._Body_abrev_dic[body]
            shortSequences.append(seq_name)
        return shortSequences
    
    def orbital_period(self,planet):
        JNSQ_Dy2s = 12*3600
        Vanilla_Dy2s = 6*3600
        period = planet.period(0)
        if self.planet_pack == "JNSQ":
            return round(period / JNSQ_Dy2s,0)
        elif self.planet_pack == "Vanilla":
            return round(period / Vanilla_Dy2s,0)            
        else :
            logger.warning("Unknown planet pack, auto-tof will use Vanilla KSP values for its guess")
            return round(period / Vanilla_Dy2s,0)

    def hohmann_transfer_time(self, departure, arrival):
        """Return the ideal half-ellipse transfer time in local KSP days."""
        departure = self._fullname_dic.get(departure, departure)
        arrival = self._fullname_dic.get(arrival, arrival)
        departure_a = float(departure.elements(0)[0])
        arrival_a = float(arrival.elements(0)[0])
        mu = float(departure.mu_central_body)
        if not np.isclose(mu, float(arrival.mu_central_body)):
            raise ValueError("Hohmann estimate requires bodies with the same central attractor")
        transfer_a = 0.5 * (departure_a + arrival_a)
        seconds = np.pi * np.sqrt(transfer_a ** 3 / mu)
        local_day_seconds = 86400.0 / float(self._Edy2Kdy)
        return float(seconds / local_day_seconds)

    def estimate_direct_leg_tof_bounds(
        self,
        seq_fullname,
        profile="relaxed",
        resonance_bounds_by_body=None,
    ):
        """Estimate independent direct-mode ToF bounds for every MGA leg.

        ``planner`` follows KSP-MGA-Planner's 0.7--1.5 Hohmann envelope
        for non-resonant legs and 1--2 body periods for resonant legs.
        ``relaxed`` keeps the same lower bound, permits 2.5 Hohmann times
        on intermediate legs and 2.0 on the terminal leg, and keeps the
        legacy Kerbin 2:1 resonance inside a 1.9--2.1 period envelope.
        """
        profile = str(profile).lower()
        if profile not in ("planner", "relaxed"):
            raise ValueError("Unknown direct ToF estimate profile: " + profile)
        if len(seq_fullname) < 2:
            raise ValueError("A sequence must contain at least two bodies")

        resonance_bounds = {"Kerbin": (1.9, 2.1)} if profile == "relaxed" else {}
        if resonance_bounds_by_body:
            resonance_bounds.update(resonance_bounds_by_body)

        bounds = []
        last_leg_index = len(seq_fullname) - 2
        for leg_index, (departure_name, arrival_name) in enumerate(
            zip(seq_fullname[:-1], seq_fullname[1:])
        ):
            departure = self._fullname_dic[departure_name]
            if departure_name == arrival_name:
                period = float(self.orbital_period(departure))
                lower_factor, upper_factor = resonance_bounds.get(
                    departure_name, (1.0, 2.0)
                )
                lower = period * float(lower_factor)
                upper = period * float(upper_factor)
                kind = "resonant"
                reference = period
            else:
                reference = self.hohmann_transfer_time(departure_name, arrival_name)
                lower_factor = 0.7
                if profile == "planner":
                    upper_factor = 1.5
                elif leg_index == last_leg_index:
                    upper_factor = 2.0
                else:
                    upper_factor = 2.5
                lower = reference * lower_factor
                upper = reference * upper_factor
                kind = "transfer"
            bounds.append(
                {
                    "leg_index": leg_index,
                    "departure": departure_name,
                    "arrival": arrival_name,
                    "kind": kind,
                    "reference_days": reference,
                    "lower_factor": float(lower_factor),
                    "upper_factor": float(upper_factor),
                    "lower_days": float(lower),
                    "upper_days": float(upper),
                }
            )
        return bounds

    def estimate_direct_tof_bounds(self, seq_fullname, profile="relaxed", **kwargs):
        """Return per-leg and summed ToF bounds for direct MGA encoding."""
        legs = self.estimate_direct_leg_tof_bounds(
            seq_fullname, profile=profile, **kwargs
        )
        return {
            "profile": str(profile).lower(),
            "legs": legs,
            "direct_bounds_days": [
                [leg["lower_days"], leg["upper_days"]] for leg in legs
            ],
            "total_lower_days": sum(leg["lower_days"] for leg in legs),
            "total_upper_days": sum(leg["upper_days"] for leg in legs),
        }

    def soi_radius(self, planet):
        if self.planet_pack == "JNSQ":
            return self._planet_pack_module.body[
                self._planet_pack_module.row[planet.name], self._planet_pack_module.col['R_soi (km)']
            ] * 1000
        return self._planet_pack_module.body[
            self._planet_pack_module.row[planet.name], self._planet_pack_module.col['R_soi (km)']
        ] * 1000
        
    def auto_tof(self,seq_fullname,debug = False):
        '''
        Improving the auto_tof feature
        1) when the sequence contains KXXT, with XX a dual assist, plan 1-2Xyr margin.
        2) For the last body, it should never be a full year, half at most, so 0.2-0.4 of period
        '''
        planet_sequence = list(map(self._fullname_dic.get, seq_fullname))
        tof_guess = 0

        for prev_planet, planet in zip(planet_sequence[0:], planet_sequence[1:]):

            if planet == prev_planet:
                tof_guess += 2*self.orbital_period(planet)
            elif planet == planet_sequence[-1] :
                tof_guess += 0.4*self.orbital_period(planet)
            elif planet != prev_planet:  
                tof_guess += self.orbital_period(planet)
            if debug:
                logger.debug("%s %s", planet.name, tof_guess)

        '''I chose tof between 0.55 and 1.1 of the sum'''
        if debug :
            logger.debug(
                "tof_lb: %s tof_ub: %s",
                int(round(tof_guess*0.5,-2)),
                int(round(tof_guess*1.0,-2)),
            )
        return int(round(tof_guess*0.5,-2)), int(round(tof_guess*1.0,-2))
   
    def recalc_results(self):
        raise NotImplementedError(
            "recalc_results was tied to the legacy dataframe datastore. "
            "A SQL recalc path should operate on stored jobs/runs explicitly."
        )
        
    def audit_results(self):
        raise NotImplementedError(
            "audit_results was tied to the legacy dataframe datastore. "
            "A SQL audit should query jobs/results from SQLite directly."
        )
     
    def optimize(
            self,
            n=5000,
            save_it=True,
            sqlite_db_path=None,
            sqlite_batch_name=None,
            sqlite_template=None,
            sqlite_generation_options=None,
            auto_workers=True,
            reserve_cores=2,
            topology=None):
        if sqlite_db_path is None:
            raise ValueError(
                "optimize is SQL-only now: provide sqlite_db_path or call optimize_sqlite(...)."
            )
        return self.optimize_sqlite(
            sqlite_db_path,
            n=n,
            batch_name=sqlite_batch_name,
            auto_workers=auto_workers,
            reserve_cores=reserve_cores,
            topology=topology,
        )
