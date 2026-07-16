"""SQLite orchestration for Tisserand/Lambert scout -> L0 -> funnel."""

from __future__ import annotations

from dataclasses import asdict, dataclass

from _LambertArcFilter import LambertArc1Config
from _SequenceScout import TisserandScoutConfig
from _SQLiteStore import SQLiteJobStore


@dataclass(frozen=True)
class SequenceScoutWorkflowConfig:
    """Serializable production policy for the SQL sequence-scout workflow."""

    t0_min_days: float = 0.0
    t0_max_days: float = 1000.0
    bin_width_days: float = 100.0
    t0_step_days: float = 10.0
    tof_samples: int = 40
    top_sequences: int = 20
    maximum_ejection_ratio: float = 1.05
    l0_seed_base: int = 42000
    l0_islands: int = 64
    l0_population: int = 8
    l0_evo_steps: int = 5
    l0_sade_gen: int = 20
    promotions_per_bin: int = 2
    funnel_islands: int = 16
    funnel_population: int = 32
    funnel_evo_steps: int = 20
    funnel_sade_gen: int = 20
    ejection_altitude_m: float = 100000.0
    injection_altitude_m: float = 100000.0
    arrival_mode: str = "flyby"
    tof_profile: str = "relaxed"

    def __post_init__(self):
        if float(self.t0_max_days) <= float(self.t0_min_days):
            raise ValueError("t0_max_days must exceed t0_min_days")
        if float(self.bin_width_days) <= 0.0:
            raise ValueError("bin_width_days must be positive")
        if float(self.t0_step_days) <= 0.0:
            raise ValueError("t0_step_days must be positive")
        for name in (
            "tof_samples", "top_sequences", "l0_islands", "l0_population",
            "l0_evo_steps", "l0_sade_gen", "promotions_per_bin",
            "funnel_islands", "funnel_population", "funnel_evo_steps",
            "funnel_sade_gen",
        ):
            if int(getattr(self, name)) < 1:
                raise ValueError(name + " must be positive")
        if int(self.l0_population) < 7 or int(self.funnel_population) < 7:
            raise ValueError("PyGMO SADE populations must be at least seven")
        if float(self.maximum_ejection_ratio) <= 0.0:
            raise ValueError("maximum_ejection_ratio must be positive")

    def to_dict(self):
        return asdict(self)


class SequenceScoutSQLWorkflow:
    """Prepare, promote, and inspect a sequence scout in one SQLite database."""

    L0_STRATEGY = "funnel_l0_screen"
    FUNNEL_STRATEGY = "funnel_l0_continuation"

    def __init__(self, wayfinder):
        self.wayfinder = wayfinder

    def _short_name(self, sequence):
        return "".join(
            self.wayfinder._Body_abrev_dic.get(body, body)
            for body in sequence
        )

    def _leg_bounds(self, sequence, lambert_tof_days, profile):
        bounds = [
            list(pair) for pair in self.wayfinder.estimate_direct_tof_bounds(
                sequence, profile=profile,
            )["direct_bounds_days"]
        ]
        hint = float(lambert_tof_days)
        bounds[0][0] = min(float(bounds[0][0]), 0.8 * hint)
        bounds[0][1] = max(float(bounds[0][1]), 1.2 * hint)
        return [[float(low), float(high)] for low, high in bounds]

    def _job_params(
        self, sequence, short_name, bin_start, bin_end, bounds, config,
        optimizer_seed, workflow_stage,
    ):
        arrival_mode, add_vinf_arr, orbit_insertion, e_target = (
            self.wayfinder._resolve_arrival_mode(config.arrival_mode)
        )
        target = self.wayfinder._fullname_dic[sequence[-1]]
        if workflow_stage == "l0":
            opt_level = "sequence_scout_l0"
            optimizer_strategy = self.L0_STRATEGY
            sade_gen = config.l0_sade_gen
            n_island = config.l0_islands
            island_pop = config.l0_population
            n_evo_steps = config.l0_evo_steps
            topology = "unconnected"
        elif workflow_stage == "funnel":
            opt_level = "sequence_scout_funnel"
            optimizer_strategy = self.FUNNEL_STRATEGY
            sade_gen = config.funnel_sade_gen
            n_island = config.funnel_islands
            island_pop = config.funnel_population
            n_evo_steps = config.funnel_evo_steps
            topology = "ring"
        else:
            raise ValueError("Unknown scout workflow stage: " + workflow_stage)
        return {
            "planet_pack": self.wayfinder.planet_pack,
            "sequence": list(sequence),
            "sequence_short_name": short_name,
            "t0_min": float(bin_start),
            "t0_max": float(bin_end),
            "tof_min": float(sum(low for low, _ in bounds)),
            "tof_max": float(sum(high for _, high in bounds)),
            "vinf": [0.8, 1.8],
            "tof_encoding": "direct",
            "tof_profile": str(config.tof_profile) + "_scout_arc1",
            "leg_tof_bounds": bounds,
            "opt_level": opt_level,
            "opt_injection": arrival_mode,
            "arrival_mode": arrival_mode,
            "ejection_altitude": float(config.ejection_altitude_m),
            "injection_altitude": float(config.injection_altitude_m),
            "rp_target": self.wayfinder.rp_target_ward(
                target, config.injection_altitude_m,
            ),
            "e_target": float(e_target),
            "add_vinf_dep": True,
            "add_vinf_arr": bool(add_vinf_arr),
            "orbit_insertion": bool(orbit_insertion),
            "multi_objective": False,
            "lambert_max_revs": 0,
            "sade_gen": int(sade_gen),
            "n_island": int(n_island),
            "island_pop": int(island_pop),
            "n_evo_steps": int(n_evo_steps),
            "optimizer_topology": topology,
            "optimizer_seed": int(optimizer_seed),
            "optimizer_strategy": optimizer_strategy,
        }

    def prepare(
        self, db_path, name, start_body, target_body, config=None,
        candidate_bodies=None, scout_config=None, lambert_config=None,
    ):
        """Persist a deterministic scout and create all selected L0 jobs."""
        config = config or SequenceScoutWorkflowConfig()
        if isinstance(config, dict):
            config = SequenceScoutWorkflowConfig(**config)
        scout_config = scout_config or TisserandScoutConfig()
        if isinstance(scout_config, dict):
            scout_config = TisserandScoutConfig(**scout_config)
        lambert_config = lambert_config or LambertArc1Config(
            parking_altitude_m=config.ejection_altitude_m,
        )
        if isinstance(lambert_config, dict):
            lambert_config = LambertArc1Config(**lambert_config)

        candidates = self.wayfinder.scout_sequences(
            start_body, target_body, candidate_bodies=candidate_bodies,
            config=scout_config,
        )
        scan = self.wayfinder.scan_scout_sequence_bins(
            candidates,
            [config.t0_min_days, config.t0_max_days],
            config=lambert_config,
            bin_width_days=config.bin_width_days,
            t0_step_days=config.t0_step_days,
            tof_samples=config.tof_samples,
            maximum_ejection_ratio=config.maximum_ejection_ratio,
        )
        ranked = scan.ranked_unique(config.top_sequences)
        rank_by_sequence = {
            row.candidate.sequence: index + 1
            for index, row in enumerate(ranked)
        }
        selected_sequences = set(rank_by_sequence)
        selected = sorted(
            (
                row for row in scan.candidates
                if row.candidate.sequence in selected_sequences
            ),
            key=lambda row: (
                row.bin_start_days, row.score, row.candidate.sequence,
            ),
        )
        persisted_config = {
            "workflow": config.to_dict(),
            "tisserand": scout_config.to_dict(),
            "lambert": lambert_config.to_dict(),
            "candidate_bodies": (
                list(candidate_bodies) if candidate_bodies is not None else None
            ),
        }
        store = SQLiteJobStore(db_path)
        try:
            scout_run_id = store.upsert_sequence_scout_run(
                name, self.wayfinder.planet_pack, start_body, target_body,
                persisted_config,
                direct_reference=scan.direct_reference.to_dict(),
                status="L0_PENDING",
                raw_candidate_count=len(candidates),
                viable_candidate_count=len(scan.candidates),
            )
            batch_name = str(name) + "__l0"
            batch_id = store.upsert_batch(
                batch_name, self.wayfinder.planet_pack,
                template={"start": start_body, "target": target_body},
                generation_options={
                    "scout_run_id": scout_run_id,
                    "workflow_stage": "l0",
                    "optimizer_strategy": self.L0_STRATEGY,
                    "config": config.to_dict(),
                },
                purpose="sequence_scout_l0",
            )
            job_ids = []
            for index, row in enumerate(selected):
                sequence = tuple(row.candidate.sequence)
                short_name = self._short_name(sequence)
                bounds = self._leg_bounds(
                    sequence, row.solution.tof_days, config.tof_profile,
                )
                sequence_id = store.upsert_sequence(
                    self.wayfinder.planet_pack, short_name, sequence,
                )
                candidate_id = store.upsert_sequence_scout_candidate(
                    scout_run_id, sequence_id,
                    row.bin_start_days, row.bin_end_days,
                    rank_by_sequence[sequence], row.score, row.ejection_ratio,
                    row.solution.t0_days, row.solution.tof_days, bounds,
                    selected_for_l0=True,
                )
                seed = int(config.l0_seed_base) + index
                job_id = store.upsert_job(
                    self._job_params(
                        sequence, short_name, row.bin_start_days,
                        row.bin_end_days, bounds, config, seed, "l0",
                    ),
                    batch_id=batch_id,
                    status="TODO",
                )
                store.link_sequence_scout_job(
                    candidate_id, job_id, "l0",
                    metadata={"global_sequence_rank": rank_by_sequence[sequence]},
                )
                job_ids.append(job_id)
            return {
                "scout_run_id": scout_run_id,
                "name": str(name),
                "l0_batch_name": batch_name,
                "l0_batch_id": int(batch_id),
                "raw_candidate_count": len(candidates),
                "viable_candidate_count": len(scan.candidates),
                "selected_sequence_count": len(selected_sequences),
                "l0_job_count": len(job_ids),
                "direct_reference": scan.direct_reference.to_dict(),
            }
        finally:
            store.close()

    def promote(
        self, db_path, scout_run_id, per_bin=None, allow_partial=False,
    ):
        """Promote the best completed L0 jobs in every departure bin."""
        store = SQLiteJobStore(db_path)
        try:
            scout_run = store.sequence_scout_run(scout_run_id=scout_run_id)
            if scout_run is None:
                raise ValueError("Unknown scout_run_id: " + str(scout_run_id))
            status = store.sequence_scout_status(scout_run_id)
            incomplete_l0 = sum(
                int(status["job_counts"].get("l0:" + job_status, 0))
                for job_status in ("TODO", "RUNNING")
            )
            if incomplete_l0 and not allow_partial:
                raise RuntimeError(
                    "Cannot promote while {} L0 jobs are still pending or "
                    "running; pass allow_partial=True only for diagnostics"
                    .format(incomplete_l0)
                )
            config = SequenceScoutWorkflowConfig(
                **scout_run["config"]["workflow"]
            )
            per_bin = (
                int(config.promotions_per_bin)
                if per_bin is None else int(per_bin)
            )
            if per_bin < 1:
                raise ValueError("per_bin must be positive")
            results = store.sequence_scout_l0_results(scout_run_id)
            by_bin = {}
            for row in results:
                by_bin.setdefault(float(row["bin_start"]), []).append(row)
            selected = []
            for bin_start in sorted(by_bin):
                selected.extend(sorted(
                    by_bin[bin_start],
                    key=lambda row: (
                        float(row["objective_dv"]), row["sequence_short_name"],
                    ),
                )[:per_bin])
            batch_name = str(scout_run["name"]) + "__funnel"
            batch_id = store.upsert_batch(
                batch_name, self.wayfinder.planet_pack,
                template={
                    "start": scout_run["start_body"],
                    "target": scout_run["target_body"],
                },
                generation_options={
                    "scout_run_id": int(scout_run_id),
                    "workflow_stage": "funnel",
                    "optimizer_strategy": self.FUNNEL_STRATEGY,
                    "promotions_per_bin": per_bin,
                    "config": config.to_dict(),
                },
                purpose="sequence_scout_funnel",
            )
            existing_funnel_jobs = {
                int(row["job_id"])
                for row in store.conn.execute(
                    """
                    SELECT sj.job_id
                    FROM sequence_scout_jobs sj
                    JOIN sequence_scout_candidates c
                      ON c.id = sj.scout_candidate_id
                    WHERE c.scout_run_id = ?
                      AND sj.workflow_stage = 'funnel'
                    """,
                    (int(scout_run_id),),
                )
            }
            promoted = []
            newly_promoted = []
            rank_by_bin = {}
            for row in selected:
                bin_start = float(row["bin_start"])
                rank_by_bin[bin_start] = rank_by_bin.get(bin_start, 0) + 1
                sequence = tuple(row["bodies"])
                seed = int(row["optimizer_seed"])
                job_id = store.upsert_job(
                    self._job_params(
                        sequence, row["sequence_short_name"], row["bin_start"],
                        row["bin_end"], row["leg_tof_bounds"], config, seed,
                        "funnel",
                    ),
                    batch_id=batch_id,
                    status="TODO",
                )
                promotion_rank = rank_by_bin[bin_start]
                store.link_sequence_scout_job(
                    row["scout_candidate_id"], job_id, "funnel",
                    parent_job_id=row["job_id"], parent_run_id=row["run_id"],
                    promotion_rank=promotion_rank,
                    metadata={
                        "l0_objective_dv": float(row["objective_dv"]),
                        "promotion_policy": "best_per_t0_bin",
                    },
                )
                promoted.append(job_id)
                if job_id not in existing_funnel_jobs:
                    newly_promoted.append(job_id)
                    existing_funnel_jobs.add(job_id)
            store.upsert_sequence_scout_run(
                scout_run["name"], scout_run["planet_pack"],
                scout_run["start_body"], scout_run["target_body"],
                scout_run["config"],
                direct_reference=scout_run["direct_reference"],
                status="FUNNEL_PENDING" if promoted else "L0_INCOMPLETE",
                raw_candidate_count=scout_run["raw_candidate_count"],
                viable_candidate_count=scout_run["viable_candidate_count"],
            )
            return {
                "scout_run_id": int(scout_run_id),
                "funnel_batch_name": batch_name,
                "funnel_batch_id": int(batch_id),
                "completed_l0_count": len(results),
                "promoted_job_count": len(promoted),
                "newly_promoted_job_count": len(newly_promoted),
                "promotions_per_bin": per_bin,
                "partial_promotion": bool(allow_partial),
            }
        finally:
            store.close()

    def status(self, db_path, scout_run_id):
        store = SQLiteJobStore(db_path)
        try:
            return store.sequence_scout_status(scout_run_id)
        finally:
            store.close()
