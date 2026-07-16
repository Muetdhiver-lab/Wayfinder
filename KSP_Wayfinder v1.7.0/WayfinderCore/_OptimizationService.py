"""Technical helpers for optimizer orchestration.

This is deliberately a narrow extraction from ``Wayfinder``. Public entry
points remain on the facade while optimizer policy moves into testable helpers.
"""

import json
import os
from dataclasses import dataclass, field

import numpy as np
import pygmo as pg

from _Optimization import alpha_leg_tofs, direct_leg_tofs


@dataclass
class StageConfig:
    """Serializable description of one optimizer funnel stage."""

    ALGORITHMS = {"hybrid", "sade", "sade_nlopt", "mbh"}
    EJECTION_MODELS = {"approximate", "vector_3d"}
    TOPOLOGIES = {
        "unconnected", "ring", "fully_connected", "fully-connected", "full",
        "split_ring_16_4",
    }

    name: str
    n_island: int
    island_pop: int
    evo_steps: int
    sade_gen: int
    ejection_model: str
    initialization: str
    algorithm: str
    options: dict = field(default_factory=dict)

    def __post_init__(self):
        self.name = str(self.name)
        self.n_island = self._positive_int("n_island", self.n_island)
        self.island_pop = self._positive_int("island_pop", self.island_pop)
        if self.island_pop < 7:
            raise ValueError("island_pop must be at least 7 for Pygmo SADE")
        self.evo_steps = self._positive_int("evo_steps", self.evo_steps)
        self.sade_gen = self._positive_int("sade_gen", self.sade_gen)
        self.ejection_model = self._choice(
            "ejection_model", self.ejection_model, self.EJECTION_MODELS,
        )
        self.initialization = str(self.initialization)
        self.algorithm = self._choice("algorithm", self.algorithm, self.ALGORITHMS)
        self.options = dict(self.options or {})
        self._validate_options()

    @staticmethod
    def _positive_int(name, value):
        value = int(value)
        if value < 1:
            raise ValueError(name + " must be positive")
        return value

    @staticmethod
    def _non_negative_int(name, value):
        value = int(value)
        if value < 0:
            raise ValueError(name + " must be non-negative")
        return value

    @staticmethod
    def _choice(name, value, allowed):
        value = str(value)
        if value not in allowed:
            raise ValueError(
                "{} must be one of: {}".format(name, ", ".join(sorted(allowed)))
            )
        return value

    def _validate_options(self):
        topology = self.options.get("topology")
        if topology is not None:
            topology = self._choice(
                "topology", self.options["topology"], self.TOPOLOGIES,
            )
            self.options["topology"] = topology
        for name in (
            "migration_rate", "archive_interval", "archive_size",
        ):
            if name in self.options:
                self.options[name] = self._non_negative_int(
                    name, self.options[name],
                )
        if "adaptive_stop" in self.options and not isinstance(
            self.options["adaptive_stop"], dict
        ):
            raise ValueError("adaptive_stop must be a dict")
        if "annealing" in self.options and not isinstance(
            self.options["annealing"], dict
        ):
            raise ValueError("annealing must be a dict")
        has_split_ring = "split_ring" in self.options
        if topology == "split_ring_16_4" and not has_split_ring:
            raise ValueError(
                "split_ring_16_4 topology requires split_ring options"
            )
        if has_split_ring and topology != "split_ring_16_4":
            raise ValueError(
                "split_ring options require split_ring_16_4 topology"
            )
        if has_split_ring:
            split_ring = self.options["split_ring"]
            if not isinstance(split_ring, dict):
                raise ValueError("split_ring must be a dict")
            split_ring = dict(split_ring)
            for name in ("current_islands", "alternative_islands"):
                if name not in split_ring:
                    raise ValueError("split_ring is missing " + name)
                split_ring[name] = self._positive_int(name, split_ring[name])
                if split_ring[name] < 3:
                    raise ValueError(
                        "split_ring partitions require at least three islands"
                    )
            bridge_weight = float(split_ring.get("bridge_weight", 0.25))
            if not 0.0 <= bridge_weight <= 1.0:
                raise ValueError("split_ring bridge_weight must be in [0, 1]")
            split_ring["bridge_weight"] = bridge_weight
            if (
                split_ring["current_islands"]
                + split_ring["alternative_islands"] != self.n_island
            ):
                raise ValueError(
                    "split_ring partitions must match n_island"
                )
            self.options["split_ring"] = split_ring

    @classmethod
    def from_dict(cls, values):
        required = {
            "name", "n_island", "island_pop", "evo_steps", "sade_gen",
            "ejection_model", "initialization", "algorithm",
        }
        missing = sorted(required - set(values))
        if missing:
            raise ValueError(
                "Missing stage configuration fields: " + ", ".join(missing)
            )
        options = {
            key: value for key, value in dict(values).items()
            if key not in required
        }
        return cls(
            name=str(values["name"]),
            n_island=values["n_island"],
            island_pop=values["island_pop"],
            evo_steps=values["evo_steps"],
            sade_gen=values["sade_gen"],
            ejection_model=str(values["ejection_model"]),
            initialization=str(values["initialization"]),
            algorithm=str(values["algorithm"]),
            options=options,
        )

    def to_dict(self):
        values = {
            "name": self.name,
            "n_island": int(self.n_island),
            "island_pop": int(self.island_pop),
            "evo_steps": int(self.evo_steps),
            "sade_gen": int(self.sade_gen),
            "ejection_model": self.ejection_model,
            "initialization": self.initialization,
            "algorithm": self.algorithm,
        }
        values.update(self.options)
        return values


@dataclass
class FunnelConfig:
    """Serializable optimizer funnel preset."""

    name: str
    exact_strategy: str
    requested: dict
    stages: list
    kind: str = "funnel"
    pressure_cascade: dict = None

    @classmethod
    def from_stage_dicts(
        cls, name, exact_strategy, requested, stage_dicts, kind="funnel",
        pressure_cascade=None,
    ):
        return cls(
            name=str(name),
            exact_strategy=str(exact_strategy),
            requested=dict(requested),
            stages=[StageConfig.from_dict(stage) for stage in stage_dicts],
            kind=str(kind),
            pressure_cascade=(
                dict(pressure_cascade) if pressure_cascade is not None else None
            ),
        )

    def stage_dicts(self):
        return [stage.to_dict() for stage in self.stages]

    def to_dict(self):
        values = {
            "kind": self.kind,
            "optimizer_strategy": self.name,
            "exact_strategy": self.exact_strategy,
            "requested": dict(self.requested),
            "stages": self.stage_dicts(),
        }
        if self.pressure_cascade is not None:
            values["pressure_cascade"] = dict(self.pressure_cascade)
        return values


class OptimizationService:
    """Build topologies and report archipelago behavior."""

    FUNNEL_STRATEGIES = {
        "funnel": "legacy",
        # SQL sequence-scout workflow. The screen is a standalone L0; the
        # continuation consumes its persisted population and starts at L1.
        "funnel_l0_screen": "l0_screen",
        "funnel_l0_continuation": "scout_archive_nm_64_mbh_between",
        "funnel_local_exact": "local",
        "funnel_hybrid_exact": "hybrid",
        "funnel_phase_elites_nm": "phase_elites_nm",
        "funnel_phase_elites_equal": "phase_elites_nm_equal",
        # Public L0-oriented names. These keep the historical
        # scout/archive implementation but make the intent explicit:
        # run a wide, shallow, unconnected level-0 before the normal funnel.
        # The unsuffixed policy keeps the production compromise at 64 islands.
        "funnel_l0_recall": "scout_archive_nm_64",
        "funnel_l0_recall_32": "scout_archive_nm_32",
        "funnel_l0_recall_64": "scout_archive_nm_64",
        "funnel_l0_recall_128": "scout_archive_nm_128",
        "funnel_scout_archive": "scout_archive_nm",
        "funnel_scout_archive_32": "scout_archive_nm_32",
        "funnel_scout_archive_64": "scout_archive_nm_64",
        "funnel_scout_archive_128": "scout_archive_nm_128",
    }

    @classmethod
    def funnel_exact_strategy(cls, optimizer_strategy):
        optimizer_strategy = str(optimizer_strategy)
        if optimizer_strategy.endswith("_pressure_cascade"):
            optimizer_strategy = optimizer_strategy[:-len("_pressure_cascade")]
        exact_strategy = cls.FUNNEL_STRATEGIES.get(optimizer_strategy)
        if exact_strategy is not None:
            return exact_strategy
        prefix = "funnel_l0_recall_"
        if str(optimizer_strategy).startswith(prefix):
            parts = str(optimizer_strategy)[len(prefix):].split("_")
            if parts and parts[0] in ("32", "64", "128"):
                exact_parts = ["scout", "archive", "nm", parts[0]]
                remaining = parts[1:]
                if remaining and remaining[0].startswith("g"):
                    generation_part = remaining.pop(0)
                    if (
                        not generation_part[1:].isdigit()
                        or int(generation_part[1:]) <= 0
                    ):
                        return None
                    exact_parts.append(generation_part)
                if remaining and remaining[0].startswith("w"):
                    wide_part = remaining.pop(0)
                    if (
                        not wide_part[1:].isdigit()
                        or int(wide_part[1:]) <= 0
                    ):
                        return None
                    exact_parts.append(wide_part)
                if remaining and remaining[0] == "mbh":
                    if len(remaining) < 2 or remaining[1] not in ("between", "l2"):
                        return None
                    exact_parts.extend(["mbh", remaining[1]])
                    remaining = remaining[2:]
                if remaining[:2] == ["pareto", "l0"]:
                    exact_parts.extend(["pareto", "l0"])
                    remaining = remaining[2:]
                elif remaining and remaining[0] == "pareto":
                    exact_parts.append("pareto")
                    remaining = remaining[1:]
                elif remaining[:2] == ["basin", "l0"]:
                    exact_parts.extend(["basin", "l0"])
                    remaining = remaining[2:]
                elif remaining[:4] == ["hill", "valley", "p2", "l0"]:
                    exact_parts.extend(["hill", "valley", "p2", "l0"])
                    remaining = remaining[4:]
                elif remaining[:4] == ["hill", "valley", "mr", "l0"]:
                    exact_parts.extend(["hill", "valley", "mr", "l0"])
                    remaining = remaining[4:]
                elif remaining[:4] == ["hill", "valley", "mr32", "l0"]:
                    exact_parts.extend(["hill", "valley", "mr32", "l0"])
                    remaining = remaining[4:]
                elif remaining[:3] == ["hill", "valley", "l0"]:
                    exact_parts.extend(["hill", "valley", "l0"])
                    remaining = remaining[3:]
                elif remaining[:4] == ["portfolio", "16", "4", "l0"]:
                    exact_parts.extend(["portfolio", "16", "4", "l0"])
                    remaining = remaining[4:]
                if remaining == ["exact"]:
                    exact_parts.append("exact")
                elif remaining:
                    return None
                return "_".join(exact_parts)
        return None

    @staticmethod
    def funnel_pressure_cascade_config(optimizer_strategy):
        """Return production L0 pressure-cascade settings for a strategy."""
        if not str(optimizer_strategy).endswith("_pressure_cascade"):
            return None
        return {
            "enabled": True,
            "source_stage_index": 1,
            "top": 32,
            "near_fraction": 0.03,
            "min_pressure_count": 2,
            "relax_fraction": 0.20,
            "min_relax_days": 5.0,
            "max_relax_days": 60.0,
            "min_leg_lower_days": 1.0,
            "adjustment_mode": "widen",
            "branch_policy": "l1_combined",
            "max_l1_improvement": 0.15,
            "min_l1_best_to_l0_median": 0.725,
            "rescue_mode": "cascade",
            "cascade_min_improvement": 0.10,
            "retry_seed_offset": 100,
        }

    @staticmethod
    def automatic_worker_count(reserve_cores=2):
        """Return a conservative worker count based on the local CPU count."""
        cpu_count = os.cpu_count() or 1
        return max(1, int(cpu_count) - int(reserve_cores))

    @classmethod
    def resolve_island_count(
        cls, requested_islands, auto_workers=True, reserve_cores=2,
    ):
        requested_islands = int(requested_islands)
        if not auto_workers:
            return requested_islands
        return min(
            requested_islands,
            cls.automatic_worker_count(reserve_cores=reserve_cores),
        )

    @staticmethod
    def make_archipelago_topology(name, n_islands):
        name = str(name).lower()
        # Pygmo adds one topology vertex per pushed island. Pre-sizing here
        # would create twice as many vertices and invalid migration edges.
        if int(n_islands) < 1:
            raise ValueError("n_islands must be positive")
        if name == "unconnected":
            return pg.unconnected()
        if name == "ring":
            return pg.ring()
        if name in ("fully_connected", "fully-connected", "full"):
            return pg.fully_connected()
        raise ValueError("Unknown archipelago topology: " + str(name))

    @staticmethod
    def make_split_ring_topology(
        current_islands=16, alternative_islands=4, bridge_weight=0.25,
    ):
        """Return two bidirectional rings joined by one weak bridge."""
        current_islands = int(current_islands)
        alternative_islands = int(alternative_islands)
        bridge_weight = float(bridge_weight)
        if current_islands < 3 or alternative_islands < 3:
            raise ValueError(
                "Each split-ring partition requires at least three islands"
            )
        if not 0.0 <= bridge_weight <= 1.0:
            raise ValueError("bridge_weight must be between 0 and 1")
        topology = pg.free_form()
        total = current_islands + alternative_islands
        for _ in range(total):
            topology.add_vertex()

        def add_ring(start, count):
            for offset in range(count):
                left = start + offset
                right = start + ((offset + 1) % count)
                topology.add_edge(left, right, 1.0)
                topology.add_edge(right, left, 1.0)

        add_ring(0, current_islands)
        add_ring(current_islands, alternative_islands)
        topology.add_edge(
            current_islands - 1, current_islands, bridge_weight,
        )
        topology.add_edge(
            current_islands, current_islands - 1, bridge_weight,
        )
        return topology

    @staticmethod
    def population_id_origins(populations):
        """Map persistent Pygmo IDs to every island containing their copies."""
        origins = {}
        for island_index, population in enumerate(populations):
            for identifier in population.get_ID():
                origins.setdefault(int(identifier), set()).add(island_index)
        return origins

    @staticmethod
    def archipelago_telemetry(
        archipelago, populations, previous_id_origins, topology_name,
    ):
        """Measure topology shape and observed migration for one evolution."""
        island_count = len(populations)
        normalized_topology = str(topology_name).lower()
        if normalized_topology == "unconnected":
            topology_vertices = island_count
            topology_edges = 0
        else:
            graph = archipelago.get_topology().to_networkx()
            topology_vertices = int(graph.number_of_nodes())
            topology_edges = int(len(graph.edges()))
        if topology_vertices != island_count:
            raise RuntimeError(
                "Archipelago topology has {} vertices for {} islands".format(
                    topology_vertices, island_count
                )
            )

        migrants_db = archipelago.get_migrants_db()
        migrants_published = sum(len(entry[0]) for entry in migrants_db)
        migration_islands_active = sum(bool(len(entry[0])) for entry in migrants_db)
        migrations_accepted = sum(
            1
            for island_index, population in enumerate(populations)
            for identifier in population.get_ID()
            if int(identifier) in previous_id_origins
            and island_index not in previous_id_origins[int(identifier)]
        )
        return {
            "island_count": island_count,
            "topology_vertices": topology_vertices,
            "topology_edges": topology_edges,
            "migrants_published": int(migrants_published),
            "migration_islands_active": int(migration_islands_active),
            # Observable lower bound: a multi-hop migrant is counted only at
            # its final island after one asynchronous evolve call.
            "migrations_accepted": int(migrations_accepted),
        }

    @staticmethod
    def balanced_annealing_schedule(problem, population_size, sade_gen):
        """Match simulated-annealing evaluations to one SADE evolve call."""
        dimension = int(problem.get_nx())
        target = max(1, int(population_size) * int(sade_gen))
        n_temperature = 2 if target >= 2 * dimension else 1
        n_range = 2 if target >= 2 * n_temperature * dimension else 1
        bin_size = max(
            1, int(round(target / (n_temperature * n_range * dimension)))
        )
        return {
            "n_T_adj": n_temperature,
            "n_range_adj": n_range,
            "bin_size": bin_size,
            "target_evaluations": target,
            "nominal_evaluations": (
                n_temperature * n_range * bin_size * dimension
            ),
        }

    @staticmethod
    def funnel_stage_plan(
        n_islands, island_pop, evo_steps, sade_gen, exact_strategy="legacy",
    ):
        """Return a strictly narrowing optimizer plan."""
        if exact_strategy == "l0_screen":
            return [{
                "name": "scout_unconnected",
                "n_island": max(1, int(n_islands)),
                "island_pop": max(7, int(island_pop)),
                "evo_steps": max(1, int(evo_steps)),
                "sade_gen": max(1, int(sade_gen)),
                "ejection_model": "approximate",
                "initialization": "random",
                "algorithm": "sade",
                "topology": "unconnected",
                "migration_rate": 0,
            }]
        n_stage_1 = max(1, int(n_islands))
        n_stage_2 = max(1, n_stage_1 // 2) if n_stage_1 > 1 else 1
        n_stage_3 = min(8, max(1, n_stage_2 // 2)) if n_stage_2 > 1 else 1
        island_pop = max(7, int(island_pop))
        focused_pop = max(7, min(island_pop, 14))
        stages = [
            {
                "name": "wide",
                "n_island": n_stage_1,
                "island_pop": island_pop,
                "evo_steps": int(evo_steps),
                "sade_gen": int(sade_gen),
                "ejection_model": "approximate",
                "initialization": "random_plus_champion",
                "algorithm": "hybrid",
                "annealing": {"Ts": 3000.0, "Tf": 10.0},
            },
            {
                "name": "intermediate",
                "n_island": n_stage_2,
                "island_pop": focused_pop,
                "evo_steps": int(evo_steps) * 2,
                "sade_gen": int(sade_gen),
                "ejection_model": "approximate",
                "initialization": "random_plus_champion",
                "algorithm": "hybrid",
                "annealing": {"Ts": 1000.0, "Tf": 1.0},
            },
            {
                "name": "exact_ejection",
                "n_island": n_stage_3,
                "island_pop": focused_pop,
                "evo_steps": 5,
                "sade_gen": 500,
                "ejection_model": "vector_3d",
                "initialization": "random_plus_champion",
                "sigma_fraction": 0.04,
                "algorithm": "hybrid",
                "annealing": {"Ts": 200.0, "Tf": 0.1},
            },
        ]
        if exact_strategy in ("local", "hybrid"):
            stages[-1].update({
                "evo_steps": 10,
                "sade_gen": 50,
                "initialization": (
                    "local" if exact_strategy == "local" else "local_global"
                ),
                "algorithm": "sade",
                "adaptive_stop": {
                    "min_steps": 2,
                    "window": 2,
                    "patience": 2,
                    "best_relative_tolerance": 0.002,
                    "average_relative_tolerance": 0.01,
                },
            })
        elif exact_strategy != "legacy":
            if exact_strategy in ("phase_elites_nm", "phase_elites_nm_equal"):
                stages[1].update({
                    "initialization": "phase_elites_mixed",
                    "elite_fraction": 0.35,
                })
                if exact_strategy == "phase_elites_nm_equal":
                    stages[1]["evo_steps"] = stages[0]["evo_steps"]
                stages[2].update({
                    "n_island": min(8, n_stage_2),
                    "evo_steps": 5,
                    "sade_gen": 100,
                    "algorithm": "sade_nlopt",
                })
            elif exact_strategy.startswith("scout_archive_nm"):
                scout_island_options = {
                    "scout_archive_nm_32": 32,
                    "scout_archive_nm": 64,
                    "scout_archive_nm_64": 64,
                    "scout_archive_nm_128": 128,
                }
                scout_generations_override = None
                wide_islands_override = None
                exact_l0_selection = False
                mbh_mode = None
                handoff_policy = None
                if exact_strategy not in scout_island_options:
                    parts = str(exact_strategy).split("_")
                    if (
                        5 <= len(parts) <= 13
                        and parts[:3] == ["scout", "archive", "nm"]
                        and parts[3] in ("32", "64", "128")
                    ):
                        scout_island_options[exact_strategy] = int(parts[3])
                        suffixes = parts[4:]
                        if suffixes and suffixes[0].startswith("g"):
                            if (
                                not suffixes[0][1:].isdigit()
                                or int(suffixes[0][1:]) <= 0
                            ):
                                raise ValueError(
                                    "Unknown scout/archive strategy: "
                                    + str(exact_strategy)
                            )
                            scout_generations_override = int(suffixes.pop(0)[1:])
                        if suffixes and suffixes[0].startswith("w"):
                            if (
                                not suffixes[0][1:].isdigit()
                                or int(suffixes[0][1:]) <= 0
                            ):
                                raise ValueError(
                                    "Unknown scout/archive strategy: "
                                    + str(exact_strategy)
                                )
                            wide_islands_override = int(suffixes.pop(0)[1:])
                        if suffixes and suffixes[0] == "mbh":
                            if (
                                len(suffixes) < 2
                                or suffixes[1] not in ("between", "l2")
                            ):
                                raise ValueError(
                                    "Unknown scout/archive strategy: "
                                    + str(exact_strategy)
                                )
                            mbh_mode = suffixes[1]
                            suffixes = suffixes[2:]
                        if suffixes[:2] == ["pareto", "l0"]:
                            handoff_policy = "pareto_l0"
                            suffixes = suffixes[2:]
                        elif suffixes and suffixes[0] == "pareto":
                            handoff_policy = "pareto_all"
                            suffixes = suffixes[1:]
                        elif suffixes[:2] == ["basin", "l0"]:
                            handoff_policy = "basin_l0"
                            suffixes = suffixes[2:]
                        elif suffixes[:4] == ["hill", "valley", "p2", "l0"]:
                            handoff_policy = "hill_valley_p2_l0"
                            suffixes = suffixes[4:]
                        elif suffixes[:4] == ["hill", "valley", "mr", "l0"]:
                            handoff_policy = "hill_valley_mr_l0"
                            suffixes = suffixes[4:]
                        elif suffixes[:4] == ["hill", "valley", "mr32", "l0"]:
                            handoff_policy = "hill_valley_mr32_l0"
                            suffixes = suffixes[4:]
                        elif suffixes[:3] == ["hill", "valley", "l0"]:
                            handoff_policy = "hill_valley_l0"
                            suffixes = suffixes[3:]
                        elif suffixes[:4] == ["portfolio", "16", "4", "l0"]:
                            handoff_policy = "portfolio_16_4_l0"
                            suffixes = suffixes[4:]
                        if suffixes == ["exact"]:
                            exact_l0_selection = True
                        elif suffixes:
                            raise ValueError(
                                "Unknown scout/archive strategy: "
                                + str(exact_strategy)
                            )
                    else:
                        raise ValueError(
                            "Unknown scout/archive strategy: " + str(exact_strategy)
                        )
                if exact_strategy not in scout_island_options:
                    raise ValueError(
                        "Unknown scout/archive strategy: " + str(exact_strategy)
                    )
                scout_islands = scout_island_options[exact_strategy]
                scout_population = 8
                scout_steps = 5
                # Keep the L0 evaluation budget constant while varying width.
                scout_generations = max(
                    1, int(round(int(sade_gen) * 64 / scout_islands))
                )
                if scout_generations_override is not None:
                    scout_generations = scout_generations_override
                if wide_islands_override is not None:
                    wide_islands = int(wide_islands_override)
                    stages[0]["n_island"] = wide_islands
                    stages[1]["n_island"] = (
                        max(1, wide_islands // 2) if wide_islands > 1 else 1
                    )
                    stages[2]["n_island"] = min(
                        8,
                        max(1, stages[1]["n_island"] // 2)
                        if stages[1]["n_island"] > 1 else 1,
                    )
                wide_step_cost = max(
                    1, stages[0]["n_island"] * island_pop * int(sade_gen)
                )
                scout_step_cost = (
                    scout_islands * scout_population * scout_generations
                )
                wide_steps_replaced = int(np.ceil(
                    scout_step_cost * scout_steps / wide_step_cost
                ))
                stages[0].update({
                    "evo_steps": max(5, int(evo_steps) - wide_steps_replaced),
                    "initialization": "scout_diverse",
                    "archive_exact": True,
                    "archive_interval": 5,
                    "archive_size": 32,
                })
                if exact_l0_selection:
                    stages[0]["selection_problem"] = "exact"
                stages[1].update({
                    "initialization": "phase_elites_mixed",
                    "elite_fraction": 0.35,
                    "archive_exact": True,
                    "archive_interval": 5,
                    "archive_size": 32,
                })
                stages[2].update({
                    "n_island": min(8, n_stage_2),
                    "evo_steps": 10,
                    "sade_gen": 100,
                    "algorithm": "sade_nlopt",
                    "use_exact_archive": True,
                    "adaptive_stop": {
                        "min_steps": 5,
                        "window": 2,
                        "patience": 2,
                        "best_relative_tolerance": 0.001,
                        "average_relative_tolerance": 0.005,
                        "require_average_plateau": False,
                    },
                })
                mbh_stage = {
                    "name": "mbh_refine",
                    "n_island": stages[1]["n_island"],
                    "island_pop": focused_pop,
                    "evo_steps": 1,
                    "sade_gen": int(sade_gen),
                    "ejection_model": "approximate",
                    "initialization": "phase_elites_mixed",
                    "elite_fraction": 0.35,
                    "algorithm": "mbh",
                    "topology": "unconnected",
                    "migration_rate": 0,
                    "archive_exact": True,
                    "archive_interval": 1,
                    "archive_size": 32,
                    "mbh_stop": 3,
                    "mbh_perturb": 0.05,
                }
                if mbh_mode == "between":
                    stages.insert(1, dict(mbh_stage))
                elif mbh_mode == "l2":
                    mbh_l2 = dict(mbh_stage)
                    mbh_l2["name"] = "intermediate_mbh"
                    mbh_l2["evo_steps"] = max(1, stages[1]["evo_steps"] // 5)
                    stages[1].update(mbh_l2)
                stages.insert(0, {
                    "name": "scout_unconnected",
                    "n_island": scout_islands,
                    "island_pop": scout_population,
                    "evo_steps": scout_steps,
                    "sade_gen": scout_generations,
                    "ejection_model": "approximate",
                    "initialization": "random",
                    "algorithm": "sade",
                    "topology": "unconnected",
                    "migration_rate": 0,
                })
                if handoff_policy is not None:
                    for stage in stages:
                        use_policy = stage["initialization"] == "scout_diverse"
                        use_policy = use_policy or (
                            handoff_policy == "pareto_all"
                            and stage["initialization"] == "phase_elites_mixed"
                        )
                        if use_policy:
                            stage["selection_policy"] = handoff_policy
                            if handoff_policy in (
                                "hill_valley_l0", "hill_valley_p2_l0",
                                "hill_valley_mr_l0", "hill_valley_mr32_l0",
                            ):
                                stage["elite_fraction"] = 0.35
                            if handoff_policy == "hill_valley_p2_l0":
                                stage["barrier_relative_tolerance"] = 2.0
                            if handoff_policy in (
                                "hill_valley_mr_l0", "hill_valley_mr32_l0",
                            ):
                                stage["valley_slot_fraction"] = 0.75
                            if handoff_policy == "hill_valley_mr32_l0":
                                stage["n_island"] = 32
                                stage["island_pop"] = 16
                            if handoff_policy == "portfolio_16_4_l0":
                                stage.update({
                                    "n_island": 20,
                                    "initialization": "preselected",
                                    "topology": "split_ring_16_4",
                                    "split_ring": {
                                        "current_islands": 16,
                                        "alternative_islands": 4,
                                        "bridge_weight": 0.25,
                                    },
                                })
            else:
                raise ValueError(
                    "Unknown exact funnel strategy: " + str(exact_strategy)
                )
        return stages

    @classmethod
    def funnel_run_config(
        cls, optimizer_strategy, n_islands, island_pop, evo_steps, sade_gen,
    ):
        """Return the complete canonical funnel configuration for one run."""
        exact_strategy = cls.funnel_exact_strategy(optimizer_strategy)
        if exact_strategy is None:
            return None
        stages = cls.funnel_stage_plan(
            n_islands,
            island_pop,
            evo_steps,
            sade_gen,
            exact_strategy=exact_strategy,
        )
        if str(optimizer_strategy) == "funnel_l0_continuation":
            if not stages or stages[0]["name"] != "scout_unconnected":
                raise ValueError(
                    "L0 continuation requires a scout_unconnected first stage"
                )
            stages = stages[1:]
        funnel = FunnelConfig.from_stage_dicts(
            str(optimizer_strategy),
            exact_strategy,
            {
                "n_islands": int(n_islands),
                "island_pop": int(island_pop),
                "evo_steps": int(evo_steps),
                "sade_gen": int(sade_gen),
            },
            stages,
            pressure_cascade=cls.funnel_pressure_cascade_config(
                optimizer_strategy,
            ),
        )
        return funnel.to_dict()

    @staticmethod
    def _percentile(sorted_values, fraction):
        if not sorted_values:
            return None
        index = int(round(float(fraction) * (len(sorted_values) - 1)))
        return sorted_values[max(0, min(len(sorted_values) - 1, index))]

    @classmethod
    def tof_boundary_pressure_actions(
        cls, bounds, tof_vectors, near_fraction=0.03, min_pressure_count=2,
    ):
        """Return one-sided pressure actions for top candidate leg TOFs."""
        actions = []
        if not tof_vectors:
            return actions
        for leg_index, (low, high) in enumerate(bounds):
            span = max(float(high) - float(low), 1e-12)
            values = sorted(float(vector[leg_index]) for vector in tof_vectors)
            near_low = sum(
                (value - float(low)) / span <= float(near_fraction)
                for value in values
            )
            near_high = sum(
                (float(high) - value) / span <= float(near_fraction)
                for value in values
            )
            p05 = cls._percentile(values, 0.05)
            p95 = cls._percentile(values, 0.95)
            low_pressed = (
                near_low >= int(min_pressure_count)
                or (
                    p05 is not None
                    and (p05 - float(low)) / span <= float(near_fraction)
                )
            )
            high_pressed = (
                near_high >= int(min_pressure_count)
                or (
                    p95 is not None
                    and (float(high) - p95) / span <= float(near_fraction)
                )
            )
            if low_pressed:
                actions.append({
                    "leg": leg_index + 1,
                    "side": "low",
                    "near_count": int(near_low),
                    "p05": p05,
                    "p95": p95,
                })
            if high_pressed:
                actions.append({
                    "leg": leg_index + 1,
                    "side": "high",
                    "near_count": int(near_high),
                    "p05": p05,
                    "p95": p95,
                })
        return actions

    @staticmethod
    def adjust_leg_tof_bounds(
        bounds,
        actions,
        relax_fraction=0.20,
        min_relax_days=5.0,
        max_relax_days=60.0,
        min_leg_lower_days=1.0,
        adjustment_mode="widen",
    ):
        """Apply one-sided widen/shift corrections to leg TOF bounds."""
        adjusted = [[float(low), float(high)] for low, high in bounds]
        deltas = {}
        modes = {}
        actions_by_leg = {}
        for action in actions:
            actions_by_leg.setdefault(int(action["leg"]), set()).add(
                str(action["side"])
            )

        for leg, sides in sorted(actions_by_leg.items()):
            index = leg - 1
            low, high = adjusted[index]
            span = max(high - low, 1e-12)
            delta = max(
                float(min_relax_days),
                min(float(max_relax_days), span * float(relax_fraction)),
            )
            if adjustment_mode == "shift" and sides == {"low"}:
                actual_delta = min(delta, low - float(min_leg_lower_days))
                if actual_delta > 0:
                    adjusted[index][0] = low - actual_delta
                    adjusted[index][1] = high - actual_delta
                else:
                    adjusted[index][0] = max(
                        float(min_leg_lower_days), low - delta
                    )
                deltas[(leg, "low")] = actual_delta if actual_delta > 0 else delta
                modes[(leg, "low")] = "shift"
            elif adjustment_mode == "shift" and sides == {"high"}:
                adjusted[index][0] = low + delta
                adjusted[index][1] = high + delta
                deltas[(leg, "high")] = delta
                modes[(leg, "high")] = "shift"
            else:
                if "low" in sides:
                    adjusted[index][0] = max(
                        float(min_leg_lower_days), low - delta
                    )
                    deltas[(leg, "low")] = delta
                    modes[(leg, "low")] = "widen"
                if "high" in sides:
                    adjusted[index][1] = high + delta
                    deltas[(leg, "high")] = delta
                    modes[(leg, "high")] = "widen"
        return adjusted, deltas, modes

    @staticmethod
    def pressure_branch_decision(
        actions,
        policy="l1_combined",
        l1_improvement=None,
        max_l1_improvement=0.15,
        l1_best_to_l0_median=None,
        min_l1_best_to_l0_median=0.725,
    ):
        """Return whether a pressured L0 should spawn an adjusted branch."""
        if not actions:
            return False
        if policy == "pressure":
            return True
        if policy == "l1_improvement":
            return (
                l1_improvement is not None
                and float(l1_improvement) <= float(max_l1_improvement)
            )
        if policy == "l1_combined":
            poor_improvement = (
                l1_improvement is not None
                and float(l1_improvement) <= float(max_l1_improvement)
            )
            high_relative_l1 = (
                l1_best_to_l0_median is not None
                and float(l1_best_to_l0_median)
                >= float(min_l1_best_to_l0_median)
            )
            return bool(poor_improvement or high_relative_l1)
        raise ValueError("Unknown pressure branch policy: " + str(policy))

    @staticmethod
    def pressure_retry_decision(
        baseline_dv,
        current_best_dv,
        relaxed_dv,
        min_improvement=0.10,
    ):
        """Return True when same-seed rescue remains suspicious."""
        if relaxed_dv is None:
            return True
        if float(baseline_dv) <= 0.0:
            return False
        improvement = (
            float(baseline_dv) - float(current_best_dv)
        ) / float(baseline_dv)
        return improvement <= float(min_improvement)

    @staticmethod
    def select_exact_diverse_seeds(problem, genes, count):
        """Select exact-fitness seeds while retaining basin diversity."""
        count = int(count)
        if count <= 0:
            return []
        scored = sorted(
            ((float(problem.fitness(gene)[0]), list(gene)) for gene in genes),
            key=lambda item: item[0],
        )
        count = min(count, len(scored))
        if count <= 1:
            return [scored[0][1]] if scored else []
        candidate_pool = scored[:max(count, min(len(scored), count * 8))]
        lower, upper = problem.get_bounds()
        span = np.maximum(np.asarray(upper) - np.asarray(lower), 1e-12)
        selected = [candidate_pool.pop(0)]
        while len(selected) < count and candidate_pool:
            selected_genes = [np.asarray(item[1]) for item in selected]
            next_index = max(
                range(len(candidate_pool)),
                key=lambda index: min(
                    np.linalg.norm(
                        (np.asarray(candidate_pool[index][1]) - gene) / span
                    )
                    for gene in selected_genes
                ),
            )
            selected.append(candidate_pool.pop(next_index))
        return [item[1] for item in selected]

    @staticmethod
    def encounter_phase_embedding(bodies_by_name, context, problem, gene):
        """Embed a candidate in physical encounter phase space."""
        bodies = json.loads(context["bodies_json"])
        n_legs = len(bodies) - 1
        if context.get("tof_encoding", "direct") == "alpha":
            leg_tofs = alpha_leg_tofs(gene, n_legs)
        else:
            leg_tofs = direct_leg_tofs(gene, n_legs)
        encounter_epochs = [float(gene[0])]
        for leg_tof in leg_tofs:
            encounter_epochs.append(encounter_epochs[-1] + float(leg_tof))

        features = []
        for body_name, epoch in zip(bodies, encounter_epochs):
            position, velocity = bodies_by_name[body_name].eph(epoch)
            position = np.asarray(position, dtype=float)
            velocity = np.asarray(velocity, dtype=float)
            features.extend(position / max(np.linalg.norm(position), 1e-12))
            features.extend(velocity / max(np.linalg.norm(velocity), 1e-12))

        if context.get("tof_encoding", "direct") == "direct":
            lower, upper = problem.get_bounds()
            for leg_index, leg_tof in enumerate(leg_tofs):
                gene_index = 5 + 4 * leg_index
                span = max(
                    float(upper[gene_index] - lower[gene_index]), 1e-12
                )
                features.append(
                    (float(leg_tof) - lower[gene_index]) / span
                )
        else:
            total_tof = max(float(sum(leg_tofs)), 1e-12)
            features.extend(float(leg_tof) / total_tof for leg_tof in leg_tofs)
        return np.asarray(features, dtype=float)

    @staticmethod
    def handoff_diversity_embedding(bodies_by_name, context, problem, gene):
        """Embed cheap physical and decision-space handoff diversity.

        The feature families are balanced so that the many Cartesian phase
        components do not drown out launch epoch, leg TOFs, DSM locations, or
        flyby controls. Flyby beta and periapsis radius are intentionally used
        as cheap proxies for B-plane / v-infinity geometry: reconstructing
        physical incoming and outgoing v-infinity vectors would repeat every
        Lambert solve during each funnel handoff.
        """
        bodies = json.loads(context["bodies_json"])
        n_legs = len(bodies) - 1
        if context.get("tof_encoding", "direct") == "alpha":
            leg_tofs = alpha_leg_tofs(gene, n_legs)
        else:
            leg_tofs = direct_leg_tofs(gene, n_legs)
        lower, upper = problem.get_bounds()
        lower = np.asarray(lower, dtype=float)
        upper = np.asarray(upper, dtype=float)

        def normalized_gene(index):
            span = max(float(upper[index] - lower[index]), 1e-12)
            return (float(gene[index]) - float(lower[index])) / span

        features = []

        def extend_balanced(values):
            values = [float(value) for value in values]
            if values:
                scale = np.sqrt(float(len(values)))
                features.extend(value / scale for value in values)

        # Raw epoch separation remains useful when two windows have similar
        # orbital phases but belong to different synodic cycles.
        extend_balanced([normalized_gene(0)])

        encounter_epochs = [float(gene[0])]
        for leg_tof in leg_tofs:
            encounter_epochs.append(encounter_epochs[-1] + float(leg_tof))
        phase_features = []
        for body_name, epoch in zip(bodies, encounter_epochs):
            position, velocity = bodies_by_name[body_name].eph(epoch)
            position = np.asarray(position, dtype=float)
            velocity = np.asarray(velocity, dtype=float)
            phase_features.extend(
                position / max(np.linalg.norm(position), 1e-12)
            )
            phase_features.extend(
                velocity / max(np.linalg.norm(velocity), 1e-12)
            )
        extend_balanced(phase_features)

        if context.get("tof_encoding", "direct") == "direct":
            tof_features = [
                normalized_gene(5 + 4 * leg_index)
                for leg_index in range(n_legs)
            ]
        else:
            total_tof = max(float(sum(leg_tofs)), 1e-12)
            tof_features = [float(leg_tof) / total_tof for leg_tof in leg_tofs]
        extend_balanced(tof_features)

        extend_balanced([
            normalized_gene(4 + 4 * leg_index)
            for leg_index in range(n_legs)
        ])

        flyby_features = []
        for flyby_index in range(max(0, n_legs - 1)):
            flyby_features.extend([
                normalized_gene(6 + 4 * flyby_index),
                normalized_gene(7 + 4 * flyby_index),
            ])
        extend_balanced(flyby_features)
        return np.asarray(features, dtype=float)

    @staticmethod
    def select_phase_diverse_elites(
        problem, genes, count, embedding_for_gene, elite_fraction=0.35,
    ):
        """Keep good candidates while maximizing encounter-phase separation."""
        candidates = sorted(
            ((float(problem.fitness(gene)[0]), list(gene)) for gene in genes),
            key=lambda item: item[0],
        )
        count = min(int(count), len(candidates))
        if not count:
            return []
        pool_size = max(count, int(np.ceil(len(candidates) * elite_fraction)))
        pool = candidates[:min(len(candidates), pool_size)]
        embeddings = np.asarray([
            embedding_for_gene(gene) for _, gene in pool
        ])
        selected_indices = [0]
        min_distances = np.sum((embeddings - embeddings[0]) ** 2, axis=1)
        min_distances[0] = -1.0
        while len(selected_indices) < count:
            index = int(np.argmax(min_distances))
            selected_indices.append(index)
            distance = np.sum((embeddings - embeddings[index]) ** 2, axis=1)
            min_distances = np.minimum(min_distances, distance)
            min_distances[selected_indices] = -1.0
        return [pool[index][1] for index in selected_indices]

    @staticmethod
    def pareto_quality_novelty_front(quality, novelty, indices=None):
        """Return the non-dominated min-quality / max-novelty indices."""
        if indices is None:
            indices = range(len(quality))
        indices = list(indices)
        front = []
        for index in indices:
            dominated = any(
                quality[other] <= quality[index]
                and novelty[other] >= novelty[index]
                and (
                    quality[other] < quality[index]
                    or novelty[other] > novelty[index]
                )
                for other in indices
                if other != index
            )
            if not dominated:
                front.append(index)
        return front

    @classmethod
    def select_pareto_diverse_elites(
        cls, problem, genes, count, embedding_for_gene,
        elite_fraction=0.35, quality_weight=0.5,
    ):
        """Select quality-gated elites on a quality/novelty Pareto front.

        The absolute champion is always retained. Each following choice uses
        its distance to the already selected set as novelty, computes the
        non-dominated front, and chooses the front point closest to the
        normalized quality/novelty utopia. Ranking quality makes the policy
        independent from sequence-specific Delta-V scales and penalties.
        """
        count = int(count)
        if count <= 0:
            return []
        quality_weight = float(quality_weight)
        if not 0.0 <= quality_weight <= 1.0:
            raise ValueError("quality_weight must be between 0 and 1")
        candidates = sorted(
            ((float(problem.fitness(gene)[0]), list(gene)) for gene in genes),
            key=lambda item: item[0],
        )
        count = min(count, len(candidates))
        if not count:
            return []
        pool_size = max(count, int(np.ceil(len(candidates) * elite_fraction)))
        pool = candidates[:min(len(candidates), pool_size)]
        embeddings = np.asarray([
            embedding_for_gene(gene) for _, gene in pool
        ], dtype=float)
        if embeddings.ndim != 2:
            raise ValueError("handoff embeddings must have a stable dimension")

        selected_indices = [0]
        available = set(range(1, len(pool)))
        min_distances = np.sum((embeddings - embeddings[0]) ** 2, axis=1)
        quality = np.arange(len(pool), dtype=float)
        if len(pool) > 1:
            quality /= float(len(pool) - 1)

        while len(selected_indices) < count and available:
            novelty_scale = max(
                max(float(min_distances[index]) for index in available),
                1e-12,
            )
            novelty = min_distances / novelty_scale
            front = cls.pareto_quality_novelty_front(
                quality, novelty, indices=available,
            )
            next_index = min(
                front,
                key=lambda index: (
                    quality_weight * quality[index] ** 2
                    + (1.0 - quality_weight) * (1.0 - novelty[index]) ** 2,
                    quality[index],
                    -novelty[index],
                    index,
                ),
            )
            selected_indices.append(next_index)
            available.remove(next_index)
            distance = np.sum(
                (embeddings - embeddings[next_index]) ** 2, axis=1
            )
            min_distances = np.minimum(min_distances, distance)
        return [pool[index][1] for index in selected_indices]

    @staticmethod
    def _quality_density_clusters(embeddings, cluster_count, iterations=5):
        """Return deterministic k-medoids-like cluster memberships."""
        embeddings = np.asarray(embeddings, dtype=float)
        if embeddings.ndim != 2:
            raise ValueError("cluster embeddings must have a stable dimension")
        n_candidates = len(embeddings)
        cluster_count = min(max(1, int(cluster_count)), n_candidates)
        medoids = [0]
        min_distances = np.sum(
            (embeddings - embeddings[0]) ** 2, axis=1
        )
        while len(medoids) < cluster_count:
            min_distances[medoids] = -1.0
            index = int(np.argmax(min_distances))
            medoids.append(index)
            distance = np.sum(
                (embeddings - embeddings[index]) ** 2, axis=1
            )
            min_distances = np.minimum(min_distances, distance)

        assignments = np.zeros(n_candidates, dtype=int)
        for _ in range(max(1, int(iterations))):
            distances = np.asarray([
                np.sum((embeddings - embeddings[index]) ** 2, axis=1)
                for index in medoids
            ]).T
            assignments = np.argmin(distances, axis=1)
            updated = []
            for cluster_index in range(cluster_count):
                members = np.flatnonzero(assignments == cluster_index)
                if not len(members):
                    updated.append(medoids[cluster_index])
                    continue
                member_vectors = embeddings[members]
                pairwise = np.sum(
                    (member_vectors[:, None, :] - member_vectors[None, :, :]) ** 2,
                    axis=2,
                )
                costs = np.sum(pairwise, axis=1)
                # Members retain global quality order, so the first minimum
                # is also the best-quality medoid when geometric costs tie.
                updated.append(int(members[int(np.argmin(costs))]))
            if updated == medoids:
                break
            medoids = updated

        distances = np.asarray([
            np.sum((embeddings - embeddings[index]) ** 2, axis=1)
            for index in medoids
        ]).T
        assignments = np.argmin(distances, axis=1)
        return [
            list(map(int, np.flatnonzero(assignments == cluster_index)))
            for cluster_index in range(cluster_count)
        ]

    @classmethod
    def select_basin_allocated_elites(
        cls, problem, genes, count, embedding_for_gene,
        elite_fraction=0.75, max_clusters=8,
    ):
        """Allocate L0 handoff slots to good, populated phase-space basins.

        Every discovered basin receives one representative. Remaining slots
        are distributed with a D'Hondt-style divisor over a rank-based basin
        mass, so several independently good candidates provide stronger
        evidence than one isolated elite without requiring sequence-specific
        Delta-V thresholds.
        """
        count = int(count)
        if count <= 0:
            return []
        candidates = sorted(
            ((float(problem.fitness(gene)[0]), list(gene)) for gene in genes),
            key=lambda item: item[0],
        )
        count = min(count, len(candidates))
        if not count:
            return []
        pool_size = max(count, int(np.ceil(len(candidates) * elite_fraction)))
        pool = candidates[:min(len(candidates), pool_size)]
        embeddings = np.asarray([
            embedding_for_gene(gene) for _, gene in pool
        ], dtype=float)
        cluster_count = min(int(max_clusters), count, len(pool))
        clusters = cls._quality_density_clusters(
            embeddings, cluster_count,
        )

        allocations = [1] * len(clusters)
        masses = [
            sum(1.0 / np.sqrt(float(rank + 1)) for rank in members)
            for members in clusters
        ]
        remaining = count - len(clusters)
        while remaining > 0:
            eligible = [
                index for index, members in enumerate(clusters)
                if allocations[index] < len(members)
            ]
            if not eligible:
                break
            cluster_index = max(
                eligible,
                key=lambda index: (
                    masses[index] / float(allocations[index] + 1),
                    -min(clusters[index]),
                    -index,
                ),
            )
            allocations[cluster_index] += 1
            remaining -= 1

        selected_indices = []
        for members, allocation in zip(clusters, allocations):
            selected_indices.extend(members[:allocation])
        selected_indices.sort()
        return [pool[index][1] for index in selected_indices[:count]]

    @staticmethod
    def normalized_decision_embedding(problem, gene, periodic_indices=()):
        """Return box-normalized genes with circular angular coordinates."""
        lower, upper = problem.get_bounds()
        periodic_indices = set(int(index) for index in periodic_indices)
        features = []
        for index, value in enumerate(gene):
            if index in periodic_indices:
                features.extend([np.cos(float(value)), np.sin(float(value))])
                continue
            span = max(float(upper[index] - lower[index]), 1e-12)
            features.append(
                (float(value) - float(lower[index])) / span
            )
        return np.asarray(features, dtype=float)

    @staticmethod
    def _interpolate_gene(
        left, right, fraction, lower, upper, periodic_indices=(),
    ):
        periodic_indices = set(int(index) for index in periodic_indices)
        result = []
        for index, (left_value, right_value) in enumerate(zip(left, right)):
            left_value = float(left_value)
            right_value = float(right_value)
            if index in periodic_indices:
                delta = (right_value - left_value + np.pi) % (2.0 * np.pi) - np.pi
                value = left_value + float(fraction) * delta
                while value < float(lower[index]):
                    value += 2.0 * np.pi
                while value > float(upper[index]):
                    value -= 2.0 * np.pi
            else:
                value = left_value + float(fraction) * (
                    right_value - left_value
                )
            result.append(float(np.clip(value, lower[index], upper[index])))
        return result

    @classmethod
    def select_hill_valley_elites(
        cls, problem, genes, count, periodic_indices=(), elite_fraction=0.35,
        max_neighbors=None, max_test_points=3,
        barrier_relative_tolerance=0.0, valley_slot_fraction=1.0,
        return_diagnostics=False,
    ):
        """Cluster a quality-gated sample with nearest-better valley tests."""
        count = int(count)
        if count <= 0:
            empty = ([], {"pool_size": 0, "cluster_count": 0,
                          "hill_valley_evaluations": 0})
            return empty if return_diagnostics else empty[0]
        elite_fraction = float(elite_fraction)
        if not 0.0 < elite_fraction <= 1.0:
            raise ValueError("elite_fraction must be in (0, 1]")
        max_test_points = max(1, int(max_test_points))
        barrier_relative_tolerance = max(
            0.0, float(barrier_relative_tolerance)
        )
        valley_slot_fraction = float(valley_slot_fraction)
        if not 0.0 < valley_slot_fraction <= 1.0:
            raise ValueError("valley_slot_fraction must be in (0, 1]")

        unique = {}
        candidate_evaluations = 0
        for gene in genes:
            normalized = tuple(float(value) for value in gene)
            if normalized in unique:
                continue
            fitness = float(problem.fitness(gene)[0])
            candidate_evaluations += 1
            unique[normalized] = (fitness, list(gene))
        candidates = sorted(unique.values(), key=lambda item: item[0])
        count = min(count, len(candidates))
        if not count:
            empty = ([], {"pool_size": 0, "cluster_count": 0,
                          "candidate_evaluations": candidate_evaluations,
                          "hill_valley_evaluations": 0})
            return empty if return_diagnostics else empty[0]
        pool_size = max(count, int(np.ceil(len(candidates) * elite_fraction)))
        pool = candidates[:min(len(candidates), pool_size)]
        embeddings = np.asarray([
            cls.normalized_decision_embedding(
                problem, gene, periodic_indices=periodic_indices,
            )
            for _, gene in pool
        ])
        dimension = max(1, int(problem.get_nx()))
        neighbor_limit = (
            dimension + 1 if max_neighbors is None else max(1, int(max_neighbors))
        )
        expected_edge_length = (1.0 / float(len(pool))) ** (1.0 / dimension)
        lower, upper = problem.get_bounds()
        clusters = [[0]]
        cluster_by_index = {0: 0}
        hill_valley_evaluations = 0
        hill_valley_tests = 0

        for index in range(1, len(pool)):
            distances = np.sum(
                (embeddings[:index] - embeddings[index]) ** 2, axis=1
            )
            nearest_better = np.argsort(distances)
            checked_clusters = set()
            assigned_cluster = None
            for neighbor in nearest_better[:min(index, neighbor_limit)]:
                neighbor = int(neighbor)
                cluster_index = cluster_by_index[neighbor]
                if cluster_index in checked_clusters:
                    continue
                checked_clusters.add(cluster_index)
                edge_length = float(np.sqrt(max(distances[neighbor], 0.0)))
                test_points = min(
                    max_test_points,
                    1 + int(np.floor(
                        edge_length / max(expected_edge_length, 1e-12)
                    )),
                )
                hill_valley_tests += 1
                same_valley = True
                threshold = max(pool[index][0], pool[neighbor][0])
                tolerance = max(
                    1e-10 * max(1.0, abs(float(threshold))),
                    barrier_relative_tolerance * max(1.0, abs(float(threshold))),
                )
                for test_index in range(1, test_points + 1):
                    fraction = test_index / float(test_points + 1)
                    test_gene = cls._interpolate_gene(
                        pool[neighbor][1], pool[index][1], fraction,
                        lower, upper, periodic_indices=periodic_indices,
                    )
                    test_fitness = float(problem.fitness(test_gene)[0])
                    hill_valley_evaluations += 1
                    if test_fitness > threshold + tolerance:
                        same_valley = False
                        break
                if same_valley:
                    assigned_cluster = cluster_index
                    break
            if assigned_cluster is None:
                assigned_cluster = len(clusters)
                clusters.append([])
            clusters[assigned_cluster].append(index)
            cluster_by_index[index] = assigned_cluster

        # Preserve every detected valley before assigning additional handoff
        # slots by pure fitness. This keeps isolated deep candidates without
        # using basin density as a quality proxy.
        valley_slots = min(
            len(clusters), count,
            max(1, int(np.ceil(count * valley_slot_fraction))),
        )
        selected_indices = [
            clusters[index][0] for index in range(valley_slots)
        ]
        selected_set = set(selected_indices)
        family_candidates = sorted(
            member
            for cluster in clusters[:valley_slots]
            for member in cluster[1:]
        )
        for index in family_candidates:
            if len(selected_indices) >= count:
                break
            selected_indices.append(index)
            selected_set.add(index)
        # If selected valleys contain too few family members, retain further
        # strict valley roots before falling back to global quality.
        for cluster in clusters[valley_slots:]:
            if len(selected_indices) >= count:
                break
            index = cluster[0]
            selected_indices.append(index)
            selected_set.add(index)
        if len(selected_indices) < count:
            selected_indices.extend(
                index for index in range(len(pool))
                if index not in selected_set
                and len(selected_indices) < count
            )
        selected_indices.sort()
        selected = [pool[index][1] for index in selected_indices]
        diagnostics = {
            "candidate_evaluations": candidate_evaluations,
            "pool_size": len(pool),
            "cluster_count": len(clusters),
            "cluster_sizes": [len(members) for members in clusters],
            "hill_valley_tests": hill_valley_tests,
            "hill_valley_evaluations": hill_valley_evaluations,
            "barrier_relative_tolerance": barrier_relative_tolerance,
            "valley_slot_fraction": valley_slot_fraction,
            "valley_slots": valley_slots,
            "family_slots": max(0, len(selected_indices) - valley_slots),
        }
        return (selected, diagnostics) if return_diagnostics else selected

    @staticmethod
    def update_exact_phase_archive(
        archive_genes, candidate_genes, max_size, select_elites,
    ):
        """Retain exact-good candidates while preserving encounter diversity."""
        unique = {}
        for gene in list(archive_genes) + list(candidate_genes):
            normalized = tuple(float(value) for value in gene)
            unique[normalized] = list(gene)
        if not unique:
            return []
        return select_elites(
            list(unique.values()), min(int(max_size), len(unique))
        )
