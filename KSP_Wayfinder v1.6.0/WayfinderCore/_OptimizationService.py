"""Technical helpers for optimizer orchestration.

This is deliberately a narrow extraction from ``Wayfinder``. Public entry
points remain on the facade while optimizer policy moves into testable helpers.
"""

import json
import os

import numpy as np
import pygmo as pg

from _Optimization import alpha_leg_tofs, direct_leg_tofs


class OptimizationService:
    """Build topologies and report archipelago behavior."""

    FUNNEL_STRATEGIES = {
        "funnel": "legacy",
        "funnel_local_exact": "local",
        "funnel_hybrid_exact": "hybrid",
        "funnel_phase_elites_nm": "phase_elites_nm",
        "funnel_phase_elites_equal": "phase_elites_nm_equal",
        "funnel_scout_archive": "scout_archive_nm",
        "funnel_scout_archive_32": "scout_archive_nm_32",
        "funnel_scout_archive_64": "scout_archive_nm_64",
        "funnel_scout_archive_128": "scout_archive_nm_128",
    }

    @classmethod
    def funnel_exact_strategy(cls, optimizer_strategy):
        return cls.FUNNEL_STRATEGIES.get(optimizer_strategy)

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
                wide_step_cost = max(
                    1, n_stage_1 * island_pop * int(sade_gen)
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
        return {
            "kind": "funnel",
            "optimizer_strategy": str(optimizer_strategy),
            "exact_strategy": exact_strategy,
            "requested": {
                "n_islands": int(n_islands),
                "island_pop": int(island_pop),
                "evo_steps": int(evo_steps),
                "sade_gen": int(sade_gen),
            },
            "stages": cls.funnel_stage_plan(
                n_islands,
                island_pop,
                evo_steps,
                sade_gen,
                exact_strategy=exact_strategy,
            ),
        }

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
