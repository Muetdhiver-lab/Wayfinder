# -*- coding: utf-8 -*-
"""SQLite datastore for Wayfinder jobs and results.

The store keeps wildcard search templates, concrete sequences and optimized
jobs separate so that results can be compared across batches and binning
strategies.
"""

import hashlib
import json
import os
import sqlite3
import socket
import uuid
from datetime import datetime, timedelta


SCHEMA_VERSION = 15


def _json_dumps(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def _json_loads(value, default=None):
    if value in (None, ""):
        return default
    return json.loads(value)


def _jsonable(value):
    if hasattr(value, "tolist"):
        return _jsonable(value.tolist())
    if isinstance(value, tuple):
        return [_jsonable(item) for item in value]
    if isinstance(value, list):
        return [_jsonable(item) for item in value]
    if isinstance(value, dict):
        return {key: _jsonable(item) for key, item in value.items()}
    return value


def _utc_now():
    return datetime.utcnow().isoformat(timespec="microseconds") + "Z"


def _utc_after(seconds):
    return (
        datetime.utcnow() + timedelta(seconds=max(1, int(seconds)))
    ).isoformat(timespec="microseconds") + "Z"


def _default_worker_id():
    return "{}:{}:{}".format(socket.gethostname(), os.getpid(), uuid.uuid4().hex)


def _sqlite_text(value):
    if value is None:
        return None
    return str(value)


def _param_hash(params):
    return hashlib.sha256(_json_dumps(params).encode("utf-8")).hexdigest()


class SQLiteJobStore:
    def __init__(self, path):
        self.path = str(path)
        self.conn = sqlite3.connect(self.path)
        self.conn.row_factory = sqlite3.Row
        self.initialize()

    def close(self):
        self.conn.close()

    def initialize(self):
        cur = self.conn.cursor()
        cur.executescript(
            """
            PRAGMA foreign_keys = ON;

            CREATE TABLE IF NOT EXISTS metadata (
                key TEXT PRIMARY KEY,
                value TEXT NOT NULL
            );

            CREATE TABLE IF NOT EXISTS batches (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                name TEXT NOT NULL,
                planet_pack TEXT NOT NULL,
                template_json TEXT,
                generation_options_json TEXT,
                purpose TEXT NOT NULL DEFAULT 'production',
                created_at TEXT NOT NULL,
                UNIQUE(name, planet_pack, template_json, generation_options_json)
            );

            CREATE TABLE IF NOT EXISTS sequences (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                planet_pack TEXT NOT NULL,
                short_name TEXT NOT NULL,
                body_path TEXT NOT NULL,
                bodies_json TEXT NOT NULL,
                start_body TEXT NOT NULL,
                target_body TEXT NOT NULL,
                flyby_count INTEGER NOT NULL,
                UNIQUE(planet_pack, body_path)
            );

            CREATE TABLE IF NOT EXISTS sequence_bodies (
                sequence_id INTEGER NOT NULL,
                position INTEGER NOT NULL,
                body_name TEXT NOT NULL,
                role TEXT NOT NULL,
                PRIMARY KEY(sequence_id, position),
                FOREIGN KEY(sequence_id) REFERENCES sequences(id) ON DELETE CASCADE
            );

            CREATE TABLE IF NOT EXISTS jobs (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                sequence_id INTEGER NOT NULL,
                planet_pack TEXT NOT NULL,
                t0_min REAL NOT NULL,
                t0_max REAL NOT NULL,
                tof_min REAL NOT NULL,
                tof_max REAL NOT NULL,
                vinf_min REAL NOT NULL,
                vinf_max REAL NOT NULL,
                tof_encoding TEXT NOT NULL,
                tof_profile TEXT,
                leg_tof_bounds_json TEXT,
                opt_level TEXT NOT NULL,
                opt_injection TEXT NOT NULL,
                arrival_mode TEXT NOT NULL DEFAULT 'legacy',
                ejection_altitude REAL NOT NULL,
                injection_altitude REAL NOT NULL,
                rp_target REAL,
                e_target REAL,
                add_vinf_dep INTEGER NOT NULL,
                add_vinf_arr INTEGER NOT NULL,
                orbit_insertion INTEGER NOT NULL,
                multi_objective INTEGER NOT NULL,
                lambert_max_revs INTEGER NOT NULL,
                sade_gen INTEGER,
                n_island INTEGER,
                island_pop INTEGER,
                n_evo_steps INTEGER,
                optimizer_topology TEXT,
                optimizer_seed INTEGER,
                status TEXT NOT NULL,
                claimed_at TEXT,
                claim_expires_at TEXT,
                worker_id TEXT,
                param_hash TEXT NOT NULL UNIQUE,
                created_at TEXT NOT NULL,
                updated_at TEXT NOT NULL,
                FOREIGN KEY(sequence_id) REFERENCES sequences(id)
            );

            CREATE TABLE IF NOT EXISTS batch_jobs (
                batch_id INTEGER NOT NULL,
                job_id INTEGER NOT NULL,
                PRIMARY KEY(batch_id, job_id),
                FOREIGN KEY(batch_id) REFERENCES batches(id) ON DELETE CASCADE,
                FOREIGN KEY(job_id) REFERENCES jobs(id) ON DELETE CASCADE
            );

            CREATE TABLE IF NOT EXISTS runs (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                job_id INTEGER NOT NULL,
                pykep_version TEXT,
                pygmo_version TEXT,
                wayfinder_version TEXT,
                status TEXT NOT NULL,
                started_at TEXT,
                ended_at TEXT,
                optimizer_topology TEXT,
                optimizer_seed INTEGER,
                effective_optimizer_seed INTEGER,
                requested_n_island INTEGER,
                actual_n_island INTEGER,
                island_pop INTEGER,
                sade_gen INTEGER,
                n_evo_steps INTEGER,
                actual_n_evo_steps INTEGER,
                adaptive_stop_json TEXT,
                optimizer_strategy TEXT,
                funnel_config_json TEXT,
                funnel_config_hash TEXT,
                code_revision TEXT,
                planet_pack_hash TEXT,
                stop_reason TEXT,
                runtime_seconds REAL,
                notes TEXT,
                error TEXT,
                FOREIGN KEY(job_id) REFERENCES jobs(id) ON DELETE CASCADE
            );

            CREATE TABLE IF NOT EXISTS results (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                run_id INTEGER NOT NULL UNIQUE,
                objective_dv REAL NOT NULL,
                result_t0 REAL NOT NULL,
                result_tof REAL NOT NULL,
                ejection_vinf REAL NOT NULL,
                arrival_vinf REAL,
                FOREIGN KEY(run_id) REFERENCES runs(id) ON DELETE CASCADE
            );

            CREATE TABLE IF NOT EXISTS genes (
                run_id INTEGER PRIMARY KEY,
                gene_json TEXT NOT NULL,
                FOREIGN KEY(run_id) REFERENCES runs(id) ON DELETE CASCADE
            );

            CREATE TABLE IF NOT EXISTS optimizer_snapshots (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                run_id INTEGER NOT NULL,
                step INTEGER NOT NULL,
                best_fitness REAL NOT NULL,
                average_fitness REAL,
                champion_fitness_json TEXT NOT NULL,
                champion_gene_json TEXT NOT NULL,
                island_count INTEGER,
                topology_vertices INTEGER,
                topology_edges INTEGER,
                migrants_published INTEGER,
                migration_islands_active INTEGER,
                migrations_accepted INTEGER,
                created_at TEXT NOT NULL,
                FOREIGN KEY(run_id) REFERENCES runs(id) ON DELETE CASCADE,
                UNIQUE(run_id, step)
            );

            CREATE TABLE IF NOT EXISTS optimizer_population_points (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                run_id INTEGER NOT NULL,
                step INTEGER NOT NULL,
                island_index INTEGER NOT NULL,
                individual_index INTEGER NOT NULL,
                fitness REAL NOT NULL,
                fitness_json TEXT NOT NULL,
                gene_json TEXT NOT NULL,
                source TEXT NOT NULL,
                created_at TEXT NOT NULL,
                FOREIGN KEY(run_id) REFERENCES runs(id) ON DELETE CASCADE,
                UNIQUE(run_id, step, island_index, individual_index, source)
            );

            CREATE TABLE IF NOT EXISTS optimizer_stages (
                run_id INTEGER NOT NULL,
                stage_index INTEGER NOT NULL,
                stage_name TEXT NOT NULL,
                n_island INTEGER NOT NULL,
                island_pop INTEGER NOT NULL,
                sade_gen INTEGER NOT NULL,
                planned_evo_steps INTEGER NOT NULL,
                actual_evo_steps INTEGER NOT NULL,
                ejection_model TEXT NOT NULL,
                initialization TEXT,
                topology_name TEXT,
                migration_rate INTEGER,
                exact_archive_size INTEGER,
                adaptive_stop_json TEXT,
                algorithms_json TEXT NOT NULL,
                topology_vertices INTEGER,
                topology_edges INTEGER,
                migrants_published INTEGER,
                migration_islands_active INTEGER,
                migrations_accepted INTEGER,
                best_fitness REAL NOT NULL,
                average_fitness REAL,
                runtime_seconds REAL NOT NULL,
                stop_reason TEXT NOT NULL,
                PRIMARY KEY(run_id, stage_index),
                FOREIGN KEY(run_id) REFERENCES runs(id) ON DELETE CASCADE
            );

            CREATE TABLE IF NOT EXISTS porkchop_samples (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                source_run_id INTEGER NOT NULL,
                sampler_name TEXT NOT NULL,
                t0 REAL NOT NULL,
                tof REAL NOT NULL,
                metric REAL,
                fitness REAL,
                ejection_vinf REAL,
                gene_json TEXT NOT NULL,
                decode_error TEXT,
                created_at TEXT NOT NULL,
                FOREIGN KEY(source_run_id) REFERENCES runs(id) ON DELETE CASCADE,
                UNIQUE(source_run_id, sampler_name, t0, tof)
            );

            CREATE TABLE IF NOT EXISTS porkchop_sampler_metadata (
                source_run_id INTEGER NOT NULL,
                sampler_name TEXT NOT NULL,
                metadata_json TEXT NOT NULL,
                created_at TEXT NOT NULL,
                updated_at TEXT NOT NULL,
                FOREIGN KEY(source_run_id) REFERENCES runs(id) ON DELETE CASCADE,
                PRIMARY KEY(source_run_id, sampler_name)
            );

            CREATE INDEX IF NOT EXISTS idx_jobs_window ON jobs(planet_pack, t0_min, t0_max);
            CREATE INDEX IF NOT EXISTS idx_sequences_route ON sequences(planet_pack, start_body, target_body);
            CREATE INDEX IF NOT EXISTS idx_results_objective ON results(objective_dv);
            CREATE INDEX IF NOT EXISTS idx_optimizer_snapshots_run ON optimizer_snapshots(run_id, step);
            CREATE INDEX IF NOT EXISTS idx_optimizer_population_run ON optimizer_population_points(run_id, step);
            CREATE INDEX IF NOT EXISTS idx_porkchop_samples_run ON porkchop_samples(source_run_id, sampler_name);
            """
        )
        self._ensure_column("batches", "purpose", "TEXT NOT NULL DEFAULT 'production'")
        self._ensure_column("jobs", "optimizer_topology", "TEXT")
        self._ensure_column("jobs", "optimizer_seed", "INTEGER")
        self._ensure_column("jobs", "tof_profile", "TEXT")
        self._ensure_column("jobs", "leg_tof_bounds_json", "TEXT")
        self._ensure_column("jobs", "arrival_mode", "TEXT NOT NULL DEFAULT 'legacy'")
        self._ensure_column("jobs", "claimed_at", "TEXT")
        self._ensure_column("jobs", "claim_expires_at", "TEXT")
        self._ensure_column("jobs", "worker_id", "TEXT")
        self._ensure_column("results", "arrival_vinf", "REAL")
        self._ensure_column("optimizer_snapshots", "average_fitness", "REAL")
        self._ensure_column("optimizer_snapshots", "island_count", "INTEGER")
        self._ensure_column("optimizer_snapshots", "topology_vertices", "INTEGER")
        self._ensure_column("optimizer_snapshots", "topology_edges", "INTEGER")
        self._ensure_column("optimizer_snapshots", "migrants_published", "INTEGER")
        self._ensure_column("optimizer_snapshots", "migration_islands_active", "INTEGER")
        self._ensure_column("optimizer_snapshots", "migrations_accepted", "INTEGER")
        self._ensure_column("runs", "optimizer_topology", "TEXT")
        self._ensure_column("runs", "optimizer_seed", "INTEGER")
        self._ensure_column("runs", "effective_optimizer_seed", "INTEGER")
        self._ensure_column("runs", "requested_n_island", "INTEGER")
        self._ensure_column("runs", "actual_n_island", "INTEGER")
        self._ensure_column("runs", "island_pop", "INTEGER")
        self._ensure_column("runs", "sade_gen", "INTEGER")
        self._ensure_column("runs", "n_evo_steps", "INTEGER")
        self._ensure_column("runs", "actual_n_evo_steps", "INTEGER")
        self._ensure_column("runs", "adaptive_stop_json", "TEXT")
        self._ensure_column("runs", "optimizer_strategy", "TEXT")
        self._ensure_column("runs", "funnel_config_json", "TEXT")
        self._ensure_column("runs", "funnel_config_hash", "TEXT")
        self._ensure_column("runs", "code_revision", "TEXT")
        self._ensure_column("runs", "planet_pack_hash", "TEXT")
        self._ensure_column("runs", "stop_reason", "TEXT")
        self._ensure_column("runs", "runtime_seconds", "REAL")
        self._ensure_column("optimizer_stages", "initialization", "TEXT")
        self._ensure_column("optimizer_stages", "topology_name", "TEXT")
        self._ensure_column("optimizer_stages", "migration_rate", "INTEGER")
        self._ensure_column("optimizer_stages", "exact_archive_size", "INTEGER")
        self._ensure_column("optimizer_stages", "adaptive_stop_json", "TEXT")
        self._ensure_column("optimizer_stages", "topology_vertices", "INTEGER")
        self._ensure_column("optimizer_stages", "topology_edges", "INTEGER")
        self._ensure_column("optimizer_stages", "migrants_published", "INTEGER")
        self._ensure_column("optimizer_stages", "migration_islands_active", "INTEGER")
        self._ensure_column("optimizer_stages", "migrations_accepted", "INTEGER")
        cur.execute(
            "CREATE INDEX IF NOT EXISTS idx_jobs_claim "
            "ON jobs(status, claim_expires_at)"
        )
        cur.execute(
            "INSERT OR REPLACE INTO metadata(key, value) VALUES (?, ?)",
            ("schema_version", str(SCHEMA_VERSION)),
        )
        self.conn.commit()

    def _ensure_column(self, table, column, definition):
        columns = {
            row["name"]
            for row in self.conn.execute(f"PRAGMA table_info({table})").fetchall()
        }
        if column not in columns:
            self.conn.execute(f"ALTER TABLE {table} ADD COLUMN {column} {definition}")

    def record_optimizer_population(self, run_id, step, populations, source="final"):
        now = _utc_now()
        values = []
        for island_index, population in enumerate(populations):
            genes = _jsonable(population.get_x())
            fitnesses = _jsonable(population.get_f())
            for individual_index, gene in enumerate(genes):
                fitness = fitnesses[individual_index]
                if isinstance(fitness, list):
                    fitness_value = float(fitness[0])
                else:
                    fitness_value = float(fitness)
                values.append((
                    int(run_id),
                    int(step),
                    int(island_index),
                    int(individual_index),
                    fitness_value,
                    _json_dumps(fitness),
                    _json_dumps(gene),
                    source,
                    now,
                ))
        self.conn.executemany(
            """
            INSERT INTO optimizer_population_points(
                run_id, step, island_index, individual_index, fitness,
                fitness_json, gene_json, source, created_at
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
            ON CONFLICT(run_id, step, island_index, individual_index, source) DO UPDATE SET
                fitness = excluded.fitness,
                fitness_json = excluded.fitness_json,
                gene_json = excluded.gene_json,
                created_at = excluded.created_at
            """,
            values,
        )
        self.conn.commit()

    def replace_porkchop_samples(self, source_run_id, sampler_name, rows, metadata=None):
        self.conn.execute(
            "DELETE FROM porkchop_samples WHERE source_run_id = ? AND sampler_name = ?",
            (int(source_run_id), sampler_name),
        )
        now = _utc_now()
        values = []
        for row in rows:
            values.append((
                int(source_run_id),
                sampler_name,
                float(row["t0"]),
                float(row["tof"]),
                None if row.get("metric") is None else float(row["metric"]),
                None if row.get("fitness") is None else float(row["fitness"]),
                None if row.get("ejection_vinf") is None else float(row["ejection_vinf"]),
                _json_dumps(row["gene"]),
                row.get("decode_error"),
                now,
            ))
        self.conn.executemany(
            """
            INSERT INTO porkchop_samples(
                source_run_id, sampler_name, t0, tof, metric, fitness,
                ejection_vinf, gene_json, decode_error, created_at
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            values,
        )
        if metadata is not None:
            self.upsert_porkchop_sampler_metadata(source_run_id, sampler_name, metadata, commit=False)
        self.conn.commit()

    def upsert_porkchop_sampler_metadata(self, source_run_id, sampler_name, metadata, commit=True):
        now = _utc_now()
        self.conn.execute(
            """
            INSERT INTO porkchop_sampler_metadata(
                source_run_id, sampler_name, metadata_json, created_at, updated_at
            ) VALUES (?, ?, ?, ?, ?)
            ON CONFLICT(source_run_id, sampler_name) DO UPDATE SET
                metadata_json = excluded.metadata_json,
                updated_at = excluded.updated_at
            """,
            (
                int(source_run_id),
                sampler_name,
                _json_dumps(_jsonable(metadata)),
                now,
                now,
            ),
        )
        if commit:
            self.conn.commit()

    def porkchop_sampler_metadata(self, source_run_id, sampler_name):
        row = self.conn.execute(
            """
            SELECT metadata_json, created_at, updated_at
            FROM porkchop_sampler_metadata
            WHERE source_run_id = ? AND sampler_name = ?
            """,
            (int(source_run_id), sampler_name),
        ).fetchone()
        if row is None:
            return None
        metadata = _json_loads(row["metadata_json"], default={})
        metadata["created_at"] = row["created_at"]
        metadata["updated_at"] = row["updated_at"]
        return metadata

    def porkchop_samples(self, source_run_id, sampler_name):
        if isinstance(sampler_name, (list, tuple, set)):
            sampler_names = list(sampler_name)
            if not sampler_names:
                return []
            placeholders = ", ".join(["?"] * len(sampler_names))
            rows = self.conn.execute(
                f"""
                SELECT
                    source_run_id,
                    sampler_name,
                    t0,
                    tof,
                    metric,
                    fitness,
                    ejection_vinf,
                    gene_json,
                    decode_error,
                    created_at
                FROM porkchop_samples
                WHERE source_run_id = ? AND sampler_name IN ({placeholders})
                ORDER BY sampler_name ASC, t0 ASC, tof ASC
                """,
                [int(source_run_id)] + sampler_names,
            ).fetchall()
        else:
            rows = self.conn.execute(
                """
            SELECT
                source_run_id,
                sampler_name,
                t0,
                tof,
                metric,
                fitness,
                ejection_vinf,
                gene_json,
                decode_error,
                created_at
            FROM porkchop_samples
            WHERE source_run_id = ? AND sampler_name = ?
            ORDER BY t0 ASC, tof ASC
            """,
                (int(source_run_id), sampler_name),
            ).fetchall()
        return [
            {
                "source_run_id": int(row["source_run_id"]),
                "sampler_name": row["sampler_name"],
                "result_t0": float(row["t0"]),
                "result_tof": float(row["tof"]),
                "result_DV": None if row["metric"] is None else float(row["metric"]),
                "fitness": None if row["fitness"] is None else float(row["fitness"]),
                "ejection_vinf": None if row["ejection_vinf"] is None else float(row["ejection_vinf"]),
                "gene": _json_loads(row["gene_json"], default=[]),
                "decode_error": row["decode_error"],
                "created_at": row["created_at"],
            }
            for row in rows
        ]

    def upsert_batch(self, name, planet_pack, template=None, generation_options=None, purpose="production"):
        template_json = _json_dumps(template) if template is not None else None
        options_json = _json_dumps(generation_options or {})
        now = _utc_now()
        self.conn.execute(
            """
            INSERT OR IGNORE INTO batches(
                name, planet_pack, template_json, generation_options_json, purpose, created_at
            ) VALUES (?, ?, ?, ?, ?, ?)
            """,
            (name, planet_pack, template_json, options_json, purpose, now),
        )
        self.conn.execute(
            """
            UPDATE batches SET purpose = ?
            WHERE name = ? AND planet_pack = ?
              AND template_json IS ? AND generation_options_json = ?
            """,
            (purpose, name, planet_pack, template_json, options_json),
        )
        row = self.conn.execute(
            """
            SELECT id FROM batches
            WHERE name = ? AND planet_pack = ?
              AND template_json IS ? AND generation_options_json = ?
            """,
            (name, planet_pack, template_json, options_json),
        ).fetchone()
        self.conn.commit()
        return row["id"]

    def upsert_sequence(self, planet_pack, short_name, bodies):
        bodies = list(bodies)
        body_path = ">".join(bodies)
        start_body = bodies[0]
        target_body = bodies[-1]
        flyby_count = max(0, len(bodies) - 2)
        self.conn.execute(
            """
            INSERT OR IGNORE INTO sequences(
                planet_pack, short_name, body_path, bodies_json, start_body, target_body, flyby_count
            ) VALUES (?, ?, ?, ?, ?, ?, ?)
            """,
            (
                planet_pack,
                short_name,
                body_path,
                _json_dumps(bodies),
                start_body,
                target_body,
                flyby_count,
            ),
        )
        seq = self.conn.execute(
            "SELECT id FROM sequences WHERE planet_pack = ? AND body_path = ?",
            (planet_pack, body_path),
        ).fetchone()
        sequence_id = seq["id"]
        self.conn.execute("DELETE FROM sequence_bodies WHERE sequence_id = ?", (sequence_id,))
        for position, body_name in enumerate(bodies):
            if position == 0:
                role = "departure"
            elif position == len(bodies) - 1:
                role = "arrival"
            else:
                role = "flyby"
            self.conn.execute(
                """
                INSERT INTO sequence_bodies(sequence_id, position, body_name, role)
                VALUES (?, ?, ?, ?)
                """,
                (sequence_id, position, body_name, role),
            )
        self.conn.commit()
        return sequence_id

    def _job_params_from_row(self, planet_pack, seq_short_name, t0_lb, tof_lb, row):
        bodies = list(row["mga_seq_fullname"])
        t0 = list(row["mga_t0"])
        tof = list(row["mga_tof"])
        vinf = list(row["mga_vinf"])
        batch_t0_bin = float(row["batch_t0_bin"])
        batch_tof_bin = float(row["batch_tof_bin"])
        return {
            "planet_pack": planet_pack,
            "sequence": bodies,
            "sequence_short_name": seq_short_name,
            "t0_min": float(t0_lb),
            "t0_max": float(t0_lb) + batch_t0_bin,
            "tof_min": float(tof_lb),
            "tof_max": float(tof_lb) + batch_tof_bin,
            "mga_t0": t0,
            "mga_tof": tof,
            "vinf": vinf,
            "tof_encoding": row["mga_tof_encoding"],
            "opt_level": row["batch_opt_level"],
            "opt_injection": row.get("batch_opt_insertion", ""),
            "ejection_altitude": float(row["mga_alt_start"]),
            "injection_altitude": float(row["mga_alt_target"]),
            "rp_target": float(row["mga_rp_target"]),
            "e_target": float(row["mga_e_target"]),
            "add_vinf_dep": bool(row["mga_add_vinf_dep"]),
            "add_vinf_arr": bool(row["mga_add_vinf_arr"]),
            "orbit_insertion": bool(row["mga_orbit_insertion"]),
            "multi_objective": bool(row["mga_multi_objective"]),
            "lambert_max_revs": int(row.get("mga_lambert_max_revs", 0)),
            "optimizer_topology": row.get("optimizer_topology", "ring"),
            "optimizer_seed": row.get("optimizer_seed", None),
        }

    def upsert_job(self, params, batch_id=None, status="TODO", versions=None, result=None):
        planet_pack = params["planet_pack"]
        seq_short_name = params["sequence_short_name"]
        bodies = list(params["sequence"])
        sequence_id = self.upsert_sequence(planet_pack, seq_short_name, bodies)
        job_hash = _param_hash(params)
        now = _utc_now()
        values = (
            sequence_id,
            planet_pack,
            params["t0_min"],
            params["t0_max"],
            params["tof_min"],
            params["tof_max"],
            float(params["vinf"][0]),
            float(params["vinf"][1]),
            params["tof_encoding"],
            params.get("tof_profile"),
            _json_dumps(params["leg_tof_bounds"]) if params.get("leg_tof_bounds") is not None else None,
            params["opt_level"],
            params["opt_injection"],
            params.get("arrival_mode", "legacy"),
            params["ejection_altitude"],
            params["injection_altitude"],
            params["rp_target"],
            params["e_target"],
            int(params["add_vinf_dep"]),
            int(params["add_vinf_arr"]),
            int(params["orbit_insertion"]),
            int(params["multi_objective"]),
            params["lambert_max_revs"],
            int(params["sade_gen"]),
            int(params["n_island"]),
            int(params["island_pop"]),
            int(params["n_evo_steps"]),
            params.get("optimizer_topology", "ring"),
            params.get("optimizer_seed", None),
            status,
            job_hash,
            now,
            now,
        )
        self.conn.execute(
            """
            INSERT INTO jobs(
                sequence_id, planet_pack, t0_min, t0_max, tof_min, tof_max,
                vinf_min, vinf_max, tof_encoding, tof_profile, leg_tof_bounds_json,
                opt_level, opt_injection, arrival_mode,
                ejection_altitude, injection_altitude, rp_target, e_target,
                add_vinf_dep, add_vinf_arr, orbit_insertion, multi_objective,
                lambert_max_revs, sade_gen, n_island, island_pop, n_evo_steps,
                optimizer_topology, optimizer_seed, status, param_hash, created_at, updated_at
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ON CONFLICT(param_hash) DO UPDATE SET
                status = CASE
                    WHEN jobs.status = 'DONE' AND excluded.status = 'TODO' THEN jobs.status
                    ELSE excluded.status
                END,
                updated_at = excluded.updated_at
            """,
            values,
        )
        job_id = self.conn.execute(
            "SELECT id FROM jobs WHERE param_hash = ?", (job_hash,)
        ).fetchone()["id"]
        job_id = int(job_id)
        if batch_id is not None:
            self.conn.execute(
                "INSERT OR IGNORE INTO batch_jobs(batch_id, job_id) VALUES (?, ?)",
                (batch_id, job_id),
            )
        if status == "DONE" and result is not None:
            self.upsert_result(
                job_id,
                result["result_DV"],
                result["result_t0"],
                result["result_tof"],
                result["result_ej_vinf"],
                result["gene"],
                versions=versions,
            )
        self.conn.commit()
        return job_id

    def upsert_job_from_row(self, planet_pack, index, row, batch_id=None, versions=None):
        seq_short_name, t0_lb, tof_lb = index
        params = self._job_params_from_row(planet_pack, seq_short_name, t0_lb, tof_lb, row)
        params.update({
            "sade_gen": int(row["job_sade_gen"]),
            "n_island": int(row["job_n_island"]),
            "island_pop": int(row["job_island_pop"]),
            "n_evo_steps": int(row["job_n_evo_steps"]),
            "optimizer_topology": row.get("optimizer_topology", "ring"),
            "optimizer_seed": row.get("optimizer_seed", None),
        })
        result = None
        if row["job_status"] == "DONE":
            result = {
                "result_DV": row["result_DV"],
                "result_t0": row["result_t0"],
                "result_tof": row["result_tof"],
                "result_ej_vinf": row["result_ej_vinf"],
                "gene": row["gene"],
            }
        return self.upsert_job(
            params,
            batch_id=batch_id,
            status=row["job_status"],
            versions=versions,
            result=result,
        )

    def update_job_status(self, job_id, status):
        if status in ("TODO", "DONE", "FAILED"):
            self.conn.execute(
                """
                UPDATE jobs
                SET status = ?, claimed_at = NULL, claim_expires_at = NULL,
                    worker_id = NULL, updated_at = ?
                WHERE id = ?
                """,
                (status, _utc_now(), int(job_id)),
            )
        else:
            self.conn.execute(
                "UPDATE jobs SET status = ?, updated_at = ? WHERE id = ?",
                (status, _utc_now(), int(job_id)),
            )
        self.conn.commit()

    def start_run(self, job_id, versions=None, optimizer_metadata=None):
        versions = versions or {}
        optimizer_metadata = optimizer_metadata or {}
        now = _utc_now()
        funnel_config = optimizer_metadata.get("funnel_config")
        funnel_config_json = (
            _json_dumps(funnel_config) if funnel_config is not None else None
        )
        funnel_config_hash = (
            hashlib.sha256(funnel_config_json.encode("utf-8")).hexdigest()
            if funnel_config_json is not None else None
        )
        cur = self.conn.execute(
            """
            INSERT INTO runs(
                job_id, pykep_version, pygmo_version, wayfinder_version,
                status, started_at, optimizer_topology, optimizer_seed,
                effective_optimizer_seed, requested_n_island, actual_n_island,
                island_pop, sade_gen, n_evo_steps, adaptive_stop_json,
                optimizer_strategy, funnel_config_json, funnel_config_hash,
                code_revision, planet_pack_hash
            ) VALUES (?, ?, ?, ?, 'RUNNING', ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                int(job_id),
                _sqlite_text(versions.get("pykep")),
                _sqlite_text(versions.get("pygmo")),
                _sqlite_text(versions.get("wayfinder")),
                now,
                _sqlite_text(optimizer_metadata.get("optimizer_topology")),
                optimizer_metadata.get("optimizer_seed"),
                optimizer_metadata.get("effective_optimizer_seed"),
                optimizer_metadata.get("requested_n_island"),
                optimizer_metadata.get("actual_n_island"),
                optimizer_metadata.get("island_pop"),
                optimizer_metadata.get("sade_gen"),
                optimizer_metadata.get("n_evo_steps"),
                _json_dumps(optimizer_metadata.get("adaptive_stop"))
                if optimizer_metadata.get("adaptive_stop") else None,
                _sqlite_text(optimizer_metadata.get("optimizer_strategy")),
                funnel_config_json,
                funnel_config_hash,
                _sqlite_text(optimizer_metadata.get("code_revision")),
                _sqlite_text(optimizer_metadata.get("planet_pack_hash")),
            ),
        )
        self.conn.commit()
        return int(cur.lastrowid)

    def record_optimizer_stage(self, run_id, summary):
        """Persist one funnel-stage summary for runtime and convergence analysis."""
        self.conn.execute(
            """
            INSERT OR REPLACE INTO optimizer_stages(
                run_id, stage_index, stage_name, n_island, island_pop,
                sade_gen, planned_evo_steps, actual_evo_steps,
                ejection_model, initialization, topology_name, migration_rate,
                exact_archive_size, adaptive_stop_json,
                algorithms_json, topology_vertices, topology_edges,
                migrants_published, migration_islands_active,
                migrations_accepted, best_fitness,
                average_fitness, runtime_seconds, stop_reason
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                int(run_id), int(summary["stage_index"]), summary["stage_name"],
                int(summary["n_island"]), int(summary["island_pop"]),
                int(summary["sade_gen"]), int(summary["planned_evo_steps"]),
                int(summary["actual_evo_steps"]), summary["ejection_model"],
                summary.get("initialization"),
                summary.get("topology"), summary.get("migration_rate"),
                summary.get("exact_archive_size"),
                _json_dumps(summary.get("adaptive_stop"))
                if summary.get("adaptive_stop") else None,
                _json_dumps(summary["algorithms"]),
                summary.get("topology_vertices"), summary.get("topology_edges"),
                summary.get("migrants_published"),
                summary.get("migration_islands_active"),
                summary.get("migrations_accepted"),
                float(summary["best_fitness"]),
                float(summary["average_fitness"]), float(summary["runtime_seconds"]),
                summary["stop_reason"],
            ),
        )
        self.conn.commit()

    def optimizer_stages(self, run_id):
        rows = self.conn.execute(
            "SELECT * FROM optimizer_stages WHERE run_id = ? ORDER BY stage_index",
            (int(run_id),),
        ).fetchall()
        result = []
        for row in rows:
            item = dict(row)
            item["algorithms"] = _json_loads(item.pop("algorithms_json"), [])
            item["adaptive_stop"] = _json_loads(item.pop("adaptive_stop_json"), None)
            result.append(item)
        return result

    def finish_run(
        self, run_id, runtime_seconds=None, actual_n_evo_steps=None,
        stop_reason=None,
    ):
        self.conn.execute(
            """
            UPDATE runs
            SET status = 'DONE', ended_at = ?, runtime_seconds = ?,
                actual_n_evo_steps = ?, stop_reason = ?, error = NULL
            WHERE id = ?
            """,
            (_utc_now(), runtime_seconds, actual_n_evo_steps, stop_reason, int(run_id)),
        )
        self.conn.commit()

    def complete_claimed_run(
        self, job_id, run_id, worker_id, claimed_at, objective_dv,
        result_t0, result_tof, ejection_vinf, gene, arrival_vinf=None,
        runtime_seconds=None, actual_n_evo_steps=None, stop_reason=None,
    ):
        """Atomically publish a result only for the current live claim."""
        job_id = int(job_id)
        run_id = int(run_id)
        now = _utc_now()
        try:
            self.conn.execute("BEGIN IMMEDIATE")
            claim = self.conn.execute(
                """
                SELECT 1 FROM jobs
                WHERE id = ? AND status = 'RUNNING' AND worker_id = ?
                  AND claimed_at = ? AND claim_expires_at > ?
                """,
                (job_id, str(worker_id), str(claimed_at), now),
            ).fetchone()
            if claim is None:
                raise RuntimeError(
                    "Job claim is no longer active; refusing stale result"
                )
            cur = self.conn.execute(
                """
                UPDATE runs
                SET status = 'DONE', ended_at = ?, runtime_seconds = ?,
                    actual_n_evo_steps = ?, stop_reason = ?, error = NULL
                WHERE id = ? AND job_id = ? AND status = 'RUNNING'
                """,
                (
                    now, runtime_seconds, actual_n_evo_steps, stop_reason,
                    run_id, job_id,
                ),
            )
            if cur.rowcount != 1:
                raise RuntimeError(
                    "Optimizer run is no longer active; refusing result"
                )
            self.conn.execute("DELETE FROM results WHERE run_id = ?", (run_id,))
            self.conn.execute(
                """
                INSERT INTO results(
                    run_id, objective_dv, result_t0, result_tof,
                    ejection_vinf, arrival_vinf
                ) VALUES (?, ?, ?, ?, ?, ?)
                """,
                (
                    run_id, float(objective_dv), float(result_t0),
                    float(result_tof), float(ejection_vinf),
                    None if arrival_vinf is None else float(arrival_vinf),
                ),
            )
            self.conn.execute("DELETE FROM genes WHERE run_id = ?", (run_id,))
            self.conn.execute(
                "INSERT INTO genes(run_id, gene_json) VALUES (?, ?)",
                (run_id, _json_dumps(list(gene))),
            )
            cur = self.conn.execute(
                """
                UPDATE jobs
                SET status = 'DONE', claimed_at = NULL,
                    claim_expires_at = NULL, worker_id = NULL, updated_at = ?
                WHERE id = ? AND status = 'RUNNING' AND worker_id = ?
                  AND claimed_at = ? AND claim_expires_at > ?
                """,
                (now, job_id, str(worker_id), str(claimed_at), now),
            )
            if cur.rowcount != 1:
                raise RuntimeError("Job claim changed while publishing result")
            self.conn.commit()
        except Exception:
            self.conn.rollback()
            raise

    def fail_claimed_run(
        self, job_id, run_id, worker_id, claimed_at, error,
        runtime_seconds=None,
    ):
        """Fail a run only if its worker still owns the same claim."""
        now = _utc_now()
        try:
            self.conn.execute("BEGIN IMMEDIATE")
            cur = self.conn.execute(
                """
                UPDATE jobs
                SET status = 'FAILED', claimed_at = NULL,
                    claim_expires_at = NULL, worker_id = NULL, updated_at = ?
                WHERE id = ? AND status = 'RUNNING' AND worker_id = ?
                  AND claimed_at = ?
                """,
                (now, int(job_id), str(worker_id), str(claimed_at)),
            )
            if cur.rowcount == 1:
                self.conn.execute(
                    """
                    UPDATE runs
                    SET status = 'FAILED', ended_at = ?, runtime_seconds = ?,
                        error = ?
                    WHERE id = ? AND job_id = ? AND status = 'RUNNING'
                    """,
                    (
                        now, runtime_seconds, str(error), int(run_id),
                        int(job_id),
                    ),
                )
            self.conn.commit()
            return cur.rowcount == 1
        except Exception:
            self.conn.rollback()
            raise

    def fail_run(self, run_id, error, runtime_seconds=None):
        self.conn.execute(
            """
            UPDATE runs
            SET status = 'FAILED', ended_at = ?, runtime_seconds = ?, error = ?
            WHERE id = ?
            """,
            (_utc_now(), runtime_seconds, str(error), int(run_id)),
        )
        self.conn.commit()

    def record_optimizer_snapshot(
        self, run_id, step, champion_fitness, champion_genes, populations=None,
        telemetry=None,
    ):
        telemetry = telemetry or {}
        champion_fitness = _jsonable(champion_fitness)
        champion_genes = _jsonable(champion_genes)
        fitness_values = [float(value[0] if isinstance(value, list) else value) for value in champion_fitness]
        best_fitness = min(fitness_values)
        population_fitness = []
        for population in populations or []:
            for value in _jsonable(population.get_f()):
                population_fitness.append(
                    float(value[0] if isinstance(value, list) else value)
                )
        average_fitness = (
            sum(population_fitness) / len(population_fitness)
            if population_fitness else None
        )
        self.conn.execute(
            """
            INSERT INTO optimizer_snapshots(
                run_id, step, best_fitness, average_fitness,
                champion_fitness_json, champion_gene_json, island_count,
                topology_vertices, topology_edges, migrants_published,
                migration_islands_active, migrations_accepted, created_at
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ON CONFLICT(run_id, step) DO UPDATE SET
                best_fitness = excluded.best_fitness,
                average_fitness = excluded.average_fitness,
                champion_fitness_json = excluded.champion_fitness_json,
                champion_gene_json = excluded.champion_gene_json,
                island_count = excluded.island_count,
                topology_vertices = excluded.topology_vertices,
                topology_edges = excluded.topology_edges,
                migrants_published = excluded.migrants_published,
                migration_islands_active = excluded.migration_islands_active,
                migrations_accepted = excluded.migrations_accepted,
                created_at = excluded.created_at
            """,
            (
                int(run_id),
                int(step),
                best_fitness,
                average_fitness,
                _json_dumps(champion_fitness),
                _json_dumps(champion_genes),
                telemetry.get("island_count"),
                telemetry.get("topology_vertices"),
                telemetry.get("topology_edges"),
                telemetry.get("migrants_published"),
                telemetry.get("migration_islands_active"),
                telemetry.get("migrations_accepted"),
                _utc_now(),
            ),
        )
        self.conn.commit()

    def upsert_result(self, job_id, objective_dv, result_t0, result_tof, ejection_vinf, gene, versions=None, run_id=None, arrival_vinf=None):
        job_id = int(job_id)
        versions = versions or {}
        now = _utc_now()
        if run_id is None:
            run = self.conn.execute(
                """
                SELECT id FROM runs
                WHERE job_id = ? AND status = 'DONE'
                ORDER BY id DESC LIMIT 1
                """,
                (job_id,),
            ).fetchone()
        else:
            run = self.conn.execute(
                "SELECT id FROM runs WHERE id = ? AND job_id = ?",
                (int(run_id), job_id),
            ).fetchone()
        if run is None and run_id is None:
            cur = self.conn.execute(
                """
                INSERT INTO runs(
                    job_id, pykep_version, pygmo_version, wayfinder_version,
                    status, started_at, ended_at
                ) VALUES (?, ?, ?, ?, 'DONE', ?, ?)
                """,
                (
                    job_id,
                    _sqlite_text(versions.get("pykep")),
                    _sqlite_text(versions.get("pygmo")),
                    _sqlite_text(versions.get("wayfinder")),
                    now,
                    now,
                ),
            )
            run_id = int(cur.lastrowid)
        elif run is None:
            raise ValueError("Unknown run_id for job: " + str(run_id))
        else:
            run_id = int(run["id"])
            self.conn.execute(
                """
                UPDATE runs
                SET status = 'DONE', ended_at = ?, error = NULL
                WHERE id = ?
                """,
                (now, run_id),
            )
        self.conn.execute("DELETE FROM results WHERE run_id = ?", (run_id,))
        self.conn.execute(
            """
            INSERT INTO results(
                run_id, objective_dv, result_t0, result_tof, ejection_vinf, arrival_vinf
            ) VALUES (?, ?, ?, ?, ?, ?)
            """,
            (
                run_id,
                float(objective_dv),
                float(result_t0),
                float(result_tof),
                float(ejection_vinf),
                None if arrival_vinf is None else float(arrival_vinf),
            ),
        )
        self.conn.execute("DELETE FROM genes WHERE run_id = ?", (run_id,))
        self.conn.execute(
            "INSERT INTO genes(run_id, gene_json) VALUES (?, ?)",
            (run_id, _json_dumps(list(gene))),
        )
        self.conn.commit()
        return run_id

    def upsert_result_from_row(self, job_id, row, versions=None, run_id=None):
        return self.upsert_result(
            job_id,
            row["result_DV"],
            row["result_t0"],
            row["result_tof"],
            row["result_ej_vinf"],
            row["gene"],
            versions=versions,
            run_id=run_id,
        )

    def import_dataframe(
        self,
        dataframe,
        planet_pack,
        batch_name,
        template=None,
        generation_options=None,
        versions=None,
        purpose="production",
    ):
        batch_id = self.upsert_batch(
            batch_name,
            planet_pack,
            template=template,
            generation_options=generation_options,
            purpose=purpose,
        )
        for index, row in dataframe.iterrows():
            self.upsert_job_from_row(
                planet_pack, index, row, batch_id=batch_id, versions=versions
            )
        self.conn.commit()
        return batch_id

    def best_results(
        self,
        planet_pack,
        batch_name=None,
        start_body=None,
        target_body=None,
        sequence_short_name=None,
        t0_range=None,
        contains_flyby=None,
        limit=10,
        include_benchmarks=False,
    ):
        where = ["j.planet_pack = ?", "j.status = 'DONE'"]
        params = [planet_pack]
        joins = []
        if batch_name is not None:
            joins.extend([
                "JOIN batch_jobs bj ON bj.job_id = j.id",
                "JOIN batches b ON b.id = bj.batch_id",
            ])
            where.append("b.name = ?")
            params.append(batch_name)
        elif not include_benchmarks:
            where.append(
                """
                (
                    NOT EXISTS (
                        SELECT 1 FROM batch_jobs bj_filter
                        WHERE bj_filter.job_id = j.id
                    )
                    OR EXISTS (
                        SELECT 1 FROM batch_jobs bj_filter
                        JOIN batches b_filter ON b_filter.id = bj_filter.batch_id
                        WHERE bj_filter.job_id = j.id
                          AND b_filter.purpose <> 'benchmark'
                    )
                )
                """
            )
        if start_body is not None:
            where.append("s.start_body = ?")
            params.append(start_body)
        if target_body is not None:
            where.append("s.target_body = ?")
            params.append(target_body)
        if sequence_short_name is not None:
            where.append("s.short_name = ?")
            params.append(sequence_short_name)
        if t0_range is not None:
            where.append("j.t0_min >= ? AND j.t0_max <= ?")
            params.extend([float(t0_range[0]), float(t0_range[1])])
        if contains_flyby is not None:
            where.append(
                """
                EXISTS (
                    SELECT 1 FROM sequence_bodies sb
                    WHERE sb.sequence_id = s.id
                      AND sb.role = 'flyby'
                      AND sb.body_name = ?
                )
                """
            )
            params.append(contains_flyby)
        params.append(int(limit))
        rows = self.conn.execute(
            f"""
            SELECT
                j.id AS job_id,
                s.short_name AS sequence_short_name,
                s.body_path,
                s.bodies_json,
                s.start_body,
                s.target_body,
                j.t0_min,
                j.t0_max,
                j.tof_min,
                j.tof_max,
                j.vinf_min,
                j.vinf_max,
                j.tof_encoding,
                j.tof_profile,
                j.leg_tof_bounds_json,
                j.opt_level,
                j.opt_injection,
                j.arrival_mode,
                j.ejection_altitude,
                j.injection_altitude,
                j.rp_target,
                j.e_target,
                j.add_vinf_dep,
                j.add_vinf_arr,
                j.orbit_insertion,
                j.multi_objective,
                j.lambert_max_revs,
                ru.id AS run_id,
                ru.started_at,
                ru.ended_at,
                ru.optimizer_topology,
                ru.optimizer_seed,
                ru.requested_n_island,
                ru.actual_n_island,
                ru.island_pop AS run_island_pop,
                ru.sade_gen AS run_sade_gen,
                ru.n_evo_steps AS run_n_evo_steps,
                ru.actual_n_evo_steps,
                ru.adaptive_stop_json,
                ru.optimizer_strategy,
                ru.stop_reason,
                ru.runtime_seconds,
                r.objective_dv,
                r.result_t0,
                r.result_tof,
                r.ejection_vinf,
                r.arrival_vinf,
                g.gene_json
            FROM jobs j
            JOIN sequences s ON s.id = j.sequence_id
            {" ".join(joins)}
            JOIN runs ru ON ru.job_id = j.id AND ru.status = 'DONE'
            JOIN results r ON r.run_id = ru.id
            LEFT JOIN genes g ON g.run_id = ru.id
            WHERE {" AND ".join(where)}
            ORDER BY r.objective_dv ASC
            LIMIT ?
            """,
            params,
        ).fetchall()
        return [dict(row) for row in rows]

    def run_context(self, run_id):
        row = self.conn.execute(
            """
            SELECT
                ru.id AS run_id,
                ru.status AS run_status,
                j.id AS job_id,
                s.short_name AS sequence_short_name,
                s.body_path,
                s.bodies_json,
                s.start_body,
                s.target_body,
                j.t0_min,
                j.t0_max,
                j.tof_min,
                j.tof_max,
                j.vinf_min,
                j.vinf_max,
                j.tof_encoding,
                j.tof_profile,
                j.leg_tof_bounds_json,
                j.opt_level,
                j.opt_injection,
                j.arrival_mode,
                j.ejection_altitude,
                j.injection_altitude,
                j.rp_target,
                j.e_target,
                j.add_vinf_dep,
                j.add_vinf_arr,
                j.orbit_insertion,
                j.multi_objective,
                j.lambert_max_revs,
                j.sade_gen,
                j.n_island,
                j.island_pop,
                j.n_evo_steps,
                j.optimizer_topology,
                j.optimizer_seed,
                ru.optimizer_topology AS run_optimizer_topology,
                ru.optimizer_seed AS run_optimizer_seed,
                ru.effective_optimizer_seed,
                ru.requested_n_island,
                ru.actual_n_island,
                ru.island_pop AS run_island_pop,
                ru.sade_gen AS run_sade_gen,
                ru.n_evo_steps AS run_n_evo_steps,
                ru.actual_n_evo_steps,
                ru.adaptive_stop_json,
                ru.optimizer_strategy,
                ru.funnel_config_json,
                ru.funnel_config_hash,
                ru.code_revision,
                ru.planet_pack_hash,
                ru.stop_reason,
                ru.runtime_seconds,
                r.objective_dv,
                r.result_t0,
                r.result_tof,
                r.ejection_vinf,
                g.gene_json
            FROM runs ru
            JOIN jobs j ON j.id = ru.job_id
            JOIN sequences s ON s.id = j.sequence_id
            LEFT JOIN results r ON r.run_id = ru.id
            LEFT JOIN genes g ON g.run_id = ru.id
            WHERE ru.id = ?
            """,
            (int(run_id),),
        ).fetchone()
        if row is None:
            return None
        return dict(row)

    def optimizer_snapshot_points(self, run_id):
        rows = self.conn.execute(
            """
            SELECT
                run_id,
                step,
                best_fitness,
                average_fitness,
                champion_fitness_json,
                champion_gene_json,
                island_count,
                topology_vertices,
                topology_edges,
                migrants_published,
                migration_islands_active,
                migrations_accepted,
                created_at
            FROM optimizer_snapshots
            WHERE run_id = ?
            ORDER BY step ASC
            """,
            (int(run_id),),
        ).fetchall()
        points = []
        for row in rows:
            fitness_values = _json_loads(row["champion_fitness_json"], default=[])
            genes = _json_loads(row["champion_gene_json"], default=[])
            for island_index, gene in enumerate(genes):
                fitness = fitness_values[island_index]
                if isinstance(fitness, list):
                    fitness = fitness[0]
                points.append({
                    "run_id": int(row["run_id"]),
                    "step": int(row["step"]),
                    "island_index": island_index,
                    "best_fitness": float(row["best_fitness"]),
                    "average_fitness": (
                        None if row["average_fitness"] is None
                        else float(row["average_fitness"])
                    ),
                    "fitness": float(fitness),
                    "gene": gene,
                    "island_count": row["island_count"],
                    "topology_vertices": row["topology_vertices"],
                    "topology_edges": row["topology_edges"],
                    "migrants_published": row["migrants_published"],
                    "migration_islands_active": row["migration_islands_active"],
                    "migrations_accepted": row["migrations_accepted"],
                    "created_at": row["created_at"],
                })
        return points

    def optimizer_convergence(self, run_id):
        """Return one best/average population fitness sample per evolution step."""
        rows = self.conn.execute(
            """
            SELECT step, best_fitness, average_fitness, island_count,
                   topology_vertices, topology_edges, migrants_published,
                   migration_islands_active, migrations_accepted, created_at
            FROM optimizer_snapshots
            WHERE run_id = ?
            ORDER BY step ASC
            """,
            (int(run_id),),
        ).fetchall()
        return [dict(row) for row in rows]

    def optimizer_population_points(self, run_id, source=None):
        where = ["run_id = ?"]
        params = [int(run_id)]
        if source is not None:
            where.append("source = ?")
            params.append(source)
        rows = self.conn.execute(
            f"""
            SELECT
                run_id,
                step,
                island_index,
                individual_index,
                fitness,
                fitness_json,
                gene_json,
                source,
                created_at
            FROM optimizer_population_points
            WHERE {" AND ".join(where)}
            ORDER BY step ASC, island_index ASC, individual_index ASC
            """,
            params,
        ).fetchall()
        points = []
        for row in rows:
            points.append({
                "run_id": int(row["run_id"]),
                "step": int(row["step"]),
                "island_index": int(row["island_index"]),
                "individual_index": int(row["individual_index"]),
                "fitness": float(row["fitness"]),
                "fitness_vector": _json_loads(row["fitness_json"], default=[]),
                "gene": _json_loads(row["gene_json"], default=[]),
                "source": row["source"],
                "created_at": row["created_at"],
            })
        return points

    def result_rows(
        self,
        planet_pack,
        batch_name=None,
        start_body=None,
        target_body=None,
        sequence_short_name=None,
        t0_range=None,
        contains_flyby=None,
        include_benchmarks=False,
    ):
        where = ["j.planet_pack = ?", "j.status = 'DONE'"]
        params = [planet_pack]
        joins = []
        if batch_name is not None:
            joins.extend([
                "JOIN batch_jobs bj ON bj.job_id = j.id",
                "JOIN batches b ON b.id = bj.batch_id",
            ])
            where.append("b.name = ?")
            params.append(batch_name)
        elif not include_benchmarks:
            where.append(
                """
                (
                    NOT EXISTS (
                        SELECT 1 FROM batch_jobs bj_filter
                        WHERE bj_filter.job_id = j.id
                    )
                    OR EXISTS (
                        SELECT 1 FROM batch_jobs bj_filter
                        JOIN batches b_filter ON b_filter.id = bj_filter.batch_id
                        WHERE bj_filter.job_id = j.id
                          AND b_filter.purpose <> 'benchmark'
                    )
                )
                """
            )
        if start_body is not None:
            where.append("s.start_body = ?")
            params.append(start_body)
        if target_body is not None:
            where.append("s.target_body = ?")
            params.append(target_body)
        if sequence_short_name is not None:
            where.append("s.short_name = ?")
            params.append(sequence_short_name)
        if t0_range is not None:
            where.append("j.t0_min >= ? AND j.t0_max <= ?")
            params.extend([float(t0_range[0]), float(t0_range[1])])
        if contains_flyby is not None:
            where.append(
                """
                EXISTS (
                    SELECT 1 FROM sequence_bodies sb
                    WHERE sb.sequence_id = s.id
                      AND sb.role = 'flyby'
                      AND sb.body_name = ?
                )
                """
            )
            params.append(contains_flyby)
        rows = self.conn.execute(
            f"""
            SELECT
                j.id AS job_id,
                s.short_name AS sequence_short_name,
                s.body_path,
                s.start_body,
                s.target_body,
                j.t0_min,
                j.t0_max,
                j.tof_min,
                j.tof_max,
                j.opt_level,
                j.opt_injection,
                ru.optimizer_topology,
                ru.optimizer_seed,
                ru.requested_n_island,
                ru.actual_n_island,
                ru.runtime_seconds,
                r.objective_dv,
                r.result_t0,
                r.result_tof,
                r.ejection_vinf
            FROM jobs j
            JOIN sequences s ON s.id = j.sequence_id
            {" ".join(joins)}
            JOIN runs ru ON ru.job_id = j.id AND ru.status = 'DONE'
            JOIN results r ON r.run_id = ru.id
            WHERE {" AND ".join(where)}
            ORDER BY s.short_name ASC, j.t0_min ASC, r.objective_dv ASC
            """,
            params,
        ).fetchall()
        return [dict(row) for row in rows]

    def benchmark_results(self, planet_pack, benchmark_name=None):
        where = ["j.planet_pack = ?", "b.purpose = 'benchmark'", "ru.status = 'DONE'"]
        params = [planet_pack]
        if benchmark_name is not None:
            where.append("b.name LIKE ?")
            params.append(str(benchmark_name) + "%")
        rows = self.conn.execute(
            f"""
            SELECT
                b.name AS batch_name,
                s.short_name AS sequence_short_name,
                s.body_path,
                j.t0_min,
                j.t0_max,
                j.tof_min,
                j.tof_max,
                j.opt_level,
                j.opt_injection,
                ru.id AS run_id,
                ru.optimizer_topology,
                ru.optimizer_seed,
                ru.requested_n_island,
                ru.actual_n_island,
                ru.island_pop AS run_island_pop,
                ru.sade_gen AS run_sade_gen,
                ru.n_evo_steps AS run_n_evo_steps,
                ru.actual_n_evo_steps,
                ru.adaptive_stop_json,
                ru.optimizer_strategy,
                ru.stop_reason,
                ru.runtime_seconds,
                r.objective_dv,
                r.result_t0,
                r.result_tof,
                r.ejection_vinf
            FROM batches b
            JOIN batch_jobs bj ON bj.batch_id = b.id
            JOIN jobs j ON j.id = bj.job_id
            JOIN sequences s ON s.id = j.sequence_id
            JOIN runs ru ON ru.job_id = j.id
            JOIN results r ON r.run_id = ru.id
            WHERE {" AND ".join(where)}
            ORDER BY b.name ASC, r.objective_dv ASC
            """,
            params,
        ).fetchall()
        return [dict(row) for row in rows]

    def _select_jobs(self, where, params, joins=None, limit=None):
        joins = list(joins or ["JOIN sequences s ON s.id = j.sequence_id"])
        query_limit = ""
        query_params = list(params)
        if limit is not None:
            query_limit = "LIMIT ?"
            query_params.append(int(limit))
        rows = self.conn.execute(
            f"""
            SELECT DISTINCT
                j.id AS job_id,
                s.short_name AS sequence_short_name,
                s.body_path,
                s.bodies_json,
                j.t0_min,
                j.t0_max,
                j.tof_min,
                j.tof_max,
                j.vinf_min,
                j.vinf_max,
                j.tof_encoding,
                j.tof_profile,
                j.leg_tof_bounds_json,
                j.opt_level,
                j.opt_injection,
                j.arrival_mode,
                j.ejection_altitude,
                j.injection_altitude,
                j.rp_target,
                j.e_target,
                j.add_vinf_dep,
                j.add_vinf_arr,
                j.orbit_insertion,
                j.multi_objective,
                j.lambert_max_revs,
                j.sade_gen,
                j.n_island,
                j.island_pop,
                j.n_evo_steps,
                j.optimizer_topology,
                j.optimizer_seed,
                j.status,
                j.claimed_at,
                j.claim_expires_at,
                j.worker_id
            FROM jobs j
            {" ".join(joins)}
            WHERE {" AND ".join(where)}
            ORDER BY j.t0_min ASC, j.tof_min ASC, s.short_name ASC
            {query_limit}
            """,
            query_params,
        ).fetchall()
        return [dict(row) for row in rows]

    def pending_jobs(self, planet_pack, limit=10, batch_name=None):
        where = ["j.planet_pack = ?", "j.status = 'TODO'"]
        params = [planet_pack]
        joins = [
            "JOIN sequences s ON s.id = j.sequence_id",
        ]
        if batch_name is not None:
            joins.extend([
                "JOIN batch_jobs bj ON bj.job_id = j.id",
                "JOIN batches b ON b.id = bj.batch_id",
            ])
            where.append("b.name = ?")
            params.append(batch_name)
        return self._select_jobs(where, params, joins=joins, limit=limit)

    def claim_pending_jobs(
        self, planet_pack, limit=10, batch_name=None, worker_id=None,
        lease_seconds=86400,
    ):
        """Atomically recover expired leases and claim pending jobs."""
        limit = int(limit)
        if limit <= 0:
            return []
        worker_id = str(worker_id or _default_worker_id())
        now = _utc_now()
        expires_at = _utc_after(lease_seconds)
        try:
            self.conn.execute("BEGIN IMMEDIATE")
            self.conn.execute(
                """
                UPDATE runs
                SET status = 'FAILED', ended_at = ?,
                    stop_reason = 'claim_expired',
                    error = COALESCE(
                        error, 'Worker claim expired before run completion'
                    )
                WHERE status = 'RUNNING'
                  AND job_id IN (
                      SELECT id FROM jobs
                      WHERE status = 'RUNNING'
                        AND claim_expires_at IS NOT NULL
                        AND claim_expires_at <= ?
                  )
                """,
                (now, now),
            )
            self.conn.execute(
                """
                UPDATE jobs
                SET status = 'TODO', claimed_at = NULL,
                    claim_expires_at = NULL, worker_id = NULL, updated_at = ?
                WHERE status = 'RUNNING'
                  AND claim_expires_at IS NOT NULL
                  AND claim_expires_at <= ?
                """,
                (now, now),
            )
            where = ["j.planet_pack = ?", "j.status = 'TODO'"]
            params = [planet_pack]
            joins = ["JOIN sequences s ON s.id = j.sequence_id"]
            if batch_name is not None:
                joins.extend([
                    "JOIN batch_jobs bj ON bj.job_id = j.id",
                    "JOIN batches b ON b.id = bj.batch_id",
                ])
                where.append("b.name = ?")
                params.append(batch_name)
            ids = [
                int(row["job_id"])
                for row in self._select_jobs(
                    where, params, joins=joins, limit=limit
                )
            ]
            if ids:
                placeholders = ",".join("?" for _ in ids)
                self.conn.execute(
                    f"""
                    UPDATE jobs
                    SET status = 'RUNNING', claimed_at = ?,
                        claim_expires_at = ?, worker_id = ?, updated_at = ?
                    WHERE status = 'TODO' AND id IN ({placeholders})
                    """,
                    [now, expires_at, worker_id, now] + ids,
                )
            self.conn.commit()
        except Exception:
            self.conn.rollback()
            raise
        if not ids:
            return []
        placeholders = ",".join("?" for _ in ids)
        return self._select_jobs(
            [
                f"j.id IN ({placeholders})",
                "j.status = 'RUNNING'",
                "j.worker_id = ?",
                "j.claimed_at = ?",
            ],
            ids + [worker_id, now],
            limit=None,
        )

    def renew_job_claim(
        self, job_id, worker_id, lease_seconds=86400, claimed_at=None,
    ):
        """Extend a live lease owned by ``worker_id``."""
        now = _utc_now()
        claim_clause = ""
        params = [
            _utc_after(lease_seconds), now, int(job_id), str(worker_id), now,
        ]
        if claimed_at is not None:
            claim_clause = " AND claimed_at = ?"
            params.append(str(claimed_at))
        cur = self.conn.execute(
            f"""
            UPDATE jobs
            SET claim_expires_at = ?, updated_at = ?
            WHERE id = ? AND status = 'RUNNING' AND worker_id = ?
              AND claim_expires_at > ?
              {claim_clause}
            """,
            params,
        )
        self.conn.commit()
        return cur.rowcount == 1

    def release_job_claim(self, job_id, worker_id, claimed_at=None):
        """Return an unstarted claimed job to the pending queue."""
        claim_clause = ""
        params = [_utc_now(), int(job_id), str(worker_id)]
        if claimed_at is not None:
            claim_clause = " AND claimed_at = ?"
            params.append(str(claimed_at))
        cur = self.conn.execute(
            f"""
            UPDATE jobs
            SET status = 'TODO', claimed_at = NULL,
                claim_expires_at = NULL, worker_id = NULL, updated_at = ?
            WHERE id = ? AND status = 'RUNNING' AND worker_id = ?
              {claim_clause}
            """,
            params,
        )
        self.conn.commit()
        return cur.rowcount == 1

    def count_rows(self, table):
        return self.conn.execute(f"SELECT COUNT(*) AS n FROM {table}").fetchone()["n"]
