# Wayfinder development handoff

Last updated: 2026-06-22

This file is the short operational handoff. The detailed workflow, class and
database diagrams, benchmark history and GUI-readiness critique are in
[`KSP_Wayfinder v1.6.0/DOC/RELEASE_ARCHITECTURE_REVIEW_1.6.0.md`](KSP_Wayfinder%20v1.6.0/DOC/RELEASE_ARCHITECTURE_REVIEW_1.6.0.md).

## Release state

- Target release: `v1.6.0`.
- Development branch: `codex/release-v1.6.0`.
- Authoritative version: root `VERSION` and `WAYFINDER_VERSION` in
  `_Wayfinder.py`, both `1.6.0`.
- SQLite schema: **14**.
- Active datastore: SQLite only; Excel files are historical fixtures under
  `legacy/`.
- Core runtime: pykep 3 / pygmo 2, Python 3.10.

The core directory follows the authoritative release version:
`KSP_Wayfinder v1.6.0`.

## Reproducible environment

From the repository root:

```powershell
conda env create -f environment.yml
conda activate wayfinder
cd 'KSP_Wayfinder v1.6.0'
python -m pytest -q
```

The validated environment pins Python 3.10.20, pykep 3.0.0, pygmo 2.19.8,
NumPy 2.2.6, SciPy 1.15.2, NetworkX 3.4.2, pandas 2.3.3,
Matplotlib 3.10.9 and seaborn 0.13.2.

On Windows, activate the environment or use `conda run`. Invoking its
`python.exe` directly without `<env>\Library\bin` in `PATH` can produce delayed
native-DLL crashes during plotting even though imports succeed.

## Current architecture

Core modules under `KSP_Wayfinder v1.6.0/WayfinderCore`:

- `_Wayfinder.py`: public façade, job generation, Pygmo orchestration,
  refinement and plots;
- `_SQLiteStore.py`: schema v14, atomic job leases, runs, results, optimizer
  telemetry and porkchop samples;
- `_Optimization.py`: ToF encoding and fitness decoration;
- `_Trajectory.py`: canonical ejection calculation, decoding and TransX;
- `planet_packs/`: Vanilla and JNSQ definitions.

The public optimizer default for 1.6 remains deliberately
`optimizer_strategy="funnel"`, the corrected three-stage reference. The newer
`funnel_scout_archive`, `_32`, `_64` and `_128` strategies remain explicit
research options. The 128-island L0 variant is rejected; 32 islands is the next
candidate to repeat before any default changes.

## Job lifecycle and schema v14

`optimize_sqlite()` no longer reads naked `TODO` rows. It atomically claims
jobs using `BEGIN IMMEDIATE` and stores:

- `status = RUNNING`;
- `claimed_at`;
- `claim_expires_at`;
- `worker_id`.

Expired leases are returned to `TODO` on the next claim. Claims may be renewed,
and terminal statuses clear ownership fields. `optimize_sqlite()` claims only
the next job it is about to start, renews the active lease between Pygmo epochs,
and atomically publishes run, result, gene and final job status only when the
same `worker_id` and `claimed_at` still own a live claim. A stale worker is
therefore fenced out after recovery.

Stage telemetry now persists `topology_name`, `migration_rate` and
`exact_archive_size`, in addition to topology dimensions and migration counts.

## Optimization decisions

- New searches should prefer direct per-leg ToF encoding and bin only `T0`.
- The default topology is ring, with migration rate 2 in the funnel.
- Pygmo algorithms only: SADE, simulated annealing and
  `pg.nlopt("neldermead")`.
- Exact ejection chooses the cheaper direct-inclined or planar-plus-SOI-split
  strategy and is shared by fitness and decoding.
- The research scout/archive funnel preserves phase-diverse exact candidates
  across L1/L2 and uses adaptive exact L3 refinement.
- Pygmo multiprocessing migration is asynchronous; a repeated seed does not
  imply a bitwise deterministic archipelago history.

## Verification

Run from `KSP_Wayfinder v1.6.0`:

```powershell
python -m pytest -q
python -m compileall -q WayfinderCore Tests
```

Current suite: **52 tests** (14 datastore and 38 regression tests). It covers
v13-to-v14 migration, atomic claims, lease recovery and stale-worker fencing,
topology/migration,
funnels, exact archive, adaptive stopping, Vanilla/JNSQ decoding, exact
ejection consistency and SQLite analysis paths.

The standalone `Tests/run_*.py` scripts are benchmarks/smoke harnesses and are
not all part of pytest.

## Release 1.6 scope audit

The two tracked SQLite fixtures are intentionally part of the release. Both
were migrated through `SQLiteJobStore` to schema v14 and pass SQLite integrity
and foreign-key checks. The six PNG files under `DOC/assets` are the benchmark
figures referenced by the architecture review. Python/pytest caches and SQLite
WAL, SHM or journal files are generated artifacts and are not release inputs.

## Known environment warning

Heyoka may warn that it cannot write its on-disk cache under the local user
profile. This has no observed numerical effect but removes effective cache
reuse and can distort cold/warm runtime comparisons.

## Next work after the core release

1. Commit and tag the intentionally reviewed release scope.
2. Extract a transactional `SQLiteRepository` and observable
   `OptimizationService` behind the existing façade.
3. Add GUI-safe progress, cancellation and crash recovery.
4. Separate plot data models from Matplotlib rendering.
5. Repeat the 32-island L0 qualification before reconsidering the optimizer
   default.
6. Develop the planar/normal `v_inf` two-dimensional ejection surrogate.
