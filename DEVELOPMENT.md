# Wayfinder Development Status

Last updated: 2026-06-18

This document is the handoff point for continuing Wayfinder development from a
new machine or a new coding session. Read it together with `README.md` and the
release changelog.

## Current release work

- Target release: `v1.6.0`
- Development branch: `codex/release-v1.6.0`
- Draft pull request: https://github.com/Muetdhiver-lab/Wayfinder/pull/1
- Last release commit at the time of writing: `90690af`
- Public `main` still represents the historical `v1.414` code until the pull
  request is merged.
- The directory name `KSP_Wayfinder v1.414` is historical. The authoritative
  version is stored in the root `VERSION` file.

## Supported runtime

- Python environment used during development:
  `C:\MyPrograms\MiniConda\envs\Wayfinder-pykep3`
- pykep: `3.0.0`
- pygmo: `2.19.8`
- Wayfinder is now pykep3-only. Reintroducing the old pykep2 wrapper or monkey
  patch is covered by a regression guard.

The exact environment path is machine-specific. A reproducible environment
file still needs to be added before `v1.6.0` is finalized.

## Architecture

Core modules live under `KSP_Wayfinder v1.414/WayfinderCore`:

- `_Wayfinder.py`: public orchestration API, batch generation, optimization
  workflow, queries, and plotting entry points.
- `_SQLiteStore.py`: versioned SQLite schema, migrations, jobs, batches, runs,
  results, genes, benchmarks, and porkchop samples. Current schema version: 7.
- `_Optimization.py`: optimization fitness decoration and related helpers.
- `_Trajectory.py`: trajectory decoding, metrics, plotting support, and TransX
  flight-plan output. It replaces the former `_Kraken_Patch.py` monkey patch.
- `planet_packs/`: Vanilla and JNSQ definitions. Keep the term `planet_packs` in
  APIs and documentation because that is the KSP community terminology.

SQLite is the only active datastore. Excel files and obsolete scripts are kept
under `legacy/` for historical reference only. Pandas remains acceptable at the
analysis and plotting boundary, but is no longer optimizer or datastore state.

## Important design decisions

- A single SQLite database may contain multiple planet packs. Every relevant
  row is scoped by `planet_pack`.
- Canonical parameter hashes deduplicate jobs while allowing different binning
  strategies to coexist.
- Benchmark batches use `purpose="benchmark"` and are excluded from normal
  result queries unless `include_benchmarks=True` is requested.
- Optimizer runs retain snapshots and final populations, not only the winning
  gene. These genes support local refinement and porkchop analysis.
- Alpha genes can be converted to direct encoding. A direct re-optimization of
  the difficult Vanilla KEKKJ solution preserved and slightly improved its
  objective (`2876.8206` to `2876.8034` m/s), validating the conversion path.
- Injection/capture metrics are separated from stored trajectory results so
  alternative injection assumptions can be evaluated without re-optimizing.
- Non-TransX diagnostics use Python logging. TransX output intentionally still
  uses direct console printing.
- Automatic worker selection uses the detected CPU count minus two reserved
  cores, while allowing an explicit island count.

## Tests

Run from `KSP_Wayfinder v1.414`:

```powershell
& 'C:\MyPrograms\MiniConda\envs\Wayfinder-pykep3\python.exe' `
  -m unittest discover -s Tests -p 'test_*.py' -v
```

Current result: 27 tests passing.

The test suite covers SQLite migrations and queries, benchmark filtering,
optimizer metadata and snapshots, plotting data selection, known Vanilla and
JNSQ trajectory decoding, time scaling, ejection geometry, alpha-to-direct
conversion, topology construction, CPU/island selection, and the pykep2 monkey
patch guard.

The eight committed SQLite files in `Tests/` are intentional reference and
benchmark snapshots. Treat them as immutable fixtures where possible; repeated
binary rewrites would unnecessarily increase Git history.

## Known local warning

Heyoka may warn that it cannot enable WAL mode for:

`C:\Users\vfave\AppData\Local\heyoka\cache\cache.db`

The cache file exists, but the local environment may prevent creation of its
WAL sidecar. This has no observed numerical or test impact; it only disables
effective use of the on-disk compilation cache.

## Optimization benchmark status

The first topology experiments compare `fully_connected`, `ring`, and
`unconnected` archipelagos on:

- JNSQ KEEMo in the known useful launch window.
- Vanilla KEKKJ using a difficult but direct-validated `v_inf` minimum.

Initial single-trial results exist in
`Tests/topology_benchmark_bunch1.sqlite`. Additional KEKKJ searches and direct
validation are stored in the other committed benchmark databases.

Do not interpret one optimizer seed as a deterministic replay: pygmo
archipelagos and multiprocessing did not reproduce identical populations from
the same global seed. Repeated runs should be treated as statistical trials.

## Next development steps

1. Add a reproducible pykep3 environment specification and verify setup from a
   clean machine.
2. Complete the topology benchmark with several trials for KEEMo and KEKKJ,
   then compare success rate, best objective, dispersion, and runtime.
3. Decide whether optimizer `seed` should be renamed to `trial_seed`, or seed
   individual islands explicitly if pygmo provides a robust path.
4. Review the optimization funnel after topology results: island/population
   balance, algorithm mix, migration policy, and phase transition criteria.
5. Continue the planned re-optimization and locally refined porkchop workflow.
6. Design the future GUI only after the optimizer and persistence APIs are
   stable. NiceGUI is currently the leading lightweight option; Qt was judged
   unnecessarily heavy for this project.

## Starting on another machine

Before the pull request is merged:

```powershell
git clone https://github.com/Muetdhiver-lab/Wayfinder.git
cd Wayfinder
git switch codex/release-v1.6.0
```

After merge, use `main` and the eventual `v1.6.0` tag instead.
