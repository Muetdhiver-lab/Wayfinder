# Wayfinder v1.6.0

Wayfinder is a Multiple Gravity Assist search tool for Kerbal Space Program.
It uses ESA's pykep and pygmo libraries to search and refine trajectories for
Vanilla KSP and JNSQ planet packs.

Version 1.6.0 targets pykep 3 and uses SQLite as its only datastore. Jobs,
optimizer runs, results, genes, benchmark metadata, and sampled porkchop data
can be queried across batches without depending on a particular binning
strategy.

## Repository layout

- `KSP_Wayfinder v1.414/Examples`: runnable search examples.
- `KSP_Wayfinder v1.414/WayfinderCore`: optimization, trajectory, plotting,
  SQLite, and planet-pack modules.
- `KSP_Wayfinder v1.414/Tests`: regression tests, benchmark scripts, and
  reference SQLite snapshots.
- `KSP_Wayfinder v1.414/legacy`: former Excel datastores and pykep2-era
  scripts retained for historical reference.

The directory name is historical; the authoritative release number is stored
in `VERSION`.

See [`KSP_Wayfinder v1.414/DOC/README.md`](KSP_Wayfinder%20v1.414/DOC/README.md)
for the usage guide and
[`KSP_Wayfinder v1.414/DOC/CHANGE LOG.md`](KSP_Wayfinder%20v1.414/DOC/CHANGE%20LOG.md)
for release history.

## References

- [pykep](https://esa.github.io/pykep/)
- [pygmo](https://esa.github.io/pygmo2/)
- [JNSQ launch window planner](https://github.com/LouisB3/ksp-lwp-jnsq)
- [Original KSP transfer planner](https://github.com/alexmoon/ksp)
