# Wayfinder v1.7.0

Wayfinder is a Multiple Gravity Assist search tool for Kerbal Space Program.
It uses ESA's pykep and pygmo libraries to search and refine trajectories for
Vanilla KSP and JNSQ planet packs.

Version 1.7.0 targets pykep 3 and uses SQLite as its only datastore. Jobs,
optimizer runs, results, genes, benchmark metadata, and sampled porkchop data
can be queried across batches without depending on a particular binning
strategy.

## Environment setup

Create the validated conda-forge environment from the repository root:

```powershell
conda env create -f environment.yml
conda activate wayfinder
cd 'KSP_Wayfinder v1.7.0'
python -m pytest -q
```

On Windows, prefer `conda activate wayfinder` or
`conda run -n wayfinder python ...`. Calling the environment's `python.exe`
directly without activation omits `<env>\Library\bin` from `PATH`. Imports may
still succeed, but native libraries such as Matplotlib can fail only when a
plot is first rendered. If direct invocation is unavoidable, prepend the
environment DLL directory for that process:

```powershell
$env:PATH="C:\path\to\env\Library\bin;$env:PATH"
C:\path\to\env\python.exe your_script.py
```

## Repository layout

- `KSP_Wayfinder v1.7.0/Examples`: runnable search examples.
- `KSP_Wayfinder v1.7.0/WayfinderCore`: optimization, trajectory, plotting,
  SQLite, and planet-pack modules.
- `KSP_Wayfinder v1.7.0/Tests`: regression tests, benchmark scripts, and
  reference SQLite snapshots.
- `KSP_Wayfinder v1.7.0/legacy`: former Excel datastores and pykep2-era
  scripts retained for historical reference.

The directory name and authoritative release number in `VERSION` are aligned.

See [`DEVELOPMENT.md`](DEVELOPMENT.md) before continuing active development
from another machine or coding session.

See [`KSP_Wayfinder v1.7.0/DOC/README.md`](KSP_Wayfinder%20v1.7.0/DOC/README.md)
for the usage guide and
[`KSP_Wayfinder v1.7.0/DOC/CHANGE LOG.md`](KSP_Wayfinder%20v1.7.0/DOC/CHANGE%20LOG.md)
for release history.

The pre-release workflow, architecture, class/DB diagrams, benchmark history
and GUI-readiness critique are consolidated in
[`KSP_Wayfinder v1.7.0/DOC/RELEASE_ARCHITECTURE_REVIEW_1.6.0.md`](KSP_Wayfinder%20v1.7.0/DOC/RELEASE_ARCHITECTURE_REVIEW_1.6.0.md).

## References

- [pykep](https://esa.github.io/pykep/)
- [pygmo](https://esa.github.io/pygmo2/)
- [JNSQ launch window planner](https://github.com/LouisB3/ksp-lwp-jnsq)
- [Original KSP transfer planner](https://github.com/alexmoon/ksp)
