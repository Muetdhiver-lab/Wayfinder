# Tisserand sequence scout

## Purpose

`SequenceScout` is an upstream search layer that proposes gravity-assist body
sequences before Wayfinder creates L0 jobs. It reduces the discrete sequence
space; it does not replace the continuous MGA-1DSM optimizer.

The implementation follows the standard preliminary-design use of Tisserand
graphs: constant-
\(v_\infty\) encounter contours connect heliocentric energy regions, while a
tree search enumerates encounter paths. This is consistent with the automatic
tree-search formulation described by de la Torre Sangrà et al. and with the
classical Tisserand parameter:

\[
T = \frac{R_p}{a}
    + 2\sqrt{\frac{a(1-e^2)}{R_p}}\cos i.
\]

References:

- [An Automatic Tree Search Algorithm for the Tisserand Graph](https://arxiv.org/abs/2006.07666)
- [ESA ICATT paper on Tisserand sequence selection](https://indico.esa.int/event/111/contributions/416/attachments/426/471/6th-ICATT-paper_02.pdf)

## Model contract

The current model is deliberately:

- circular and coplanar;
- patched-conic and unpowered at flybys;
- independent of epoch, time of flight and planetary phase;
- limited by each body's configured safe periapsis radius;
- bounded and deterministic through multi-resolution state pruning.

It can reject energy-incompatible chains and recover useful sequence families.
It cannot prove that a sequence is phaseable or that a low-DSM trajectory
exists. Every returned record is tagged
`tisserand_circular_coplanar_unphased` for that reason.

## Workflow

```mermaid
flowchart LR
    Goal["Start, target, planet pack"]
    TG["SequenceScout<br/>Tisserand tree"]
    L1["Lambert arc-1 filter<br/>real ephemerides"]
    L0["L0 recall"]
    Funnel["Pygmo funnel"]
    Goal --> TG --> L1 --> L0 --> Funnel
```

For each sampled departure \(v_\infty\) and direction, the scout constructs a
heliocentric orbit. If that orbit crosses a candidate body's circular reference
radius, it computes the incoming planet-relative velocity. The flyby may rotate
that vector by at most:

\[
\delta_{max}=2\arcsin\left(\frac{1}{1+r_p v_\infty^2/\mu_p}\right).
\]

The outgoing directions create the next set of Tisserand contours. States are
deduplicated in semi-major-axis, eccentricity, direction and \(v_\infty\) bins.
Each body sequence retains both its lowest proxy-cost states and a deterministic
spread across those bins, preventing narrow energy corridors from disappearing
through greedy pruning.

By default, candidate bodies are radially bounded between departure and target.
The nearest planet on the opposite side of the departure orbit is also included
because an initial inward assist can enable an outward mission, as in KEKKJ.

## API

```python
plans = Wayfinder(planet_pack="Vanilla")
candidates = plans.scout_sequences(
    start="Kerbin",
    target="Jool",
    limit=20,
)

for candidate in candidates:
    print(candidate.sequence, candidate.proxy_cost_mps)
```

`TisserandScoutConfig` exposes the maximum sequence length, departure and flyby
sampling, state-retention budget, encounter \(v_\infty\) ceiling, permitted
repeat visits and optional departure \(v_\infty\) bounds. It serializes to a
plain dictionary for future GUI and database integration.

Each candidate contains:

- the body sequence;
- departure and terminal \(v_\infty\);
- proxy cost and maximum turn utilization;
- terminal heliocentric periapsis/apoapsis;
- per-flyby incoming \(v_\infty\), incoming angle, selected/maximal turn,
  required/safe periapsis and Tisserand parameter;
- an explicit model-fidelity tag.

The proxy is only a stable ordering heuristic: departure \(v_\infty\), a small
flyby-count cost and a penalty for using the edge of the available turn. It is
not a predicted mission delta-v.

## Initial validation

With the default configuration:

| Pack and goal | Runtime | Candidates | Known sequence |
|---|---:|---:|---:|
| Vanilla Kerbin → Jool | about 1.2 s | 58 | KEKKJ rank 4 |
| JNSQ Kerbin → Moho | about 0.6 s | 20 | KEEMo rank 2 |

Unit tests also verify the equivalence between the geometric Tisserand formula
and \(T=3-(v_\infty/v_p)^2\) in the circular model, as well as the hyperbolic
turn/periapsis relation.

## Lambert arc-1 filter

`LambertArc1Filter` introduces real PyKEP ephemerides and the requested launch
window, but only for the first leg. For each Tisserand candidate it:

1. scan a coarse T0/TOF grid for the first encounter;
2. solve zero-revolution PyKEP Lambert arcs;
3. retain candidates with at least one acceptable first-arc departure energy;
4. record the best surviving T0/TOF point as a hypothesis for L0;
5. derives its comparison scale from the best direct transfer instead of a
   sequence-specific manual delta-v threshold.

Lambert solutions are computed once per distinct first encounter body and then
matched against the first-flyby \(v_\infty\) of each Tisserand candidate. The
default acceptance limits are dimensionless ratios to the direct reference and
the Tisserand \(v_\infty\), not hard-coded mission delta-v thresholds.

```python
assessments = plans.filter_scout_sequences_arc1(
    candidates,
    t0_bounds_days=[0, 2000],
    accepted_only=False,
)
```

On the initial Vanilla Kerbin→Jool check, 44 of 58 Tisserand candidates survive
the coarse first-arc filter and KEKKJ survives with a departure-energy ratio of
about 0.46 relative to the best direct Kerbin→Jool grid point.

## Fixed T0-bin scan

`scan_scout_sequence_bins()` implements the first validation harness. It uses a
single direct reference optimized over the complete requested T0 horizon, then
evaluates each Tisserand sequence independently in fixed-width bins.

The comparison now uses the canonical finite-SOI parking-orbit ejection model,
including the cheaper of direct inclined ejection and planar ejection plus an
SOI normal correction. A candidate survives a bin only if:

- a zero-revolution Lambert first arc exists in that bin;
- its arrival \(v_\infty\) matches the retained Tisserand witness;
- its exact first-leg ejection is no more than the configured ratio to the
  globally best direct ejection.

The initial Vanilla Kerbin→Jool scan used J0–2000, 100-day bins, a 10-day T0
grid, 40 TOF samples and a 1.05 ejection ratio. It produced 612 sequence/bin
rows covering 52 unique assisted sequences in about 3.9 seconds, after a
1.2-second Tisserand scout. The global direct reference was 1965.3 m/s from a
100 km Kerbin parking orbit. KEKKJ survived in 14 bins; its best coarse entry
was T0 J550, first-leg TOF 116 days and 1203.2 m/s ejection.

This is still a recall-oriented scan. Multiple Tisserand witnesses and finer
within-bin refinement remain possible future improvements.

## SQLite workflow

Schema v16 persists the complete scout definition, its direct-transfer
reference, every selected sequence/bin candidate, the L0 jobs, and the
L0-to-funnel lineage. Optimizer strategy is a job parameter, so workers can run
the batch without an out-of-band strategy choice.

```python
prepared = plans.prepare_sequence_scout_sqlite(
    "wayfinder.sqlite",
    name="kerbin_jool_0_1000",
    start_body="Kerbin",
    target_body="Jool",
)
plans.optimize_sequence_scout_stage_sqlite(
    "wayfinder.sqlite", prepared["scout_run_id"], stage="l0",
)
promotion = plans.promote_sequence_scout_sqlite(
    "wayfinder.sqlite", prepared["scout_run_id"],
)
plans.optimize_sequence_scout_stage_sqlite(
    "wayfinder.sqlite", prepared["scout_run_id"], stage="funnel",
)
```

The default policy creates all L0 jobs belonging to the 20 best unique scout
sequences, then promotes the two best completed L0 jobs in each 100-day T0
bin. Promotion is refused while L0 jobs remain pending or running. The promoted
run loads the persisted final L0 populations, applies the normal
phase-diversity selection, and starts at `wide`; it does not pay for L0 twice.

```text
sequence_scout_runs
  -> sequence_scout_candidates
  -> sequence_scout_jobs(stage=l0) -> jobs -> runs/results/populations
  -> sequence_scout_jobs(stage=funnel, parent_run_id=...)
  -> promoted jobs -> continuation runs
```

`sequence_scout_status_sqlite()` reports job counts for both stages and an
effective workflow status. `allow_partial=True` exists only for diagnostics and
smoke tests.
