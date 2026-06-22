# Search architecture foundations

## Work-unit contract

A direct MGA-1DSM job is identified by one concrete body sequence, one
launch-epoch (`T0`) bin, one lower/upper time-of-flight bound for every leg,
one arrival objective, and one optimizer configuration.

`T0` is the primary binning axis. Per-leg ToF bounds form the remaining
dimensions of the same optimization box. Total ToF is the sum of the leg
times and is recorded as a result; it is not independently pixelated. This
avoids overlapping jobs that describe different projections of the same
high-dimensional domain.

The SQLite `tof_min` and `tof_max` fields remain useful denormalized sums for
filtering and backward compatibility. For `tof_encoding = direct`,
`leg_tof_bounds_json` is authoritative.

## ToF estimation

The initial direct estimator supports two profiles:

- `planner`: a narrow KSP-MGA-Planner-like scouting envelope;
- `relaxed`: wider intermediate and terminal legs, plus a 1.9--2.1 period
  envelope for repeated Kerbin encounters.

These are search priors, not feasibility proofs. Later sequence scouting may
produce the same per-leg contract without changing the optimizer or storage
layers.

## Arrival objectives

- `vinf`: include arrival hyperbolic excess speed in the objective.
- `flyby`: minimize departure and DSM costs only. Arrival `v_inf` is still
  decoded and stored so solutions remain physically interpretable.
- `circular` / `elliptical`: include the corresponding capture model.

This separation prevents a displayed flyby cost from being confused with a
rendezvous or capture cost.

## Convergence telemetry

Every evolution step stores both the global best objective and the average
objective across every individual on every island. The convergence plot uses
these two series. A falling best with a flat average indicates isolated lucky
solutions; both curves falling indicates population-wide convergence.

Adaptive stopping is optional. It compares relative improvement over a rolling
window. The generic controller may require both champion and population-average
plateaus. The experimental exact L3 intentionally follows the exact champion
only: population-average noise otherwise prevented every early stop. Minimum
and maximum evolution-step limits remain hard safeguards. Runs store the
controller settings, actual step count, and stop reason in SQLite.

## Optimization funnel

For release 1.6, `optimize_sqlite()` deliberately keeps the corrected
three-stage `funnel` strategy as its public default. The newer scout/archive
chain is promising but remains explicit until its 32-island L0 qualification
is repeated. The flat SADE run remains available with
`optimizer_strategy="flat"` for controlled comparisons.

1. `wide`: alternating SADE and simulated-annealing islands explore the full
   job box with the fast approximate ejection model.
2. `intermediate`: the best island champions seed an archipelago with half as
   many islands, smaller populations and twice the evolution-step budget.
3. `exact_ejection`: the best intermediate champions seed a final archipelago
   with at most half as many islands again (and never more than eight). This
   stage re-ranks candidates with the canonical exact parking-orbit ejection
   model. The default retains the broad 500-generation pass.
   Two explicit experimental policies are also available:
   `funnel_local_exact` exact-scores the intermediate population and initializes
   local SADE clouds; `funnel_hybrid_exact` combines one local and one global
   exact island. Both experimental policies use adaptive 50-generation blocks.

The experimental `funnel_scout_archive` family prepends an unconnected SADE
L0, transfers phase-diverse champions to L1, exact-rescores a bounded diverse
archive during L1/L2, and initializes an alternating Pygmo SADE/NLopt-NM L3
from the archive plus the final L2 population. L0 variants expose 32, 64 and
128 islands at an equal 128,000-evaluation budget. The 128-island variant is
rejected; 32 remains the next qualification candidate.

Island counts decrease monotonically whenever more than one worker is
available. Seeding replaces one random individual rather than increasing the
configured population size. The first two stages alternate SADE and simulated
annealing correctly; the v1.5 prototype accidentally inserted the annealing
island twice. PyKEP's configured arrival/capture objective remains active at
all stages. SQLite schema v14 stores per-stage island/population settings,
initialization, topology name, migration rate, exact archive size, adaptive
controller, algorithms, best/average fitness, runtime and stop reason.

## Job ownership

Optimization workers atomically claim jobs with `BEGIN IMMEDIATE`. A claim
sets `RUNNING`, `claimed_at`, `claim_expires_at` and `worker_id`. Expired leases
are recovered on the next claim, live claims may be renewed, and terminal
statuses clear ownership. The runner claims one job at a time, renews after
each Pygmo epoch, and atomically publishes run/result/job only if the original
`worker_id` and `claimed_at` still own a live claim. This contract is the
concurrency boundary for future CLI, service and GUI runners.
