# Optimizer audit — 2026-06-19

## Scope

This audit covers the Pygmo archipelago construction, migration error handling,
funnel population sizing, stage transitions, optimizer persistence and the
controlled stage-3 SADE/Nelder–Mead benchmark.

## Critical finding: topology had twice as many vertices as islands

Wayfinder created `pg.ring(n)` or `pg.fully_connected(n)` and then asked the
archipelago constructor (or `push_back`) to add `n` islands. Pygmo adds a new
topology vertex for every added island, producing `2n` vertices for `n` migrant
databases. Migration then raised `get_migrants` index errors.

The optimizer used `wait()`, which clears worker errors without rethrowing them.
Consequently a run could be stored as `DONE` after migration failed.

Corrections:

- topology factories now return empty topologies (`pg.ring()` and
  `pg.fully_connected()`), which Pygmo grows with the islands;
- optimizer synchronization now uses `wait_check()` so worker failures fail the
  run and are persisted as such;
- regression tests assert both island and topology vertex counts and execute a
  migration-capable evolution.

Validation with eight islands:

- ring: 8 vertices and 16 directed edges;
- fully connected: 8 vertices and 56 directed edges;
- after four ring evolution/migration calls, 37/96 final population IDs came
  from another island; unconnected produced 0/96.

## High finding: invalid SADE population sizes

The debug preset used one individual per island and the funnel allowed focused
populations as small as two. Pygmo SADE requires at least seven individuals.

Corrections:

- the debug preset now requests seven individuals;
- execution clamps old/custom SQLite jobs to the SADE minimum of seven;
- all funnel stages enforce the same invariant;
- a regression test evolves SADE with the resolved population.

## Historical benchmark impact

The workspace contains 1,014 completed runs across seven databases using the
malformed ring/fully-connected construction. Their decoded trajectories may
remain physically valid, but conclusions about topology, migration efficiency,
cooperative funnel behaviour and runtime-to-convergence are invalid.

In particular, the previous conclusion that ring was the best topology must be
rebenchmarked. Unconnected runs were not affected by the vertex mismatch.

## Stage-3 eight-island benchmark

Protocol: eight native Pygmo multiprocessing islands, population 14, five
migration epochs, maximum 7,000 evaluations per island; branches were 100%
`pg.sade`, alternating `pg.sade`/`pg.nlopt("neldermead")`, and 100%
`pg.nlopt("neldermead")`. All branches shared identical stage-2 genes.

The source stage-2 populations predate the topology fix. This remains a valid
controlled comparison of stage-3 portfolios, but it is not an end-to-end
validation of the corrected funnel.

Results across KEKKJ and KEEMo (five seeds each), using optimizer fitness:

- 100% SADE: 0/10 wins;
- alternating: 6/10 wins;
- 100% Nelder–Mead: 4/10 wins.

Using decoded objective DV instead gives a 5/10 split between alternating and
100% Nelder–Mead. Median wall times are nearly identical because all eight
islands now execute in the real Pygmo process pool.

## Stage-2 phase-elite transition and end-to-end benchmark

Stage 1 now exports its complete populations. Stage 2 ranks all candidates,
keeps a pool covering the best 35%, and greedily selects 112 candidates by
max-min distance. The distance features are normalized body position and
velocity at every encounter plus normalized leg times. The selected candidates
fill 8 populations of 14 without random replacement candidates.

Two independent isolated stage-2 repetitions (five seeds per scenario) found
the mixed distribution to be the most robust of the tested random, mixed,
anchored and hybrid transitions. The selected production candidate was then
tested end to end against the corrected legacy funnel with paired Pygmo seeds
and the correct flyby objective (arrival v-infinity reported, not penalized):

- KEKKJ median decoded DV: 2,999.7 -> 2,949.5 m/s; phase funnel won 4/5 seeds;
- KEEMo median decoded DV: 2,156.7 -> 2,258.3 m/s; phase funnel won 2/5 seeds;
- KEKKJ median runtime: 5.14 -> 2.67 s;
- KEEMo median runtime: 4.39 -> 2.38 s.

The phase transition is therefore consistently faster and strong on KEKKJ,
but is not uniformly better on KEEMo at the small benchmark budget. Its best
KEEMo result was nevertheless better (2,055.0 versus 2,108.5 m/s).

A higher-budget qualification of the complete selected chain used 16 -> 8 -> 8
islands, 20/40/5 epochs and about 792,000 nominal evaluations per run:

- KEKKJ: best 2,706.1, median 2,925.1, worst 3,116.2 m/s; 4/5 below 3,100;
- KEEMo: best 1,996.0, median 2,075.5, worst 2,471.8 m/s;
- median runtime was 12.85 s for KEKKJ and 11.26 s for KEEMo.

This qualification also exposed and corrected a benchmark-harness error where
new jobs had temporarily been recreated with `arrival_mode="vinf"` instead of
`arrival_mode="flyby"`. The invalid artifacts were overwritten by corrected
runs.

The exact-ejection audit then found that the stage-3 legacy vector model used a
fixed inertial reference while trajectory decoding used the departure body's
local prograde direction. Exact fitness and decoding now call one canonical
function, and exact mode no longer adds the approximate inclination penalty a
second time. Invalid trial ejection geometries receive a 99,999 m/s rejection
cost instead of raising inside a multiprocessing worker.

After the fix, final exact fitness equals decoded objective DV for every one of
the ten qualification runs (observed difference 0.000000 m/s). The large and
inconsistent KEKKJ transition jumps disappeared: the remaining stage-2 to
stage-3 changes range from -3.3 to +124.6 m/s and represent the real difference
between approximate and local-frame ejection costs. One KEEMo basin still has a
+370.4 m/s correction, which is now a legitimate approximation error rather
than a coordinate-frame mismatch.

The final model explicitly minimizes two physically distinct strategies:

- one direct inclined burn from an equatorial prograde parking orbit;
- one zero-inclination escape using the planar v-infinity projection, followed
  by an absolute vertical-v-infinity correction at the SOI.

Both are closed-form calculations using asymptotic Pykep v-infinity. Fitness
does not call an ephemeris; the departure ephemeris is evaluated only when
decoding the final trajectory angle. A 100,000-call microbenchmark measured
65.6 microseconds per ejection evaluation. Final qualification selected one
SOI-split solution and four direct solutions in each scenario.

Final full-chain qualification:

- KEKKJ: best 2,366.4, median 2,873.5, worst 2,920.8 m/s; 5/5 below 3,100;
- KEEMo: best 1,960.7, median 1,990.2, worst 2,393.3 m/s;
- median runtime is 16.89 s and 14.96 s respectively;
- final fitness minus decoded objective is 0.000000 m/s for all ten runs.

### Stage-budget experiment: 20/40 versus 25/25

The equal-stage variant increases nominal evaluations from about 792k to 836k.
Its stage-1 result improved on every paired seed, confirming that stage 1 was
short. However, reducing stage 2 from 40 to 25 epochs harmed final robustness:

- 25/25 won only 2/5 paired final runs in KEKKJ and 2/5 in KEEMo;
- KEKKJ median was effectively unchanged: 2,873.5 versus 2,873.9 m/s;
- KEEMo median worsened from 1,990.2 to 2,016.8 m/s;
- runtime medians were similar within normal machine variance.

Conclusion: stage 1 benefits from more work, but stage 2 cannot simply be cut
to the same epoch count. The next useful candidate is a longer stage 1 with an
adaptive stage-2 ceiling, or a fixed 25/35 plan.

### Controlled level-1 algorithm benchmark

Level 1 was rerun apples-to-apples with one fixed T0 bin per scenario, the
same initial populations for each paired seed, a 16-island directed ring,
population 32, 25 epochs and equal fitness-evaluation budgets. All algorithms
are native Pygmo: 100% `pg.sade`, alternating `pg.sade` / `pg.simulated_annealing`,
and 100% `pg.simulated_annealing`. Five seeds were used and execution order was
rotated to reduce thermal bias.

- KEKKJ, T0 1100-1200: SADE won 4/5 seeds and had median exact DV 2,603.5 m/s;
  the mixed branch won 1/5 with median 2,679.2 m/s; SA won 0/5 with median
  3,332.6 m/s.
- KEEMo, T0 400-600: the mixed branch won 4/5 seeds and had median exact DV
  2,166.7 m/s; SADE won 1/5 with median 2,299.7 m/s; SA won 0/5 with median
  2,593.4 m/s.
- Runtime medians were effectively equal because evaluation budgets were
  balanced: about 7.8 s for KEKKJ and 7.0 s for KEEMo.
- The mixed branch retained more useful phase diversity than pure SADE. Pure
  SA sometimes retained still more diversity, but at substantially worse DV:
  diversity without basin quality is not sufficient.

There is no universal winner between pure SADE and the 50/50 portfolio at
level 1: KEKKJ favours SADE while KEEMo favours the portfolio. Pure SA is not a
viable level-1 default. A 75/25 SADE/SA portfolio is the next sensible
compromise to test, followed by an end-to-end comparison because level 1 must
be judged partly by the basins it supplies to level 2.

The subsequent 75/25 test used 12 regularly spaced SADE islands and four SA
islands under the same protocol. It did not provide the expected compromise:

- KEKKJ: median exact DV 2,643.6 m/s and 0/5 paired wins;
- KEEMo: median exact DV 2,222.7 m/s and 1/5 paired wins;
- runtime remained identical to the other branches.

The complete rerun also exposed an important reproducibility limitation. The
seeded control branches did not reproduce their first-run values exactly.
Pygmo's multiprocessing archipelago migrates asynchronously, so worker timing
changes which migrants are available when another island evolves. Seeds and
identical initial populations therefore control the stochastic algorithms but
do not make the complete migration history deterministic. Future claims about
systematic convergence must include repeated runs of each seed (or use a
synchronous research harness), especially when tuning migration parameters.
On present evidence 75/25 should not replace either pure SADE or 50/50.

The topology audit found 16 vertices and 32 directed edges in every run. A
subsequent migration-policy audit uncovered that policies passed to an empty
`pg.archipelago` are not inherited by islands later added with `push_back()`.
Those islands silently used Pygmo's default migration rate of one. Policies
are now attached directly to every `pg.island`; telemetry then correctly
reports 0, 16, 32 and 64 published individuals per epoch for rates 0, 1, 2 and
4 on a 16-island ring. Earlier funnel and level-1 results described in this
document therefore used an effective rate of one even when the harness stated
`select_best(2)`.

Two implementation defects were corrected during this audit:

- SADE now uses `memory=True`, preserving its adaptive state across the 25
  separate `evolve()` epochs instead of restarting adaptation every epoch.
- migration-origin telemetry maps an individual ID to every island already
  containing it. The previous one-island map falsely counted a migrant copied
  to multiple islands as a new migration on subsequent epochs.

### Controlled level-1 migration-rate benchmark

After correcting policy ownership, rates 0, 1, 2 and 4 were tested with
matching `select_best(rate)` and `fair_replace(rate)`. The protocol retained
the fixed bins, paired initial populations, 16-island ring, 25 epochs, five
seeds and equal evaluation budgets. Only pure SADE and alternating SADE/SA
were included.

- KEKKJ strongly requires migration but favours a modest rate. Pure SADE at
  rate 1 reached median exact DV 2,557.7 m/s; 50/50 at rate 2 reached 2,559.1
  m/s. Rate 0 degraded to 2,952.1 and 3,081.9 m/s respectively.
- KEEMo benefits from more aggressive sharing. Pure SADE at rate 4 reached
  median 1,919.4 m/s, versus 2,220.2 at rate 2 and 2,400.6 without migration.
  The 50/50 branch also preferred rate 4, but was weaker at 2,157.9 m/s.
- Zero migration retained the largest phase-distance diversity, but this
  diversity remained trapped in weaker basins. Migration trades raw diversity
  for propagation of useful basins.
- There is no universal migration-rate optimum. Rate 1-2 is robust for KEKKJ;
  rate 4 is compelling for KEEMo. A production default of 2 is a reasonable
  neutral setting, but adaptive or sequence-dependent rates deserve testing.

### SADE parameter screening

Eight native Pygmo SADE configurations were screened for 10 epochs with five
paired seeds, both fixed-bin scenarios, 16 islands and migration rate 2. The
mutation variants were rand/1/exp (2), rand/1/bin (7), rand/2/exp (5) and
rand/2/bin (10), each with jDE (adaptation 1) and iDE (adaptation 2).

- Variant 2 / adaptation 1 was best on KEKKJ: median exact DV 3,052.5 m/s.
- Variant 2 / adaptation 2 was best on KEEMo: median 2,205.6 m/s and the best
  mean rank across scenarios.
- Variant 7 / adaptation 1 ranked second on KEEMo and retained the steepest
  late KEKKJ slope (-265 m/s per epoch), so it must be retained as a possible
  late-converging configuration.
- Rand/2 variants 5 and 10 were generally weaker; adaptation 2 was especially
  poor for variant 5 and ranked last in both scenarios.

The 25-epoch finalists are therefore v2/a2, v2/a1 and v7/a1. The screening
rank is deliberately provisional because every leading curve was still
descending at epoch 10.

The 25-epoch qualification produced a scenario-dependent result:

- KEKKJ: v2/a1 was decisively the most robust, with median 2,571.2 m/s and a
  narrow 2,546.4-2,616.3 range. Both v2/a2 and v7/a1 failed into weak basins on
  two seeds, producing medians near 2,728 and worst cases above 3,400 m/s.
- KEEMo: v2/a2 dominated with median 1,849.2 m/s, best 1,840.0 and four seeds
  near or below 2,028. v7/a1 was second at 2,015.6; v2/a1 reached 2,176.0.

Therefore v2/a1 remains the safest general default, while v2/a2 is a strong
KEEMo/JNSQ specialist. The next inexpensive architectural test is an island
portfolio alternating v2/a1 and v2/a2, which may retain KEKKJ robustness while
allowing the a2 basins to propagate on KEEMo.

### Simulated-annealing parameter screening

Before testing that SADE portfolio, the SA half of the existing 50/50
SADE/SA archipelago was screened using only native
`pg.simulated_annealing`. SADE remained fixed at v2/a1 and migration at rate 2.
Three starting temperatures (1,000, 3,000 and 9,000) were crossed with initial
mutation ranges 1.0 and 0.5. The temperature ratio remained constant at
`Ts/Tf = 300`; five paired seeds and 10 epochs were used.

- T1000/R1 was best on KEKKJ at median exact DV 3,055.3 m/s.
- T1000/R0.5 was best on KEEMo at 2,425.6 m/s and had the best combined rank.
- The existing T3000/R1 setting ranked fifth of six in both scenarios, at
  3,362.6 and 2,586.4 m/s respectively.
- T9000/R1 ranked third overall and is retained as a deliberately hotter
  control; reducing its range to 0.5 was not useful.

The SA qualification candidates are therefore T1000/R0.5, T1000/R1 and
T9000/R1. As with the SADE screen, their curves are still descending at epoch
10, so the production SA default must not change before the 25-epoch result.

A direct 25-epoch combination test then compared the former 50/50 defaults
(SADE v2/a1 plus SA T3000/R1) with SADE v2/a2 plus the screening winner SA
T1000/R0.5. Migration remained fixed at rate 2 so only algorithm parameters
changed.

The tuned combination did not generalize:

- KEKKJ median worsened from 2,582.8 to 2,657.9 m/s; the former default won
  3/5 paired seeds.
- KEEMo median worsened from 2,184.1 to 2,269.6 m/s; the former default won
  4/5 seeds.
- KEEMo champion phase distance collapsed from 0.238 to 0.037 despite faster
  approximate-fitness convergence.

Thus independently strong SADE and SA settings are not additive inside a
migrating portfolio. The tuned pair over-concentrates the archipelago into a
basin that looks attractive to the approximate objective but rescored worse
under the exact ejection model. The former defaults remain preferable for the
50/50 portfolio pending individual 25-epoch SA qualification or a less
homogenizing portfolio design.

### Ten-seed full-chain qualification

The retained former level-1 portfolio (SADE v2/a1 plus SA T3000/R1) was then
run through the complete phase-elite funnel without a comparison branch. All
ten seeds used the same controlled T0 bin per scenario: 1100-1200 for KEKKJ
and 400-600 for KEEMo. The funnel used 16/8/8 islands, 20/40/5 epochs,
phase-diverse elites at level 2 and alternating Pygmo SADE/NLopt-NM at level 3.

- KEKKJ: best 2,364.9, median 2,452.5, worst 2,881.8 m/s; 10/10 below 3,100,
  9/10 below 2,700 and 6/10 below 2,500. Median runtime was 12.74 s.
- KEEMo: best 1,925.9, median 2,044.4, worst 2,305.3 m/s; 10/10 below 2,400,
  8/10 below 2,200 and 6/10 below 2,100. Median runtime was 11.32 s.
- Median stage best values were 2,524.3 -> 2,376.7 -> 2,452.5 m/s for KEKKJ
  and 2,082.8 -> 1,895.2 -> 2,044.4 m/s for KEEMo. The level-3 increase is
  the expected transition from approximate to exact ejection fitness, not a
  loss measured under one identical objective.

This larger fixed-bin sample confirms that the complete chain is operational
and consistently reaches useful basins. The remaining spread is concentrated
in one weaker KEKKJ seed (2,881.8 m/s) and a few KEEMo exact-model corrections,
rather than wholesale convergence failures.

### Wide scout, exact archive and adaptive exact refinement

The experimental `funnel_scout_archive` strategy addresses three distinct
failure modes without replacing the useful intermediate level:

- an unconnected 64-island level 0 performs a wide, shallow SADE scout;
- 16 phase-diverse scout champions initialize level 1;
- levels 1 and 2 periodically rescore only island champions and a few diverse
  elites with the exact ejection model;
- a bounded archive preserves candidates that are simultaneously good under
  exact scoring and dispersed in phase space;
- level 3 is initialized from that archive plus the final level-2 population
  and may stop adaptively after the exact champion plateaus.

The scout is budget-neutral with respect to the former level 1: level 0 uses
64 x 8 x 5 x 50 evaluations and level 1 uses 16 x 32 x 15 x 50, whose sum is
the former 16 x 32 x 20 x 50 budget. The extra runtime therefore comes mainly
from periodic exact rescoring and the longer, adaptive exact level 3, not from
silently granting the broad search a larger approximate-fitness budget.

On the same fixed bins and ten paired seeds, the complete strategy produced:

- KEKKJ: best 2,365.2, median 2,403.8 and worst 2,555.3 m/s, versus 2,364.9,
  2,452.5 and 2,881.8 for the reference funnel. It won 6/10 paired seeds and
  removed the reference funnel's weak tail, including a 407 m/s improvement
  on its worst seed.
- KEEMo: best 1,841.8, median 1,958.5 and worst 2,014.8 m/s, versus 1,925.9,
  2,044.4 and 2,305.3. It won all 10 paired seeds.
- Median runtime increased from 12.74 to 19.39 s on KEKKJ and from 11.32 to
  17.47 s on KEEMo, approximately 52-54%.

The exact archive recovered a candidate better than every member of the final
level-2 population in 3/10 KEKKJ runs (about 34, 65 and 35 m/s) and 1/10 KEEMo
runs (about 1 m/s). This validates the archive mechanism independently of the
overall gains: level 2 can converge away from an earlier candidate that becomes
preferable once the ejection is evaluated exactly.

The initial adaptive condition required both champion and population-average
plateaus and therefore consumed all ten exact epochs. The production condition
now follows the exact champion only, after a five-epoch minimum, while retaining
a ten-epoch ceiling. Replaying the stored histories projects earlier stops in
7/20 runs. Keeping the ceiling is justified because epochs 6-10 occasionally
recovered another 52 m/s on KEKKJ or 31 m/s on KEEMo even though their median
gain was only 9 and 6 m/s respectively.

A two-dimensional ejection surrogate indexed by planar and normal v-infinity
is retained as priority 3. It should replace the coarse inclination penalty
only after the wide scout and archive architecture are tuned; it is expected
to approximate the exact parking-orbit ejection much more closely without
putting exact geometry in every level-1/level-2 fitness evaluation.

#### Equal-budget level-0 width screen

The unconnected scout was screened at 32, 64 and 128 islands on ten fixed-bin
seeds per scenario. All variants used populations of eight and five epochs;
SADE generations were respectively 100, 50 and 25, holding level-0 evolution
at 128,000 evaluations. Levels 1-3 were unchanged.

- 32 islands was strongest and tightest on KEKKJ (median 2,384.4, worst
  2,509.0 m/s), but required about 33-36 s/run and produced one 2,307 m/s
  KEEMo outlier.
- 64 islands remained the best compromise across both problems. Its repeated
  qualifications show that asynchronous downstream migration introduces
  meaningful run-to-run variation even with identical algorithm seeds; the
  deterministic level-0 medians were identical while later-stage medians
  diverged.
- 128 islands was rejected: KEKKJ worst case reached 2,621.1 m/s and KEEMo
  2,277.6 m/s. Its 25-generation islands are too shallow, while process
  oversubscription adds runtime variance.

The adaptive exact stop was not the general cause of the weak tail: nearly all
outliers in the 64- and 128-island runs consumed the full ten exact epochs.
One 32-island KEEMo outlier did stop at epoch six, so the stopping rule still
needs a basin-quality guard. Production should retain 64 islands for now;
32 remains a promising quality-oriented KEKKJ variant, not a general default.

## Remaining limitations and recommendations

- The initial SciPy synchronous harness has been retired from the reference
  protocol. NLopt Nelder–Mead is available through Pygmo and runs correctly in
  the native Windows multiprocessing archipelago.
- NLopt may stop on an internal numerical limit before consuming `maxeval`.
  This occurred in two KEEMo seed-3 branches and is recorded via actual fitness
  evaluations rather than hidden or padded artificially.
- Migration counts and topology dimensions are not persisted in optimizer
  telemetry. **Resolved in schema v13:** every snapshot stores island count,
  topology vertices/edges, published migrants, active publishing islands and
  observed accepted migrations. Funnel stages store the corresponding totals.
- A run now fails immediately if its topology vertex count differs from its
  population island count.
- Two islands make ring and fully-connected equivalent. Topology comparisons
  should use at least eight islands.
- Re-run the end-to-end topology and funnel benchmarks after this fix before
  choosing production defaults.

## Verification

- 44 unit/regression tests pass;
- the SQL-only debug funnel smoke test passes;
- an eight-island SQL/Wayfinder flat smoke test completes with two snapshots,
  64 final population points and no worker error;
- all 35 Python files compile.
