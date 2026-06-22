
### v1.0 :
- working prototype
- data saved as pickle
    
        
### v1.1 :
    
- transition to xlsx for storage to 1) avoid safety issues 2) make the saved data "legible"
- Combinatorial method for sequence generation 
- reset methods for optimisation
- possibility to add new jobs to and existing file
- added edit method to allow re-reunning a job with a different opt level
- various code cleanup

### v1.2 :

- Added plot function to compare sequences
- default set of injection parameters (circular, highly elliptical, Vinfinite or none)
- injection and ejection altitudes are now fully fledged job parameters, and are accepted in the edit method too.
- auto Tof mode added
- added * in sequence for optional swing by
- added proper Ejection DV calculation to display on the flight plan.
	   
        
### v1.3 :

- optimization now works with a decorated MGA problem to help tweak the optimisation, for example : added an ejection burn estimation in the fitness function,
    will force a 2Kyr period for KK sequences except if the sequence is starting with KK. Experimental MoMo braking sequence helper. Curation of high inclination ejections.
- replaced the two kraken patches by a single one with planet pack toggle.
- added a T0 window selector in the sequence filter method
- added on \epsilon0 to ward off some rare errors when computing ejection angles and inclination (minimal precision impact)
- added test with a full ejection burn calculation in the optimisation function and proof it with direct transfers (success)
- improved transX style flight plan display, removal of useless parameters.
            
### v1.4 :

- when a job completes, t0, tof and DV values are stored in the file. This will make plotting much easier
- added DV vs T0 plot for aritrary sequences in all possible forms, including combinatorial. Very nice.
- added a few tests

### v1.6.0

- Fixed Pygmo topology construction: ring and fully-connected archipelagos no
  longer receive twice as many graph vertices as optimizer islands.
- Added the experimental `funnel_scout_archive` strategy: a wide unconnected
  level-0 SADE scout feeds phase-diverse candidates into the retained
  three-level optimization funnel.
- Added a bounded multi-fidelity archive. Levels 1 and 2 periodically rescore
  champions and diverse elites with the exact ejection model, and level 3 is
  seeded from both this archive and the final level-2 population.
- Added champion-based adaptive stopping to the exact level 3, with a minimum
  five epochs and a ten-epoch safety ceiling.
- Qualified the combined scout/archive strategy on ten fixed-bin seeds. Median
  decoded DV improved from 2,452.5 to 2,403.8 m/s for KEKKJ and from 2,044.4
  to 1,958.5 m/s for KEEMo; the latter improved on all ten paired seeds.
- Recorded a planar/normal v-infinity lookup surrogate for parking-orbit
  ejection as the next lower-priority model improvement.
- Replaced optimizer `wait()` calls with `wait_check()` so migration/worker
  failures cannot be silently persisted as successful runs.
- Enforced SADE's minimum population of seven for debug, legacy and custom
  SQLite jobs.
- Added topology-dimension, migration and minimum-population regression checks.
- SQLite schema v13 persists topology dimensions, published migrant counts and
  observed accepted migrations for every optimizer snapshot and funnel stage.
- SQLite schema v14 adds atomic job leases (`claimed_at`, `claim_expires_at`,
  `worker_id`), expired-claim recovery and persisted per-stage topology name,
  migration rate and exact archive size.
- `optimize_sqlite()` now claims one job immediately before execution, renews
  its lease between Pygmo epochs, and atomically publishes the run, result,
  gene and terminal job status only for the current claim. Stale workers are
  fenced out after an expired job is recovered.
- Kept the corrected three-stage `funnel` as the explicit 1.6 public default;
  scout/archive variants remain opt-in until the 32-island result is repeated.
- Added a pinned conda-forge `environment.yml` and documented the Windows
  Conda `Library\bin` DLL requirement for direct Python invocation.
- Optimizer runs now fail immediately when topology vertex and island counts
  differ, instead of relying on a later migration error.
- Stage-3 portfolio benchmarks now use only native Pygmo algorithms:
  `pg.sade` and `pg.nlopt("neldermead")` in the real multiprocessing ring.
- Added the `funnel_phase_elites_nm` strategy: stage 2 is initialized from
  fitness elites selected by max-min diversity in physical encounter-phase
  space; stage 3 alternates `pg.sade` and `pg.nlopt("neldermead")` on 8 islands.
- Funnel populations and stochastic Pygmo algorithms now receive reproducible,
  paired seeds from `optimizer_seed`.
- SADE retains its adaptive state across successive archipelago epochs via
  `memory=True`; previously each monitoring epoch restarted its adaptation.
- Migration-origin telemetry now retains every island containing an ID, so a
  copied migrant is not falsely recounted when it later migrates again.
- Fixed migration-policy ownership for incrementally constructed
  archipelagos. `select_best` and `fair_replace` are now attached to every
  island; policies supplied only to an empty archipelago were silently lost
  when islands were added later, leaving the Pygmo default rate of one.
- Added a fixed-T0-bin level-1 benchmark comparing native Pygmo SADE, a 50/50
  SADE/simulated-annealing ring and simulated annealing alone over five paired
  seeds with equal evaluation budgets.
- Extended the benchmark with a regularly interleaved 75/25 SADE/SA branch.
  It did not improve either scenario, and the rerun documented residual
  timing variance from asynchronous Pygmo migration despite seeded algorithms.
- Added a checkpointed fixed-bin migration-rate benchmark for rates 0, 1, 2
  and 4. KEKKJ preferred rate 1-2 while KEEMo strongly favoured pure SADE at
  rate 4; zero migration was substantially worse in both scenarios.
- Added a checkpointed SADE successive-halving screen across four mutation
  variants and both Pygmo adaptation schemes. Variants 2/adaptation 1,
  2/adaptation 2 and 7/adaptation 1 advance to the 25-epoch qualification.
- The 25-epoch SADE qualification retained variant 2/adaptation 1 as the most
  robust general default and identified variant 2/adaptation 2 as a markedly
  stronger KEEMo specialist but an unreliable KEKKJ default.
- Added a six-configuration Pygmo simulated-annealing screen in the 50/50
  SADE/SA archipelago. Starting temperature 1,000 outperformed the current
  3,000 setting; ranges 0.5 and 1.0 advance with a hotter 9,000/1.0 control.
- A 25-epoch combined test of SADE v2/a2 with SA T1000/R0.5 underperformed the
  former v2/a1 plus T3000/R1 portfolio on both scenarios after exact rescoring,
  despite faster approximate convergence, and sharply reduced phase diversity.
- Qualified the retained portfolio through the complete three-level funnel on
  ten seeds using one fixed T0 bin per scenario. All 20 runs succeeded; median
  decoded DV was 2,452.5 m/s for KEKKJ and 2,044.4 m/s for KEEMo.
- Added the 16 -> 8 -> 8 `benchmark_funnel` preset and an end-to-end benchmark
  harness retaining runtimes, convergence, stages and migrations in SQLite.
- Added a ~792k-evaluation `qualification_funnel` preset and full-chain quality
  harness. Benchmark jobs now explicitly use the planner-compatible `flyby`
  objective (arrival v-infinity reported but not included in objective DV).
- Unified exact ejection fitness and trajectory decoding through one canonical
  local-frame calculation. Exact mode no longer double-counts the approximate
  inclination penalty, and invalid trial geometries receive a finite rejection
  cost instead of crashing a multiprocessing island.
- The canonical ejection model now evaluates two closed-form strategies:
  a direct inclined burn and a zero-inclination escape followed by a normal SOI
  correction. It reports the selected strategy and both component costs, uses
  asymptotic Pygmo/Pykep v-infinity consistently, and performs no ephemeris call
  in the optimizer fitness loop.
- Added an explicit equal-stage phase-elite strategy and a 25/25 qualification
  preset for controlled stage-budget experiments. It remains experimental; the
  20/40 reference was more robust in the paired ten-run comparison.
- SQLite is now the only Wayfinder datastore for jobs, runs, results and genes.
- `add_batch` and `optimize` now require an explicit SQLite path, with `add_batch_sqlite` and `optimize_sqlite` as the preferred workflow.
- SQL result display/query helpers now support `batch_name` filtering to compare batches without mixing binning strategies.
- Regression references moved to a SQLite fixture.
- The old dataframe file load/save path has been removed from the core API.
- The core is now pykep3-only; pykep2 compatibility shims were removed and a regression guard prevents reintroducing them.
- Old Excel datastores were moved under `legacy/`.
- Planet definitions moved under `WayfinderCore/planet_packs/`.
- `_Kraken_Patch.py` was renamed to `_Trajectory.py` now that trajectory decoding/display is standalone and no longer monkey patches pykep.
- Non-flight-plan diagnostics now use the standard `logging` module; `transx` keeps direct console output for readable flight plans.
- Wayfinder optimization fitness decoration was moved out of `optimize` into a dedicated optimization module.
- SQLite now stores optimizer run snapshots with per-step champion fitness and genes.
- SQLite now stores final optimizer populations so non-winning genes can be reused for local porkchop analysis.
- Added a local porkchop plot helper for stored optimizer genes.
- Added explicit optimizer porkchop and sampled local porkchop helpers, with sampled grids stored in SQLite.
- Restored the v1.5 three-stage optimization funnel as the default SQLite strategy, with flat SADE retained as an explicit benchmark option.
- Corrected the staged hybrid island mix (the v1.5 prototype inserted simulated annealing twice), made island counts monotonically narrow, and kept seeded population sizes exact.
- Added a final full-vector 3D ejection-cost stage while preserving the configured PyKEP arrival/capture objective.
- SQLite schema v11 records configuration, convergence, runtime and stop reason for every funnel stage.
- SQLite schema v12 replaces the 500-generation global exact-ejection pass with locally seeded 50-generation SADE blocks, exact champion re-ranking and adaptive stopping.
- The exact stage now scores the whole intermediate population and balances one tight local island with one diverse global island, avoiding the quality loss of a purely local refinement.

### v1.5.0

- added a SQLite datastore layer for jobs, batches, sequences, runs, results and genes.
- wildcard batch templates are stored separately from concrete generated sequences.
- concrete sequences are normalized with start/target/flyby metadata, allowing route-wide queries such as best known Kerbin to Moho result across batches.
- jobs are deduplicated with a canonical parameter hash so different binning strategies can coexist cleanly.
- pykep3 compatibility added while keeping pykep2 regression coverage.

### v1.414

- added examples of use for different cases
- re-structures folders to make things neatier
- smarter auto_tof. Will use 0.5-1 orbital per body in the sequence, except for :
	1) the starting one (not counted)
	2) if "planet-planet" is found in the sequence (but not the first two), 1.0-2.0 "planet" periods are added
	3) for the last body, only 0.2-0.4 period is counted.


