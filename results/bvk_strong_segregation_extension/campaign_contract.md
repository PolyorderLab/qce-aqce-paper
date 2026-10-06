# Strong-Segregation Lamellar Benchmark Extension

## Scientific objective

Extend the independently optimized AB-diblock lamellar comparison to
\(\chi N=40,45,50\) at \(f_A=0.20,0.25,\ldots,0.50\). The extension tests
the approach of the adaptive BVK stiffness from its weak-segregation value
\(1/6\) toward its limiting value \(1/4\). It is reported separately from the
published Liu-2019 domain, \(\chi N\le 35\).

## Claim and accepted outputs

- Claim kind: `stress_free_lamellar_observable` at one composition and one
  segregation strength.
- Canonical states: 21 composition--segregation pairs.
- Required models: SCFT, statewise-mapped OPF, UD, BURP, fixed-stiffness BVK,
  and BVK2. The Liu-2019 OK regression is excluded because its published and
  implemented coefficient domain ends at \(\chi N=35\); extrapolating that fit
  would not be a controlled comparator.
- Accepted outputs: per-state summary and profile tables, a 21-state extension
  ledger, and separate strong-segregation aggregate statistics.
- The existing 133-state common cohort and 146-state Liu-domain ledger remain
  unchanged.

## Acceptance criteria

The SCFT, OPF, UD, BURP, and BVK2 criteria are identical to those recorded
for the accepted Liu-domain campaign, except that OK is not evaluated outside
its fitted domain. BVK2 and its fixed-stiffness control use
\(n_x=256\), twofold Fourier oversampling, \(k_c/k_*=12\), and the same field,
cell-stress, one-period morphology, composition, and local-period-minimum
gates. A zero process exit is not scientific acceptance.

## Resources and retries

- Pilot: \(f_A=0.50\) at \(\chi N=40,45,50\), 12 Julia threads.
- Full extension: independent compositions may run concurrently, with no more
  than seven composition workers and 96 total Julia threads.
- Fixed-stiffness roots may run concurrently after their accepted BVK2 and
  SCFT sources exist.
- Retry budget: two numerically distinct retries per rejected or provisional
  state. Every retry must reuse a compatible accepted profile or change a
  documented numerical hypothesis.
- After the initial retry budget was exhausted at
  \((f_A,\chi N)=(0.50,45)\), the user explicitly authorized reopening this
  point and resuming the full campaign. The renewed retry uses repeated
  VariableCell--BB polishing passes from the best finite in-memory state; it
  does not relax any scientific acceptance threshold.

## Stop condition

Stop when all 21 states contain accepted rows for every required model, the
fixed-stiffness controls pass the matched production gates, Figures 3 and 4
show the extension distinctly from the Liu-domain data, and Figure 5 reports
the original and strong-segregation cohorts separately.

## Decision history

- `pilot-f050-40-45-50`: infrastructure/contract failure after the accepted
  SCFT and OPF calculations at \((f_A,\chi N)=(0.50,40)\). The workflow then
  invoked the Liu-2019 OK regression, whose explicit domain gate rejects
  \(\chi N>35\), before the case checkpoint was written. No scientific row was
  promoted. The successor removes only this out-of-domain comparator.
- `pilot-f050-40-45-50-r1`: numerical retry. The process completed and all
  reduced-model rows passed, as did the SCFT references at \(\chi N=40\) and
  50. At \(\chi N=45\), however, `cell_solve!` returned an acceptable rather
  than successful termination with normalized cell stress
  \(1.80328\times10^{-5}\), above the \(10^{-5}\) acceptance threshold. The
  successor recomputes only this SCFT state in a staging directory, using
  16000 cell iterations, unit cell-update blocks, and a stricter
  \(5\times10^{-6}\) stress target. Accepted pilot checkpoints are retained.
- `pilot-f050-chi45-scft-r2`: rejected numerical retry. Updating the cell
  after every field iteration destabilized the BB relaxation, produced a
  nonfinite cell edge for all viable starts, and emitted no scientific
  artifact despite the wrapper's zero exit code. The final permitted retry
  restores the validated five-iteration cell block and warm-continues from a
  freshly converged \(\chi N=50\) state to \(\chi N=45\), retaining the
  enlarged 16000-iteration cell budget.
- `pilot-f050-chi45-scft-r3`: blocked. The intended \(\chi N=50\) warm-start
  calculation became cell-unstable and produced no checkpoint. The target
  \(\chi N=45\) calculation converged its density field to a residual of
  \(4.654\times10^{-15}\), with acceptable composition, resolution, and
  morphology, but terminated `Polyorder.Acceptable()` with cell stress
  \(2.8903\times10^{-5}\), above the \(10^{-5}\) gate. Both permitted
  numerical retries for this provisional state are exhausted. No staged row
  is promoted and the full extension remains unlaunched.
- `pilot-f050-chi45-scft-polish-r4`: user-authorized renewed retry. Each pass
  begins from the best finite SCFT field and cell produced by the preceding
  pass, while retaining the production Anderson(10)+SD(0.1) field update,
  VariableCell--BB(1.0) cell relaxation, block size 5, spatial and contour
  resolution, and the original \(10^{-5}\) cell-stress gate.
- `pilot-f050-chi45-scft-polish-r5`: r4 was rejected after preserving a provisional
  state with residual \(1.211\times10^{-14}\) and stress
  \(3.9164\times10^{-5}\). Subsequent default-BB restarts became nonfinite.
  The successor retains VariableCell--BB but bounds the polish cell step by
  \(0.05R_g\) and enlarges only the first polish displacement so that the BB
  secant curvature is resolved on the nearly flat stress branch. No field,
  resolution, or acceptance setting is changed.
- `pilot-f050-chi45-scft-polish-r6`: r5 preserved a finite provisional state
  but its second pass raised `SingularException(10)`. The successor follows
  the user's explicit constraint to use only VariableCell--BB with
  Anderson(10)+SD(0.1). Each polish pass now constructs a fresh SCFT and fresh
  Anderson updater, then copies the converged auxiliary fields and cell from
  the preceding pass; solver-history matrices are not reused.
- `pilot-f050-chi45-scft-polish-r7`: r6 showed that the singularity recurs
  with a fresh updater and therefore arises from rank loss in the current
  Anderson(10) least-squares history. The successor keeps
  VariableCell--BB(0.05) with Anderson(10)+SD(0.1) and enables Polyorder's
  built-in Anderson history-column dropping at condition number \(10^{12}\)
  during polish passes. No alternative cell optimizer is permitted or used.
- `pilot-f050-chi45-scft-continuation-r8`: r7 was interrupted before output
  after the user clarified the intended continuation. The successor first
  recomputes the accepted \(\chi N=40\) stress-free period and auxiliary
  fields, transfers both to a fixed-cell \(\chi N=45\) SCFT solve, and then
  invokes VariableCell--BB with Anderson(10)+SD(0.1). No direct-45 polishing
  and no stress-guided optimization are used.
- `pilot-f050-chi45-scft-continuation-r9`: r8 validates the continuation
  branch and its \(\chi N=40\) seed, but the \(\chi N=45\) VariableCell solve
  reaches 8000 iterations with stress \(1.8033\times10^{-5}\). The successor
  changes only the cell-iteration budget to 16000 and retains the same
  fixed-cell continuation, VariableCell--BB, Anderson(10)+SD(0.1), and
  publication gates.
- `pilot-f050-chi45-scft-continuation-r10-bb1`: r9 reproduces the r8 period
  and stress to numerical precision after 16000 iterations, establishing BB2
  step-size stagnation rather than an insufficient iteration budget. The
  successor restores a 1000-iteration cap and changes only the BB secant
  variant from BB2 to BB1. BB1 is less susceptible to an artificially small
  step when the cell-stress secant difference is dominated by residual field
  relaxation. The accepted \(\chi N=40\) fixed-cell continuation,
  Anderson(10)+SD(0.1), cell-update block 5, resolution, and all publication
  gates are unchanged.
- `pilot-f050-chi45-scft-continuation-r11-block10`: BB1 is rejected because
  it reaches the 1000-iteration cap with stress
  \(3.4796\times10^{-4}\), substantially worse than BB2. The field residual,
  incompressibility, composition, resolution, and one-period morphology still
  pass, isolating the defect to the cell update. The successor restores BB2
  and changes only the field-relaxation block from 5 to 10 so that each BB
  secant pair is evaluated after a better-relaxed field state. The
  1000-iteration cap, accepted \(\chi N=40\) continuation,
  Anderson(10)+SD(0.1), and publication gates remain fixed.
- `pilot-f050-chi45-models-r12`: the block-10 BB2 continuation is accepted.
  It converges in 156 VariableCell iterations to
  \(L_0/R_g=4.91413864\), field residual \(9.991\times10^{-7}\), and cell
  stress \(6.774\times10^{-6}\); all field, cell, resolution, composition,
  one-period morphology, and local-minimum gates pass. In accordance with the
  restart policy, this accepted \(\chi N=45\) SCFT state is reused directly.
  The successor performs no SCFT optimization and computes only the reduced
  model rows needed to complete the pilot reference at this state.
- `pilot-f050-fixed-{40,45,50}`: the r12 model-only continuation is accepted;
  SCFT, OPF, UD, BURP, and BVK2 all pass their model-specific gates at
  \((f_A,\chi N)=(0.50,45)\). The accepted case replaces the earlier
  provisional pilot checkpoint. The only remaining three-state pilot
  observable is the fixed-stiffness BVK control. Its three independent roots
  are launched concurrently using the accepted adaptive BVK2 profiles as
  production-resolution seeds; no SCFT or other reduced model is recomputed.
- `full-f{020,025,030,035,040,045}`: all three fixed-stiffness pilot roots
  pass source-fingerprint, field-stationarity, primitive-lamella, cell-stress,
  stress-orientation, and local-minimum gates. The complete \(f_A=0.50\)
  pilot is therefore accepted. The six remaining compositions are launched as
  independent workers, each continuing internally through
  \(\chi N=40\rightarrow45\rightarrow50\) with the validated
  VariableCell--BB2 block-10 SCFT protocol and a 1000-iteration cap. Accepted
  \(\chi N=35\) checkpoints seed each reduced-model history. Total allocation
  is 72 Julia threads, below the 96-thread campaign limit.
- `full-f025`, `full-f035`, and `full-f040`: accepted.  Every requested
  \(\chi N=40,45,50\) state contains accepted SCFT, OPF, UD, BURP, and BVK2
  rows, and all field, cell, resolution, composition, morphology, and
  local-period-minimum gates pass.  No successor is justified for these
  composition branches.
- `full-f020`: numerical retry.  The \(\chi N=40\) and 45 states pass, and
  the \(\chi N=50\) SCFT solve is otherwise converged
  (residual \(8.00\times10^{-7}\), stress \(1.63\times10^{-7}\)), but its
  incompressibility RMS is \(2.38\times10^{-6}\), above the
  \(10^{-6}\) composition gate.  The successor recomputes only this state
  with the same BB2 VariableCell solver, a 20-step field-relaxation block,
  a \(2\times10^{-7}\) internal residual target, and the unchanged
  1000-iteration cap and publication gates.
- `full-f030`: numerical retry.  The \(\chi N=45\) and 50 states pass.  At
  \(\chi N=40\), the one-period branch remains finite but reaches the
  1000-iteration cap with residual \(3.14\times10^{-6}\) and stress
  \(3.06\times10^{-5}\).  The successor recomputes only this state with
  BB2, a 20-step field-relaxation block, and a \(2\times10^{-7}\) internal
  residual target; neither the iteration cap nor any acceptance threshold is
  relaxed.
- `full-f045`: numerical retry.  The \(\chi N=40\) state passes.  The
  \(\chi N=45\) branch has a fully relaxed field but stalls at stress
  \(4.30\times10^{-5}\), and its continuation to \(\chi N=50\) becomes
  nonfinite before producing an artifact.  The successor therefore repairs
  only \(\chi N=45\), retaining BB2 and the 1000-iteration cap while
  increasing the field-relaxation block to 20.  The \(\chi N=50\) state will
  be launched only after the repaired 45 state passes.
- `full-f020-chi50-block20-r1`, `full-f030-chi40-block20-r1`, and
  `full-f045-chi45-block20-r1`: accepted.  The repaired SCFT references
  converge in 49, 10, and 14 VariableCell iterations, respectively.  Their
  residuals are \(7.00\times10^{-8}\), \(1.92\times10^{-7}\), and
  \(1.67\times10^{-7}\); their cell stresses are all below
  \(3.5\times10^{-8}\); and their incompressibility RMS values are below
  \(1.3\times10^{-7}\).  Every SCFT, OPF, UD, BURP, and BVK2 publication
  gate passes.  These accepted staged cases replace only their provisional
  counterparts.  The remaining independent fixed-stiffness roots are
  launched in parallel, while \((f_A,\chi N)=(0.45,50)\) is continued from a
  freshly reconstructed accepted \(\chi N=45\) state in the same process.
- `fixed-f{020,025,030,035,040}-chi{40,45,50}` and
  `fixed-f045-chi{40,45}`: accepted.  All 17 fixed-stiffness roots match the
  exact source-summary and source-profile hashes, contain 256 finite
  primitive-lamella samples, and pass the field, composition, morphology,
  cell-stress, stress-orientation, and two-sided local-minimum gates.  Across
  this slice, the largest projected-force RMS is
  \(3.23\times10^{-9}\), the largest absolute root stress is
  \(5.45\times10^{-7}\), and every \(-1\%\)/\(+1\%\) period perturbation
  brackets the stress zero with the required orientation.  No fixed-control
  successor is justified for these states.
- `full-f045-chi50-block20-r1`: accepted SCFT/OPF continuation.  Starting
  from a freshly reconstructed \(\chi N=45\) state, the \(\chi N=50\) SCFT
  reference converges in 275 VariableCell iterations to
  \(L_0/R_g=5.02630877\), residual \(1.85\times10^{-7}\), cell stress
  \(3.35\times10^{-7}\), and incompressibility RMS
  \(3.40\times10^{-7}\).  SCFT and OPF pass all declared gates.  The
  smallest successor reuses this accepted SCFT checkpoint and computes only
  the remaining UD, BURP, and BVK2 rows (with OPF regenerated for a complete
  atomic case artifact).
- `full-f045-chi50-models-r2`: accepted.  The reused SCFT reference and the
  OPF, UD, BURP, and BVK2 rows all pass their field, cell, resolution,
  composition, morphology, and local-period-minimum gates.  BVK2 gives a
  signed period error of \(-0.812\%\) and profile RMS 0.00909.  The complete
  five-model case is promoted atomically; the only remaining state-level
  calculation is its matched fixed-stiffness BVK root.
- `fixed-f045-chi50`: accepted.  Exact source-summary and source-profile
  hashes match the promoted adaptive BVK2 case; the 256-point profile is
  finite and primitive-periodic; projected-force RMS is
  \(5.06\times10^{-10}\); root stress is \(-1.28\times10^{-8}\); and the
  \(-1\%\)/\(+1\%\) period perturbations give stresses \(-0.08409\) and
  \(+0.08409\).  All 21 requested states now have accepted main-model and
  fixed-stiffness rows.  No further simulation is justified; the successor
  is the deterministic ledger assembler and publication-artifact rebuild.
- `strong-extension-assemble`: accepted.  The deterministic assembler records
  all 21 composition--segregation states, 105 accepted SCFT/OPF/UD/BURP/BVK2
  rows, and 21 matched fixed-stiffness roots.  Independent inventory, gate,
  finiteness, and uniqueness checks found no failures.  The strong-segregation
  aggregate gives mean absolute period errors of 6.35%, 15.09%, 18.88%,
  5.48%, and 1.08% for OPF, UD, BURP, fixed-stiffness BVK, and BVK2,
  respectively; the corresponding profile RMS values are 0.10168, 0.01738,
  0.01763, 0.01379, and 0.01569.  No numerical successor is justified.  The
  smallest remaining successor is deterministic integration of this accepted
  extension into Figures 3--5.
- `strong-extension-render`: stop.  Figures 3 and 4 now show the accepted
  \(\chi N=40,45,50\) extension in a visually distinct shaded region,
  and Figure 5 reports the 133-state Liu-domain cohort and the separate
  21-state strong-segregation cohort in adjacent panels.  Visual inspection,
  byte-identical repeat rendering, and all 12 targeted renderer tests pass.
  The campaign stop condition is satisfied, so no successor is launched.
