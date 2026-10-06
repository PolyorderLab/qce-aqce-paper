# Complete Integer Strong-Segregation Lamellar Campaign

## Objective and scope

Complete the stress-free AB-diblock lamellar benchmark at every integer
\(\chi N=36,37,\ldots,50\) for
\(f_A=0.20,0.25,\ldots,0.50\).  The claim kind is a stress-free lamellar
observable at one composition and segregation strength.  The accepted model
set is SCFT, statewise-mapped OPF, literal UD, BURP, BVK2, and the matched
fixed-stiffness BVK control.  The Liu-2019 OK regression is not extrapolated
beyond its published coefficient domain, \(\chi N\leq35\).

The already accepted \(\chi N=40,45,50\) states are retained.  Only the 84
missing composition--segregation states are recomputed for the main models,
followed by 84 independently optimized fixed-stiffness controls.

## Numerical protocol and acceptance

- Polyorder SCFT uses \(\Delta x\simeq0.1R_g\), \(\Delta s=0.01\), a
  fixed-cell presolve, and `cell_solve!` with
  `VariableCell(BB2, Anderson(10)+SD(0.1))`.  The Anderson warmup is 100
  iterations.  The cell loop is capped at 1000 iterations; increasing that cap
  is not an authorized retry strategy.
- Each composition is continued in ascending \(\chi N\).  Accepted neighboring
  periods and model profiles at \(\chi N=35,40,45,50\) are used as history;
  newly accepted states become the seed for the next integer.
- The internal SCFT residual target is \(2\times10^{-7}\), with 20 field
  iterations per BB cell update.  Publication acceptance remains residual
  below \(10^{-6}\), cell-stress norm at most \(10^{-5}\), incompressibility
  and mean-composition errors at most \(10^{-6}\), and one primitive,
  nonuniform lamellar period.
- OPF, UD, BURP, BVK2, and fixed-stiffness BVK retain the model-specific field,
  cell-stress, resolution, composition, morphology, and two-sided local-period
  gates of the accepted \(\chi N=40,45,50\) campaign.  BVK2 and its fixed
  control use \(n_x=256\), twofold oversampling, and \(k_c/k_*=12\).
- A zero process exit is not scientific acceptance.  Status vocabulary is
  `accepted`, `provisional`, or `rejected`.

## Resources, retries, and stop condition

- Seven independent composition workers may run concurrently, each with 12
  Julia threads; total allocation is 84 threads, below the 96-thread limit.
- Fixed-stiffness roots may be batched by composition and run concurrently
  only after their accepted adaptive BVK2 sources exist.
- Retry budget: two numerically distinct retries per rejected or provisional
  state.  A retry must reuse the latest compatible state and change a
  documented numerical hypothesis; identical repetition and cell-iteration
  caps above 1000 are disallowed.
- Canonical outputs are the per-state case directories, a 105-state integer
  ledger, aligned profiles, a matched fixed-stiffness ledger, and updated
  Figures 3--5.
- Stop when all 105 states have accepted SCFT/OPF/UD/BURP/BVK2 rows, all 105
  matched fixed-stiffness roots pass, the deterministic assembler validates
  the complete inventory, and the figures render and pass visual and
  structural checks.

## Decision history

- `integer-f{020,025,030,035,040,045,050}`: launched as seven independent
  composition continuations.  Each requests only the missing integers
  36--39, 41--44, and 46--49 while reading the accepted neighboring history
  from the canonical extension directory.  No accepted 40, 45, or 50 state is
  recomputed.
- `integer-f{020,025,030,035,040,045,050}`: accepted after completion with
  exit code 0.  Each branch contains 12 complete five-model state records and
  their aligned and native profiles.  All field, cell, resolution,
  composition, morphology, and local-period gates pass.  Across the 84 new
  SCFT states, the largest residual, stress norm, incompressibility error, and
  mean-composition error are respectively (1.999\times10^{-7}),
  (7.096\times10^{-6}), (9.349\times10^{-7}), and
  (3.442\times10^{-15}), all within the declared limits.
- `fixed-integer-f{020,025,030,035,040,045,050}`: selected as the smallest
  successors.  Each batch computes only the 12 missing matched
  fixed-stiffness BVK controls from the newly accepted adaptive-BVK2 sources;
  assembly and figure rendering remain deferred until these controls pass.
- `fixed-integer-f{020,025,030,035,040,045,050}`: accepted after completion
  with exit code 0.  All 84 roots retain exact source hashes and accepted
  adaptive-BVK2 provenance.  The maximum projected-force RMS and component
  are (5.668\times10^{-6}) and (1.470\times10^{-5}), respectively, and
  the maximum cell-stress norm is (3.776\times10^{-7}).  Every profile has
  256 finite samples, one periodic A-rich domain, and dominant Fourier mode
  one.
- `integer-extension-assemble`: selected as the single smallest successor
  after all seven fixed-stiffness branches passed.  The assembler reads the
  105 canonical case directories directly, verifies the main and fixed
  inventories and fixed-root source fingerprints, and atomically publishes
  the complete integer-
  \(\chi N=36,\ldots,50\) ledgers and validation report.  Figure rendering
  remains deferred until this assembled artifact passes.
- `integer-extension-assemble`: accepted after completion with exit code 0.
  Independent validation reproduces the report from 525 main-model rows, 105
  fixed-stiffness rows, 134400 aligned main-profile samples, and 107520 aligned
  fixed-profile samples.  The 105-state inventories and aggregate statistics
  are complete and finite.
- `integer-extension-render`: selected as the smallest successor.  It
  regenerates only the stress-free-period comparison, the profile comparison,
  and the aggregate error bars (main Figures 3--5) from the accepted integer
  ledger.  No additional numerical calculation is performed.
- `integer-extension-render`: numerical, structural, and regression gates
  passed, as did the visual checks for Figures 3 and 5.  The first rendering
  was assigned `retry` because repeated endpoint tick labels overlapped between
  adjacent top-row profile panels in Figure 4.
- `integer-extension-render-r1`: selected as a layout-only retry.  It suppresses
  the repeated left endpoint labels in the second through fourth profile
  panels and rerenders the unchanged accepted ledger.
- `integer-extension-render-r1`: the endpoint-label repair passed, but the
  render was assigned `retry` because the shaded extension region began at
  \(\chi N=37.5\), inconsistent with the first accepted extension state at
  integer \(\chi N=36\).
- `integer-extension-render-r2`: selected as the final semantic-layout retry.
  It moves the shaded boundary to \(\chi N=35.5\), leaving all numerical data
  and axes unchanged.
- `integer-extension-render-r2`: assigned stop after completion with exit
  code 0.  Figures 3--5 pass structural and visual inspection and reproduce
  byte-for-byte on a second render.  The manuscript and SI report the
  105-state integer-\(\chi N=36,\ldots,50\) cohort separately from the
  133-state published-domain aggregate.  Targeted figure and manuscript
  regressions pass, and no campaign job remains active.
