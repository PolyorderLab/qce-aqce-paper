# BVK2 fixed-stiffness Figure 6 campaign contract

- **Objective:** isolate the effect of the BVK2 adaptive stiffness by comparing
  the production adaptive functional with the matched fixed-stiffness control
  (`adaptive=false`, hence `K_psi=K_psi0`) at `(f_A, chi N)=(0.50,20)` and
  `(0.35,30)`. Both controls use the same filtered UD-theta discretization.
- **Claim kind:** paired lamellar stress-free roots and their locally relaxed
  energy/stress branches.
- **Accepted outputs:** one accepted adaptive root and one accepted
  fixed-stiffness root per state; the corresponding root profiles and branch
  samples; the regenerated
  `results/lamellar_cell_stress_mechanism/{scan.csv,summary.csv,lamellar_cell_stress_mechanism.svg,README.md}`.
- **Functional/discretization provenance:** `c2=0.16`, `nx=256`,
  `oversample=2`, `sensor_filter_ratio=12`, UD-theta coordinate, analytic
  filtered cell stress. The adaptive and fixed controls must differ only in
  the `adaptive` flag.
- **Field gates:** finite energy; certified fixed-cell convergence;
  projected-force RMS `<=1e-5`; projected-force maximum `<=2e-5`;
  `|mean(phi_A)-f_A|<=2e-8`.
- **Morphology gates:** ordered one-period lamella; contrast `>1e-4`; energy
  below the homogeneous state by `>1e-10`.
- **Cell/root gates:** analytic stress magnitude `<=1e-6`; oriented stress
  bracket with logarithmic width `<=1e-5`; root is not factor-boundary limited;
  both `+/-1%` local checks are valid, have the correct stress orientation,
  and do not lower the energy density by more than `1e-8`.
- **Branch gates:** every plotted point is independently field-relaxed with
  the same functional and discretization; every point passes the field,
  composition, and morphology gates; each promoted root is interior to the
  plotted branch and is its local energy minimum within `2e-6` after the
  within-model normalization used by Figure 6.
- **Provenance gates:** record source adaptive root/profile, command, code
  revision/worktree state, parameters, log, terminal artifacts, and scientific
  decision. Diagnostic or failed attempts are not promoted.
- **Resources:** at most two state jobs concurrently, one Julia thread per
  job, isolated state/output directories, and at most 2 h wall time per job.
- **Retry budget:** one numerically changed retry per state after diagnosis;
  no identical rerun and no additional scientific state without a revised
  contract.
- **Stop condition:** both states have accepted paired roots and branches,
  Figure 6 outputs are regenerated and visually/structurally validated, all
  targeted contracts pass, and no campaign job remains active. Otherwise stop
  with the unresolved state marked provisional or rejected and the reason
  preserved.

## Completion record

Status: **accepted; stop condition reached**.

| state | fixed-stiffness `D/R_g` | fixed error | adaptive `D/R_g` | adaptive error |
| --- | ---: | ---: | ---: | ---: |
| `(0.50, 20)` | 3.93828481949 | -2.655% | 4.05353516058 | +0.194% |
| `(0.35, 30)` | 4.25328992441 | -4.266% | 4.41125338548 | -0.711% |

Both fixed roots passed on the first attempt. Their analytic root stresses were
`-5.56e-9` and `-1.73e-8`; the independent `+/-1%` checks had oriented stress
pairs `(-0.05310,+0.05325)` and `(-0.06531,+0.06544)`. The promoted paired
branch table contains 14 independently relaxed points for each BVK2 control at
each state, and every point passes the field-force, composition, and one-period
morphology gates. The final ten-row Figure 6 summary accepts all five models at
both states. No retry or campaign extension was required.

The matched result supports the bounded causal statement that activating the
adaptive stiffness moves the cell-stress zero toward SCFT at both displayed
states. It does not establish that adaptive stiffness uniformly improves every
density profile or alone determines the nonlamellar phase diagram topology.

## Publication hardening rebuild

The final publication rebuild regenerates all five model branches at both
states under `lamellar-cell-stress-mechanism-v2`; mixed-schema partial input is
rejected. The Julia stage validates the complete model/state inventory,
one-domain morphology, field-force gates, source settings, relative provenance,
interior roots, stress orientation, and energy minima before atomically
promoting `scan.csv`, `summary.csv`, and `README.md`. The Python renderer is the
sole writer of the SVG and validates the same inventory before atomic
replacement. This rebuild does not expand the scientific state set or retry
budget.

The first hardening rebuild failed closed because its legacy UD root tolerance
excluded an already accepted asymmetric root. The next diagnostic rebuild
exposed incorrect force-gate labels for the UD and BURP rows during independent
post-validation; those artifacts were superseded and are not the publication
outputs. The corrected rebuild regenerated 140 scan rows (14 for each of 10
model/state branches) and 10 accepted summary rows, all under the v2 schema
with unique inventory and passing force, morphology, and provenance gates. It
was independently validated and accepted; no numerical root retry or
scientific-scope extension was used.

The fixed-root publisher now writes candidate root and profile tables to an
isolated staging directory, verifies the source resolution gate and SHA-256
provenance, and promotes only an accepted pair. Rejection and simulated
interruption after the first file replacement both leave the prior canonical
pair intact. Figure 6 data promotion likewise uses an in-progress marker and a
content-addressed generation manifest for `scan.csv`, `summary.csv`, and
`README.md`; the renderer refuses an incomplete or fingerprint-mismatched
generation.

The SCFT reference loader also fails closed: the unique reference row must be
accepted with passing field and cell gates and a finite positive period, and
its profile must contain exactly 256 finite samples on the complete,
nonduplicated periodic grid with matching state and period metadata. Mutation
tests cover rejected and failed-gate rows, nonfinite periods and profiles, and
duplicate or missing profile indices and cardinality.
