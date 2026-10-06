# SCFT Low-Composition, Low-Segregation Resolution Audit

## Objective and claim scope

Verify that the SCFT stress-free periods and density profiles used as the
low-end references in Figures 3 and S4 are accurate enough for the reported
model errors. The claim kind is one stress-free lamellar SCFT observable at
one \((f_A,\chi N)\) state.

The audit covers the first three integer \(\chi N\) values above the reported
RPA spinodal for each composition:

- \(f_A=0.20\): \(\chi N=25,26,27\);
- \(f_A=0.25\): \(\chi N=19,20,21\);
- \(f_A=0.30\): \(\chi N=15,16,17\).

## Numerical comparison

The canonical references use requested spacing \(0.1R_g\), contour step
\(ds=0.01\), Anderson(10) with 100 SD(0.1) warmup iterations, and
VariableCell(BB(1.0)) relaxation. The audit retains the field and cell
updaters while refining both the requested spacing to \(0.05R_g\) and the
contour step to \(ds=0.005\). Each state is solved from cell edges 0.98 and
1.02 times its canonical stress-free period.

## Acceptance gates

Each refined start must satisfy:

- `Polyorder.Successful()` convergence;
- residual norm below \(10^{-6}\);
- absolute cell stress no larger than \(10^{-5}\);
- incompressibility RMS and mean-composition error no larger than \(10^{-6}\);
- one primitive A-rich lamellar domain and dominant Fourier mode 1;
- finite free energy and nonuniform density profile.

The state is resolution-accepted when both starts pass, their optimized
periods agree within 0.05%, and the selected fine result differs from the
canonical reference by at most 0.15% in period and 0.001 in cyclically aligned
profile RMS. A failed observable triggers only the additional spatial-only or
contour-only refinement needed to identify its source.

## Resources and stop condition

At most nine independent jobs may run concurrently, with one Julia thread and
one BLAS thread per job. One numerically changed retry is allowed per state.
The campaign stops when all nine states are accepted or when a named failed
gate identifies the smallest required successor. Process exit status alone
does not promote a result.

## Outcome

The canonical references pass the resolution gates at all six audited states
with (f_A=0.25) and 0.30. The three (f_A=0.20) references do not pass:
refining to (dx\simeq0.05R_g) and (ds=0.005) changes their stress-free
periods by 0.257--0.367%, and the (chi N=25) profile changes by an aligned
RMS of 0.00152.

Separate refinements identify contour discretization as the dominant source.
At fixed (ds=0.01), halving the spatial spacing changes the (f_A=0.20)
periods by at most 0.0130%. At fixed spatial spacing, halving the contour step
to (ds=0.005) changes them by 0.257--0.364%. A further fully refined
(dx\simeq0.05R_g, ds=0.0025) calculation changes the (ds=0.005) results
by at most 0.0917% in period and 0.000517 in profile RMS. Thus
(dx\simeq0.05R_g, ds=0.005) is an adequate replacement protocol for these
three states, whereas the current (ds=0.01) canonical values are not.

The verification campaign does not replace the canonical reference table or
regenerate dependent model errors. Such promotion requires recomputing all
reported quantities that use these three SCFT periods and profiles.
