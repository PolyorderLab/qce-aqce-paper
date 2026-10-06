# QCE profile-error resolution check at low A fraction

This auxiliary audit tests the nonmonotonic QCE profile errors at low \(f_A\).
The profile-error comparison is shown in Figure 4. This spatial-refinement
check uses nine calculations covering the three-point
neighborhoods around the sharpest features at \(f_A=0.20\), 0.25, and 0.30.
The campaign contract is recorded in `campaign_contract.md`, and the initial
\(n_x=64\to128\) comparison is stored in `resolution_summary.csv`.

## Result

All nine calculations are accepted. Every state passes field stationarity,
phase identity, the cell-stress root, and the interior local-period-minimum
checks. The largest \(n_x=64\to128\) change is 0.156% in the stress-free period
and 0.00125 in the SCFT-aligned profile RMS, below the predeclared limits of
0.25% and 0.0025.

At \(n_x=128\), the local shape of each feature was:

- \(f_A=0.20\): the RMS values at \(\chi N=25,26,27\) are 0.01565,
  0.01865, and 0.01568 at \(n_x=128\);
- \(f_A=0.25\): the corresponding values at \(\chi N=26,27,28\) are
  0.00787, 0.00978, and 0.00890;
- \(f_A=0.30\): the shallow reversal at \(\chi N=28,29,30\) remains 0.01410,
  0.01401, and 0.01433.

The strongest near-spinodal structure is also common to independent OPF,
BURP, and BVK2 calculations. For \(f_A=0.20\), the adjacent SCFT-profile RMS
change is 0.0565 from \(\chi N=25\) to 26, then 0.0278 from 26 to 27, consistent
with rapid finite-amplitude evolution immediately above the RPA spinodal at
\(\chi N=24.613\). A profile error is a distance between two evolving
solutions and is not required to vary monotonically with segregation.

## Final decision

The nested direct-density \(n_x=256\) check contracts for the periods but not
for the profile RMS at selected \(f_A=0.25\) and 0.30 states, even though every
root remains converged. Six filtered, oversampled production-formulation tests
remove the prominent \(f_A=0.25\) peak. The complete 146-state fixed-stiffness
map was therefore recomputed at \(n_x=256\) with twofold oversampling and the
same physical sensor filter used by BVK2. All 146 roots pass the field,
composition, primitive-lamella, cell-stress, stress-orientation, and local
minimum gates. The publication curves now use
`../bvk_fixed_stiffness_period_map_nx256/summary.csv` exclusively.

The recomputation shows that the prominent \(f_A=0.25\) reversal was a
discretization/alignment artifact. The broad \(f_A=0.20\) maximum persists,
as does only a shallow variation at \(f_A=0.30\); these remaining features are
resolved model-to-SCFT behavior. The campaign is accepted and stopped; no
further numerical successor is justified.
