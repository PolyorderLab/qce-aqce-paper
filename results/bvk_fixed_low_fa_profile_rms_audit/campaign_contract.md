# Fixed-Stiffness BVK Low-Composition Profile-RMS Audit

## Objective and claim scope

This campaign tests whether the nonmonotonic profile-RMS features in the
fixed-stiffness BVK panel at (f_A=0.20), 0.25, and 0.30 arise from the
(n_x=64) diagnostic discretization or from the model-to-SCFT comparison.
The claim kind is a resolution audit of the stress-free lamellar observable;
the results are diagnostic until every gate below passes and the canonical
figure source is explicitly promoted.

## Pilot states

The smallest discriminating set contains the three-point neighborhoods of the
largest interior reversals:

- (f_A=0.20), \(\chi N=25,26,27\);
- (f_A=0.25), \(\chi N=26,27,28\);
- (f_A=0.30), \(\chi N=28,29,30\).

Each point is recomputed independently at (n_x=128) with (c_2=0), which is
algebraically the fixed-stiffness limit of the diagnostic implementation.

## Accepted outputs and gates

Each point must produce one terminal `candidate_rows.csv` row with:

- `status=accepted` and all stationarity, phase-identity, cell-stationarity,
  convergence, and local-period-minimum gates true;
- a non-boundary-limited primitive lamellar state;
- projected-force RMS no larger than (10^{-5}), maximum force component no
  larger than (2\times10^{-5}), and absolute cell stress no larger than
  (10^{-6});
- finite period and profile RMS, with complete provenance for (c_2=0) and
  (n_x=128).

The pilot is accepted as resolution-stable only if the \(n_x=64\to128\)
change is smaller than 0.25% in period and 0.0025 in SCFT-aligned profile RMS
at every state. These tolerances are discriminating rather than production
promotion thresholds: a failure triggers a production-discretization audit,
not automatic replacement of the publication dataset.

## Resources, retries, and stop condition

- At most nine concurrent processes, one Julia thread and one BLAS thread per
  process; no nested outer parallelism.
- Maximum 10,000 field iterations and 30 period iterations per point.
- Retry budget: two attempts per state, and every retry must change a diagnosed
  numerical cause or reuse a valid terminal state.
- First classify all nine \(n_x=128\) points. If the sharp features persist,
  repeat only these nine states at \(n_x=256\) because an RMS change that is
  small in absolute terms can still matter to the ranking among the best
  models. Accept grid convergence when every \(n_x=128\to256\) change
  contracts relative to its \(n_x=64\to128\) change, remains below 0.08% in
  period and 0.00075 in profile RMS, and all solver gates pass. Escalate to
  the filtered production fixed-stiffness formulation only if this nested
  confirmation fails.

Scientific artifact status is restricted to `accepted`, `provisional`, or
`rejected`; orchestration decisions are recorded separately by the campaign
controller.

## Pilot decision

All nine jobs completed successfully and were recorded as `accept`. The
largest period change was 0.156%, and the largest profile-RMS change was
0.00125. The local extrema persisted at \(n_x=128\). Because that RMS change
is comparable with the separation between the best aggregate models, a
nested \(n_x=256\) confirmation over the same nine states is the smallest
justified extension.

## Production-discretization decision

The nested direct-density audit contracts for all periods and for the
\(f_A=0.20\) profile RMS values, but not for the \(f_A=0.25\) and 0.30
profile RMS values. Six filtered, oversampled fixed-stiffness calculations at
those affected states pass every production gate. They remove the spurious
\(f_A=0.25\) peak while retaining only the much smaller \(f_A=0.30\) reversal.
Because fixed-stiffness BVK is now displayed as a complete panel and enters
the aggregate comparison, mixing these production results with the older
\(n_x=64\) diagnostic is not acceptable. The justified successor is therefore
the complete 146-state fixed-stiffness map with the same \(n_x=256\),
twofold-oversampled production formulation. Up to 96 independent one-thread
processes may run concurrently on the 112-core host. Promotion requires all
146 roots to pass the production field, composition, primitive-lamella,
cell-stress, stress-orientation, and local-minimum gates, followed by
consistent 1024-point cyclic profile alignment against SCFT.

## Final campaign decision

All 146 production-formulation roots passed the promotion gates. The accepted
map is assembled in `../bvk_fixed_stiffness_period_map_nx256/`, with
`summary.csv`, `profiles.csv`, and `validation_report.json` providing the
publication source and its validation record. The fixed-stiffness panel and
aggregate bars were regenerated from that source. The decision is `accept`;
the stop condition is satisfied and no successor is launched.
