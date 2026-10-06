# ABA LAM--DIS \(\chi N=21\) point contract

- **Objective:** add the matched SCFT and BVK2 ordered-lamellar branch point at
  \(f_A=0.50\), \(\chi N=21\) to the canonical Figure 7 source table.
- **Claim kind:** `ordered_branch_point`.
- **Accepted outputs:** `ldis_crossing.csv`, `ldis_crossing_summary.csv`, and
  `decision.csv` in this directory. Promotion into the canonical crossing table
  occurs only after acceptance.
- **Acceptance:** exactly one SCFT and one BVK2 row; finite free-energy
  difference, amplitude, and period; negative LAM--DIS free-energy difference;
  positive ordered-state amplitude; successful field and cell solves; SCFT
  cell-stress norm no greater than \(10^{-5}\).
- **Resources:** one local job, two Julia threads, no parameter sweep, and a
  30-minute wall-time target.
- **Retry budget:** one retry, only after diagnosis and with a changed numerical
  hypothesis.
- **Scientific status:** accepted; both model rows passed on the first attempt.
- **Stop condition:** both model rows are accepted, promoted without changing
  the analytical bifurcation coordinates, and Figure 7 plus its provenance
  tests are rebuilt successfully.
