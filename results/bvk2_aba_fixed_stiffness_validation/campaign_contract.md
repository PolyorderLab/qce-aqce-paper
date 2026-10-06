# Fixed-Stiffness ABA Figure 6 Campaign Contract

- **Objective:** compute the fixed-stiffness BVK control for the seven
  symmetric ABA lamellar states plotted in Figure 6(a--c).
- **Claim kind:** stress-free lamellar observables (period and aligned profile),
  one independently optimized state per \(\chi N=20,22,25,30,35,40,45\).
- **Accepted outputs:** one `root.csv` and one `profile.csv` under each
  `cases/f0.50_chiN*/` directory; after all states pass, assemble the canonical
  fixed-stiffness `summary.csv` and `profiles.csv` used by Figure 6.
- **Functional:** the ABA BVK functional with `c2=0.16` and
  `adaptive=false`, so \(K_{\psi,2}=K_{\psi,0}\). The Gaussian coefficients,
  nonlinear local terms, grid (`nx=128`), and field/cell workflow match the
  accepted adaptive ABA benchmark.
- **Acceptance:** converged field; projected-force RMS at most `1e-4` and
  maximum at most `2e-4`; mean composition error at most `2e-8`; nonuniform
  one-domain lamella with dominant Fourier mode 1; finite interior period
  minimum with accepted two-sided local check; finite aligned SCFT profile.
- **Resources:** seven independent one-thread Julia processes, at most seven
  concurrent jobs, no BLAS oversubscription, and no additional SCFT or
  adaptive-BVK calculations.
- **Retry budget:** one numerically changed retry per state after diagnosis.
- **Stop:** all seven states are accepted, the aggregate sources pass their
  inventory/provenance contract, Figure 6 is regenerated deterministically,
  and manuscript/SI captions and tests agree with the plotted model set.
