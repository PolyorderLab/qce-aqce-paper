# Low-\(f_A\) SCFT Reference Promotion

## Objective

Promote the accepted (f_A=0.20, chi N=25,26,27) SCFT roots from the
resolution audit into the canonical lamellar benchmark. Refit and recompute
OPF at these states because its coefficients are mapped to the SCFT density
and stress. Preserve all other reduced-model solutions and recompute only
their SCFT-relative period and profile metrics.

## Acceptance criteria

- The staged SCFT rows and profiles must match the accepted audit SHA-256
  fingerprints and retain the (dx\simeq0.05R_g, ds=0.005) protocol.
- Each OPF state must pass the existing field, cell, resolution, composition,
  morphology, and local-minimum gates under statewise force-and-stress
  mapping.
- The staged three-case tables must have complete, unique model and profile
  inventories, and all stored errors must reproduce from native profiles.
- After promotion, the full 146-state map validator, fixed-stiffness
  assembler, figure renderers, and manuscript-facing regression tests must
  pass.

## Resources and stop condition

The three OPF states are independent and may run concurrently, with one Julia
thread and one BLAS thread per state. No reduced model other than OPF may be
recomputed. One changed retry is allowed for a rejected OPF state. The
campaign stops only after the canonical tables and dependent figures are
rebuilt successfully, or after a named failed gate identifies a blocker.

## Accepted outcome

The three refined SCFT references were promoted with periods 3.3772815752,
3.6724051029, and 3.8282917975 \(R_g\) at \(\chi N=25\), 26, and 27,
respectively. The statewise mapped OPF calculations pass every production
gate. An independent \(n_x=128\rightarrow256\) check gives maximum changes of
\(4.80\times10^{-7}\) in relative period, \(4.37\times10^{-5}\) in the
SCFT-referenced profile RMS, and \(8.35\times10^{-5}\) in the directly aligned
OPF profile. The full 146-state validator passes after promotion; no successor
calculation is justified.
