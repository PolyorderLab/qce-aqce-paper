# AQCE convergence and morphology diagnostics

These auxiliary results document representative LAM, CYL, BCC, and GYR
calculations using the [phase-crossing endpoints](../bvk2_boundary_crossing_panels/crossing_endpoints.csv).
They support the numerical checks rather than form part of the current two-figure
Supporting Information. The supplied tables and field sections can be plotted
without rerunning the simulations.

| Phase | Accepted endpoint | Grid | Cell stress | Shell | Stop reason |
|---|---|---:|---:|---:|---|
| BCC | S/C at `chiN=47.5, fA=0.2085` | 48x48x48 | 0.000115211 | 110 | lbfgs_iteration_cap |
| CYL | G/C at `chiN=54, fA=0.3630` | 96x168 | 0.00119227 | 11 | lbfgs_converged |
| GYR | G/C at `chiN=54, fA=0.3630` | 112x112x112 | -0.00026751 | 211 | iteration_cap |
| LAM | G/L at `chiN=54, fA=0.3630` | 1024 | 0.00024709 | 1 | lbfgs_iteration_cap |

## Selection contract

- GYR, CYL, and LAM are the accepted endpoint states at the shared
  `chiN=54, fA=0.363` coordinate. The common-state comparison is already used
  to certify the nonempty GYR pocket.
- BCC is the accepted `chiN=47.5, fA=0.2085` S/C endpoint from the
  representative bracket where all endpoint stresses satisfy the `1e-3`
  target.
- Every selected row has `endpoint_status=accepted`,
  `phase_status=accepted`, and a passing morphology-shell audit. The script
  verifies the canonical endpoint source SHA-256 before reading the artifact.

## Files

- `representative_states.csv`: endpoint energy, stress, shell, plateau,
  residual, optimizer-stop, and content-hash provenance.
- `convergence_traces.csv`: raw objective trace normalized by cell volume plus
  the selected plateau record. Projected-theta norms are explicitly marked
  `diagnostic_only`.
- `morphology_sections.csv`: the complete 1-D LAM profile, complete 2-D CYL
  cell, and deterministic maximum-variance `xy` sections for BCC and GYR.
- `bvk2_si_convergence_morphology.svg`: publication-oriented convergence and
  morphology overview.
- `bvk2_si_convergence_morphology.png`: inspection preview rendered from the
  same tables. The SVG remains the authoritative publication asset.

## Interpretation and limitations

Acceptance is a boundary/phase-endpoint statement, not the value of the
optimizer's `converged` flag. Some exact-theta states stop at a bounded
iteration cap after an accepted energy plateau. The projected-theta R2/Rinf
values are therefore shown as diagnostics and are not silently converted into
a universal gate.

The canonical energy can use a higher post-solve quadrature factor than the
optimization trace. `accepted_energy_density` is the canonical endpoint value;
`trace_objective_energy_density` is the objective represented by the
convergence trace. Their difference is retained rather than hidden.

The BCC and GYR panels are scalar-field sections selected by maximum in-plane
variance (BCC `k=1`, GYR
`k=40`). They illustrate the accepted density
fields but do not by themselves prove 3-D topology; the reciprocal-shell
fingerprints (`110` and `211`) provide the phase-identity audit.

## Reproduction

To plot the supplied diagnostic tables, run from the repository root after
setting up the [plotting environment](../../REPRODUCING.md):

```bash
uv run python -c 'from scripts.render_macromolecules_figures import render_convergence_morphology; render_convergence_morphology()'
```

No simulation is launched. Schema: `bvk2-si-convergence-morphology-v1`.
