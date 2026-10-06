# AQCE boundary free-energy crossings

These tables provide representative numerical evidence for the AQCE phase
boundaries. Each example compares phase free energies at the ends of a
fixed-`chiN` bracket. The supplied tables can be plotted without rerunning
field calculations.

| Boundary | χN | Bracket | Root | Conservative uncertainty | Method |
|---|---:|---:|---:|---:|---|
| S/DIS | 14 | [0.301782173, 0.302407173] | 0.302094673 | ±0.0003125 | survival_bracket_midpoint |
| S/C | 47.5 | [0.207500000, 0.208500000] | 0.207827655 | ±0.0005000 | signed_endpoint_secant |
| G/C | 54 | [0.362000000, 0.363000000] | 0.362421154 | ±0.0005000 | signed_endpoint_secant |
| G/L | 54 | [0.363000000, 0.364000000] | 0.363343478 | ±0.0005000 | signed_endpoint_secant |

`crossing_endpoints.csv` records the two accepted phase energies, their
difference, stress, morphology shell, grids, audit factor, source artifact, and
SHA-256 provenance for every panel. `crossing_summary.csv` records the bracket,
root, uncertainty, secant diagnostic, and canonical publication-ledger source.

## Interpretation

- S/C, G/C, and G/L are ordinary strict signed free-energy crossings. Their
  plotted roots are endpoint secants and the reported uncertainty is half the
  accepted bracket width.
- S/DIS is an order-disorder survival/bifurcation edge. Its canonical value is
  the accepted bracket midpoint, not the nearby secant diagnostic. The dashed
  connector in the figure makes this different semantics explicit.
- At `chiN=54`, `fA=0.363`, the shared GYR energy is
  `-6.46152100806203`. It is below CYL by
  `0.000416211` and below LAM by
  `0.000620284` in energy-density
  units. This common-state comparison directly establishes a nonempty GYR
  pocket without over-interpreting the separation between two interpolated
  roots.

The selected S/C slice is `chiN=47.5`: all four endpoint cell stresses are
inside the `1e-3` target and the independent factor-4 energy evaluation gives
strict opposite signs. The G/C and G/L panels use direct GYR `112^3`, CYL
`96x168`, and LAM `1024` states with factor-3 energy evaluation.

## Reproduction

To redraw the supplied crossing data, run from the repository root after
setting up the [plotting environment](../../REPRODUCING.md):

```bash
uv run python -c 'from scripts.render_macromolecules_figures import render_crossings; render_crossings()'
```

Data and figure files:

- `bvk2_boundary_crossings.svg`
- `crossing_endpoints.csv`
- `crossing_summary.csv`

Schema: `bvk2-boundary-crossing-panels-v1`.
