# Lamellar cell-stress mechanism

This publication artifact shows why stress-free period is a stronger discriminator than a visually similar density profile. Figure 2 plots UD, BURP, QCE, and AQCE, while the validated source ledger also retains OPF for the broader comparison in Section 3.3. At each fixed period, the density is independently relaxed in the one-period LAM basin. The plotted energy is the relaxed free-energy density `F/L` after subtracting the accepted own-cell value and dividing by `max(abs((F/L)_root), 1)`. The cell stress is the common centered derivative

```math
\sigma_L = \frac{d(F/L)}{d\ln L},
```

evaluated at fixed relaxed nodal density. Within the source ledger, OPF, QCE, and AQCE use their analytic isotropic stress; UD and BURP use a centered derivative with `delta ln L = 0.0002`. By the envelope theorem this frozen-field derivative is the derivative of the minimized branch when the field stationarity gate passes. Normalization is within each model, so neither panel ranks absolute free energies across different functionals.

The black dashed marker is the independently cell-solved Polyorder SCFT period. Filled symbols mark each density functional's accepted own-cell period.

Branch stationarity is checked against each solver's accepted residual convention: UD and BURP use RMS/max tolerances `5e-4/2.5e-3`, mapped OPF uses `1e-5/1e-4`, and QCE/AQCE use `1e-5/2e-5`. Every source profile must also contain exactly one periodic A-rich domain.

## Validated source roots

| state | model | `L/L_SCFT` | period error | normalized root stress | internal minimum |
| --- | --- | ---: | ---: | ---: | --- |
| `(0.50, 20)` | Uneyama-Doi | 0.865281 | -13.472% | 1.453e-03 | pass |
| `(0.50, 20)` | OPF | 1.004187 | +0.419% | 1.339e-07 | pass |
| `(0.50, 20)` | BURP | 0.819145 | -18.085% | 7.345e-04 | pass |
| `(0.50, 20)` | QCE | 0.973450 | -2.655% | -5.276e-09 | pass |
| `(0.50, 20)` | AQCE | 1.001937 | +0.194% | -1.005e-09 | pass |
| `(0.35, 30)` | Uneyama-Doi | 0.857722 | -14.228% | 3.066e-03 | pass |
| `(0.35, 30)` | OPF | 0.969712 | -3.029% | 3.204e-09 | pass |
| `(0.35, 30)` | BURP | 0.822203 | -17.780% | 1.865e-03 | pass |
| `(0.35, 30)` | QCE | 0.957336 | -4.266% | -7.876e-09 | pass |
| `(0.35, 30)` | AQCE | 0.992891 | -0.711% | -6.664e-14 | pass |

## Interpretation

- The matched QCE control differs from AQCE only by setting `adaptive=false`, so that `K_{psi,2}=K_{psi,0}` while the Gaussian kernel, nonlinear local terms, grid, oversampling, and physical sensor filter remain unchanged.
- Activating the adaptive stiffness moves the cell-stress zero toward the SCFT marker at both the held-out symmetric state and the asymmetric training state.
- BURP and Uneyama-Doi retain valid interior minima, but their zero-stress cells are systematically too short. Their period error is therefore a constitutive cell-stress error, not a failed period minimizer.

## Reproduction

```bash
julia --project=. scripts/write_lamellar_cell_stress_mechanism.jl
```

Inputs: `results/unified_lamellar_benchmark/summary.csv` and `results/unified_lamellar_benchmark/profiles.csv`. Statewise OPF coefficients and native profiles are read from the accepted Liu-2019 stress-free-period map. Both QCE variants use frozen `c2=0.16`, `nx=256`, oversampling factor 2, and sensor-filter ratio 12.0; only the adaptive flag changes. The SCFT markers are reused only from rows with successful field and cell gates in the unified benchmark.

The Julia generator validates and writes `scan.csv`, `summary.csv`, and `README.md`. The publication SVG is written exclusively by `scripts/render_macromolecules_figures.py::render_cell_stress` after those data products pass all gates.
