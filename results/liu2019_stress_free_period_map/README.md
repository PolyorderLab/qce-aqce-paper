# Stress-free lamellar benchmark (Figures 3–5)

This directory contains the primary weak- and intermediate-segregation data
used to compare stress-free lamellar periods and profiles for SCFT, OPF,
Ohta–Kawasaki, Uneyama–Doi, BURP, and AQCE. The states cover the published
Liu et al. composition–segregation domain; `fA=0.15` is absent because its
spinodal lies above the plotted `chiN` range.

## Figures and data

- [liu2019_stress_free_period_comparison.svg](liu2019_stress_free_period_comparison.svg)
  is Figure 3, the stress-free period comparison.
- [liu2019_stress_free_profile_comparison.svg](liu2019_stress_free_profile_comparison.svg)
  is Figure 4, the representative aligned-profile comparison.
- [aggregate_error_bars.svg](aggregate_error_bars.svg) is Figure 5, the period
  and profile error summary for the common set of states and the strong-segregation
  extension.
- [summary.csv](summary.csv) contains one row per state and model, including
  the optimized period, SCFT reference, errors, solver diagnostics, and status.
- [profiles.csv](profiles.csv) contains phase-aligned profiles on a normalized
  periodic coordinate. These are comparison profiles, not native solver grids.
- [validation_report.json](validation_report.json) records validation totals.

Figures 3–5 also consume the bundled QCE control in
[../bvk_fixed_stiffness_period_map_nx256/summary.csv](../bvk_fixed_stiffness_period_map_nx256/summary.csv)
and
[../bvk_fixed_stiffness_period_map_nx256/profiles.csv](../bvk_fixed_stiffness_period_map_nx256/profiles.csv).
The strong-segregation portion of Figures 3 and 5 comes from
[../bvk_strong_segregation_extension/summary.csv](../bvk_strong_segregation_extension/summary.csv),
[../bvk_strong_segregation_extension/fixed_stiffness_summary.csv](../bvk_strong_segregation_extension/fixed_stiffness_summary.csv),
and [../bvk_strong_segregation_extension/aggregate.csv](../bvk_strong_segregation_extension/aggregate.csv).

All accepted rows passed the recorded field, composition, morphology, cell,
and local-period checks. The AQCE sensor/filter choice is validated for this
lamellar comparison but should not be read as an absolute-energy validation
between different morphologies.

## Reproduce

From the repository root, redraw these figures from the bundled CSV files:

```bash
uv run python -c 'from scripts.render_macromolecules_figures import render_liu2019_stress_free_period_map, render_lamellar_aggregate_errors; render_liu2019_stress_free_period_map(); render_lamellar_aggregate_errors()'
```

Recomputing the SCFT references and numerical tables is a separate operation
that requires the private Polyorder environment:

```bash
julia --project=. scripts/write_liu2019_stress_free_period_map.jl
```

See [REPRODUCING.md](../../REPRODUCING.md) for environment setup, dependency
limits, and the complete figure workflow.
