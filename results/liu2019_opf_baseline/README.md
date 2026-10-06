# Liu et al. 2019 OPF baseline

**Historical diagnostic.** The stored tables below predate the complete
[2019 correction to OPF Table 2](https://doi.org/10.1021/acs.macromol.9b01332).
They are retained for reference and are not the OPF data used in the manuscript.
The current source includes all four corrected entries. Rerunning this
global-regression diagnostic will therefore give different predictions.

The manuscript's OPF benchmarks use direct per-state force-and-stress fits
from eq 18, independently of Table 2. All 146 weak-to-intermediate and 105
strong-segregation OPF conditions are unaffected by the correction. The OK
comparisons inherit only unchanged B3/B4 entries and are also unaffected.
No benchmark tables or figures have been changed.

Run the focused coefficient checks from the repository root with
`julia --project=. test/opf_coefficient_correction_tests.jl`.

This diagnostic evaluated the optimized phase-field (OPF)
regressions of Liu et al., Macromolecules 52, 2878--2888 (2019),
against the existing stress-free Polyorder SCFT lamellar references.
Each OPF row is relaxed in its own optimized cell; no SCFT period is
imposed on the OPF profile.

This is an auxiliary implementation check. For the comparisons in Figures 3–5,
start with the [main AB benchmark](../liu2019_stress_free_period_map/README.md).
The implementation is in [`src/opf.jl`](../../src/opf.jl).

| fA | chiN | OPF L/Rg | SCFT L/Rg | signed period error | profile RMS | local minimum | accepted |
| ---: | ---: | ---: | ---: | ---: | ---: | --- | --- |
| 0.5 | 12 | 3.472606 | 3.426727 | +0.01339 | 0.003143 | true | true |
| 0.5 | 15 | 3.810655 | 3.713562 | +0.02615 | 0.0042378 | true | true |
| 0.5 | 20 | 4.170747 | 4.0457 | +0.03091 | 0.0041662 | true | true |
| 0.5 | 25 | 4.437186 | 4.288402 | +0.03469 | 0.0064435 | true | true |
| 0.5 | 30 | 4.649607 | 4.48235 | +0.03731 | 0.0075828 | true | true |
| 0.35 | 30 | 4.622977 | 4.442837 | +0.04055 | 0.036572 | true | true |

Representative mean absolute period error: `0.0304997`.
Representative mean profile RMS: `0.0103575`.

Acceptance requires solver convergence, an interior two-sided cell
minimum, mean-composition error below `1e-8`, and finite energy.
Period/profile errors are reported predictions, not acceptance gates.

The SCFT source tables are supplied in
[`diblock_stress_free_lamella_comparison/`](../diblock_stress_free_lamella_comparison/).
