# Data guide

Use the index below to find each manuscript figure and the numerical tables
behind it. CSV files are plain text and can be opened with a spreadsheet,
Python, Julia, or another analysis tool. No simulation software is needed to
inspect them.

## Figures and source data

| Figure | Subject and supplied graphic | Where to find the data |
|---|---|---|
| 1 | [Quadratic response, spinodal, and weak-period comparison](results/bvk2_model_hierarchy_graphics/bvk2_model_hierarchy_weak_response.svg) | [Kernel samples](results/diblock_kernel_comparison/kernel_samples.csv) and [weak-response map](results/diblock_weak_response_map/weak_response_map.csv) |
| 2 | [Lamellar free energy and cell stress](results/lamellar_cell_stress_mechanism/lamellar_cell_stress_mechanism.svg) | [Scanned curves](results/lamellar_cell_stress_mechanism/scan.csv) and [equilibrium-period summary](results/lamellar_cell_stress_mechanism/summary.csv) |
| 3 | [AB lamellar period errors](results/liu2019_stress_free_period_map/liu2019_stress_free_period_comparison.svg) | [AB benchmark](results/liu2019_stress_free_period_map/README.md), [QCE control](results/bvk_fixed_stiffness_period_map_nx256/), and [strong-segregation extension](results/bvk_strong_segregation_extension/) |
| 4 | [AB density profiles and profile errors](results/liu2019_stress_free_period_map/liu2019_stress_free_profile_comparison.svg) | `summary.csv` and `profiles.csv` in the AB benchmark and QCE-control directories |
| 5 | [Aggregate period and profile errors](results/liu2019_stress_free_period_map/aggregate_error_bars.svg) | The AB benchmark and QCE-control summaries, plus [strong-segregation aggregates](results/bvk_strong_segregation_extension/aggregate.csv) |
| 6 | [AB phase diagram](results/bvk2_publication_phase_diagram/phase_diagram.svg) | [AQCE plot data](results/bvk2_publication_phase_diagram/phase_diagram_plot_data.csv), [SCFT reference curves](results/bvk2_publication_phase_diagram/scft_reference_boundaries.csv), and [phase-diagram notes](results/bvk2_publication_phase_diagram/phase_diagram.md) |
| 7 | [Symmetric ABA profiles and errors](results/bvk2_aba_comprehensive_validation/macromolecules_figure7.svg) | [Benchmark summary](results/bvk2_aba_comprehensive_validation/benchmark_summary.csv), [profiles](results/bvk2_aba_comprehensive_validation/benchmark_profiles.csv), and [QCE ABA control](results/bvk2_aba_fixed_stiffness_validation/README.md) |
| 8 | [ABA lamellar–disordered bifurcation](results/bvk2_aba_comprehensive_validation/macromolecules_figure8.svg) | [Ordered-branch data](results/bvk2_aba_comprehensive_validation/ldis_crossing.csv) and [bifurcation summary](results/bvk2_aba_comprehensive_validation/ldis_crossing_summary.csv) |
| S1–S2 | [Homopolymer validation](results/homopolymer_urp_vk1_validation/README.md) | Source plots and numerical inputs in that directory's `source/` folder |
| TOC | [Phase behavior and architecture transfer](results/bvk2_model_hierarchy_graphics/bvk2_macromolecules_toc.svg) | The same phase-boundary and ABA period tables used for Figures 6 and 7 |

QCE results are stored separately from the multi-model benchmark.
Read the QCE-control tables alongside the main summaries when comparing all
models. The strong-segregation directory likewise separates
[`summary.csv`](results/bvk_strong_segregation_extension/summary.csv) from
[`fixed_stiffness_summary.csv`](results/bvk_strong_segregation_extension/fixed_stiffness_summary.csv).
The plotting scripts combine these inputs and apply the figure's state selection.

## Reading the tables

Column names vary slightly between datasets. The most common fields are:

| Field | Meaning |
|---|---|
| `f` or `fA` | Overall A-segment fraction, $f_A$ |
| `chiN` | Segregation parameter, $\chi N$ |
| `case_id` | State identifier containing composition and segregation strength |
| `model` | Machine-readable model identifier; see the key below |
| `period_rg`, `scft_period_rg` | Model and SCFT periods in units of $R_g$ |
| `signed_period_error` or `period_relative_error` | $(D-D_{\mathrm{SCFT}})/D_{\mathrm{SCFT}}$, stored as a fraction |
| `abs_period_error` | Absolute value of the fractional period error |
| `profile_rms` | Root-mean-square difference between the aligned model and SCFT density profiles |
| `s` | Position within a period, $x/D$ |
| `x_rg` | Position in units of $R_g$, where provided |
| `phi_a`, `phi_b` | A- and B-segment volume fractions |
| `status`, `accepted`, or `*_gate_pass` | Numerical validation status or individual checks from the calculation |

Multiply a fractional period error by 100 to obtain percent. Columns explicitly
named `*_percent`, such as `mean_abs_period_error_percent`, are already in
percent. Profile RMS is a dimensionless density difference, not a relative
percentage error.

Join summary and profile tables using composition, `chiN`, and `model`.
`case_id` strings can differ in decimal formatting between datasets, so do not
assume that identifiers from different directories have identical spelling.
Each model's profile is evaluated in its own periodic cell; use `s` to compare
profile shapes and the period columns to compare domain spacing. Empty fields
and `NaN` entries should not be interpreted as zero.

## Model names

Use this key to identify the models in tables and code:

| Identifier | Model |
|---|---|
| `scft` | Self-consistent field theory (SCFT), the quantitative reference |
| `bvk2_fixed` | QCE, the fixed-stiffness density functional |
| `bvk2` | AQCE, the adaptive-stiffness density functional |
| `uneyama_doi` | Uneyama–Doi (UD) model |
| `liu2019_opf` | Optimized phase-field (OPF) model |
| `burp_ti` | BURP, the full-RPA comparison model |
| `ohta_kawasaki` | Ohta–Kawasaki (OK) model |
| `bvk1` | Nonsmooth adaptive precursor in auxiliary results |
| `bvk1_bvk2` | Shared reduced-response entry in the quadratic-response tables |

For QCE/AQCE calculations, check the `model` and `adaptive` fields rather than
inferring the model from the directory name alone.

## Interpreting phase-boundary data

For the plotted AQCE boundaries, start with
[`phase_diagram_plot_data.csv`](results/bvk2_publication_phase_diagram/phase_diagram_plot_data.csv).
The SCFT reference curves are separate literature-derived data, not AQCE
calculations. Their source is documented in the
[reference-extraction notes](results/jcp2021_ab_phase_diagram_digitization/README.md).

[`accepted_roots/`](accepted_roots/) retains the numerical boundary records and
supporting evidence. A numerical bracket width measures the interval used to
locate a boundary; it is not by itself a total uncertainty including spatial
resolution. The [resolution audit](results/bvk2_resolution_uncertainty_table/README.md)
explains the distinction and identifies unresolved checks.

Use the linked tables above as plotting inputs. Paths in `source` columns
identify the underlying calculations and may refer to files outside this
repository.
