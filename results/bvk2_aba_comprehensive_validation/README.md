# Symmetric ABA validation (Figures 7 and 8)

This directory tests transfer of the frozen AQCE coefficient (`c2=0.16`) to
linear symmetric ABA lamellae. ABA SCFT data are validation references; they
were not used to refit the coefficient.

## Figures and data

- [macromolecules_figure7.svg](macromolecules_figure7.svg) compares ABA
  profiles at `chiN=20, 30, 45` and the signed period and profile-RMS errors for
  seven states from `chiN=20` to `45`. Each profile uses its model's optimized
  one-period cell. The models shown are SCFT, UD, BURP, QCE, and AQCE.
- [macromolecules_figure8.svg](macromolecules_figure8.svg) compares the SCFT
  and AQCE LAM/DIS bifurcations and the ordered-branch behavior above them.
- [benchmark_summary.csv](benchmark_summary.csv) contains stress-free periods,
  SCFT-relative errors, profile metrics, and acceptance fields.
- [benchmark_profiles.csv](benchmark_profiles.csv) contains aligned profiles
  on the common normalized coordinate.
- [scft_reference_summary.csv](scft_reference_summary.csv) and
  [scft_reference_profiles.csv](scft_reference_profiles.csv) retain the
  independently relaxed SCFT references.
- [ldis_crossing_summary.csv](ldis_crossing_summary.csv) stores the two analytic
  bifurcation coordinates; [ldis_crossing.csv](ldis_crossing.csv) stores the
  accepted ordered-branch points above the transition.
- [cell_stress_summary.csv](cell_stress_summary.csv) and
  [cell_stress_scan.csv](cell_stress_scan.csv) are supporting mechanism data,
  not Figure 7 inputs.

The QCE series in Figure 7 is read from the matched fixed-stiffness dataset in
[../bvk2_aba_fixed_stiffness_validation/](../bvk2_aba_fixed_stiffness_validation/).
Rows with `model=bvk1` describe the nonsmooth precursor and are not plotted in
Figure 7.

The symmetric ABA LAM/DIS transition is continuous. Accordingly, Figure 8 uses
the Gaussian/RPA bifurcation coordinates—`chiN=17.995811`, `D/Rg=2.385748`
for SCFT and `chiN=17.688844`, `D/Rg=2.386984` for AQCE—rather than a
transversal free-energy crossing. The finite-amplitude points above the
bifurcation demonstrate the expected decrease in ordered-branch free energy
and growth of the lamellar amplitude; they do not redefine the transition.

## Reproduce

From the repository root, redraw Figures 7 and 8 from the bundled tables:

```bash
uv run python scripts/render_macromolecules_figure8.py
```

Regenerating the underlying SCFT and reduced-model tables is separate and
requires the private Polyorder environment:

```bash
julia --project=. scripts/write_bvk2_aba_comprehensive_validation.jl --recompute=true
```

Omit `--recompute=true` only when using the generator's accepted cached
numerical data. See [REPRODUCING.md](../../REPRODUCING.md) for the full setup.
