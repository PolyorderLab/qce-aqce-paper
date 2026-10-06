# Fixed-stiffness QCE control for symmetric ABA lamellae

This dataset supplies the QCE control in **Figure 7**.

The calculation uses the same symmetric ABA coefficients, nonlinear local
terms, `c2=0.16`, `nx=128`, and stress-free-period workflow as the AQCE
benchmark, while setting `adaptive=false` so that
\(K_{\psi,2}=K_{\psi,0}\). All seven states pass the acceptance checks recorded
in [campaign_contract.md](campaign_contract.md). Their mean absolute period
error is 3.83%, and their mean aligned profile RMS is 0.02267.

## Files

- [summary.csv](summary.csv) contains the seven accepted stress-free roots and
  SCFT-relative period/profile errors.
- [profiles.csv](profiles.csv) contains 256 aligned samples per state.
- `cases/<state>/root.csv` and `cases/<state>/profile.csv` retain the accepted
  case-level inputs used by the assembler.

From the repository root, rebuild the aggregate CSV files from the bundled
case data:

```bash
uv run python scripts/assemble_bvk2_aba_fixed_stiffness_validation.py
```

Figure 7 itself is then rendered from this dataset and the ABA benchmark with:

```bash
uv run python scripts/render_macromolecules_figure8.py
```

Rerunning the numerical cases is separate and requires the private Polyorder
environment; the case entry point is
`scripts/run_bvk2_aba_fixed_stiffness_case.jl`. See
[REPRODUCING.md](../../REPRODUCING.md) for the supported environment and full
workflow.
