# QCE control for AB lamellar comparisons

This directory supplies the fixed-stiffness QCE results for Figures 3–5.
Read these tables alongside the [multi-model AB benchmark](../liu2019_stress_free_period_map/README.md).
The strong-segregation QCE results are stored separately in
[`fixed_stiffness_summary.csv`](../bvk_strong_segregation_extension/fixed_stiffness_summary.csv).

The calculation imposes \(K_{\psi,2}=K_{\psi,0}\) while retaining the AQCE
grid (`nx=256`), twofold Fourier oversampling, and physical sensor filter
(`kc/kstar=12`). Profiles are cyclically aligned to the SCFT references on
1024 normalized-cell samples.

## Files and numerical checks

- [`summary.csv`](summary.csv): periods, errors, and numerical checks, with
  `model=bvk2_fixed` identifying QCE.
- [`profiles.csv`](profiles.csv): aligned density profiles.
- [`cases/`](cases/): individual root and profile files.
- [`validation_report.json`](validation_report.json): aggregate validation results.

All 146 states pass the field, composition, primitive-lamella, cell-stress, stress-orientation, and local-period-minimum gates. The largest projected-force RMS is `7.176668e-06`, the largest projected-force component is `1.954148e-05`, and the largest absolute cell stress is `6.400122e-07`.

- Full 146-state mean absolute period error: 3.643018%.
- Full 146-state mean profile RMS: 0.01282003.
- Common 133-state mean absolute period error: 3.616440%.
- Common 133-state mean profile RMS: 0.01310302.

The figures use the state selection defined by the plotting scripts. The
full-map and common-subset averages above refer to different selections.

## Reproduction

For figure generation, follow the [repository instructions](../../REPRODUCING.md).
To reassemble these tables from the supplied case files, run from the repository root:

```sh
uv run python scripts/assemble_bvk2_fixed_stiffness_period_map.py
```

The assembler checks the root/profile inventory, model settings, stationarity,
cell selection, and source hashes before replacing the summary tables. Use a
working copy if you want to preserve the distributed files unchanged.
