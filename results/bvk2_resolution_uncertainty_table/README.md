# AQCE phase-boundary resolution and uncertainty audit

This audit describes the spatial resolution and bracket uncertainty of
**106 phase-boundary points**: 33 S/DIS, 27 S/C,
22 G/C, and 24 G/L. The coordinates come from
`results/bvk2_publication_phase_diagram/accepted_phase_boundaries.csv` and the
exact campaign/correction ledgers named in each row.

| Boundary / contract | Points | Resolutions | Energy treatment | Coordinate half-width | Audit status |
|---|---:|---|---|---:|---|
| `GC_FILTERED_R6_GYR96_CYL128x224` | 9 | GYR 96x96x96; CYL 128x224 | direct filtered-sensor bracket | 0.001 | current_unified_direct_contract |
| `GC_RICHARDSON_64x64x64_TO_96x96x96_48x84_TO_128x224` | 8 | GYR 64x64x64 -> 96x96x96; CYL 48x84 -> 128x224 | second-order Richardson on both phases | 0.000458–0.000984 | resolution_shift_exceeds_stored_half_width |
| `GC_THETA_GYR112_CYL96x168_AUDIT3` | 1 | GYR 112x112x112; CYL 96x168 | direct; no Richardson correction | 0.0005 | current_unified_direct_contract |
| `GC_THETA_GYR96_CYL48x84_AUDIT3` | 4 | GYR 96x96x96; CYL 48x84 | direct; no Richardson correction | 0.0005 | legacy_accepted_endpoint_exceeds_current_stress_maximum;legacy_direct_current_stress_gate_pass |
| `GL_MATCHED_GYR112_LAM46_AUDIT3` | 7 | GYR 112x112x112; LAM 46 | matched-spacing direct factor-3 bracket | 0.00025–0.0005 | current_unified_matched_spacing_contract |
| `GL_MATCHED_GYR64_LAM26_AUDIT3` | 2 | GYR 64x64x64; LAM 26 | matched-spacing direct factor-3 bracket | 0.00075 | current_unified_matched_spacing_contract |
| `GL_MATCHED_GYR96_LAM40_AUDIT3` | 11 | GYR 96x96x96; LAM 40 | matched-spacing direct factor-3 bracket | 0.000476–0.00089 | current_unified_matched_spacing_contract |
| `GL_THETA_GYR112_LAM1024_AUDIT3` | 3 | GYR 112x112x112; LAM 1024 | direct; no Richardson correction | 0.0005 | current_direct_contract |
| `GL_THETA_GYR64_LAM1024_AUDIT3` | 1 | GYR 64x64x64; LAM 1024 | direct; no Richardson correction | 0.000312 | current_direct_contract |
| `SC_FIXED_FA_BCC48_CYL48x84` | 6 | BCC 48x48x48; CYL 48x84 | direct grid; no fine-grid correction | 0.0075–0.0226 | legacy_phase_specific_acceptance_not_uniformly_recertified |
| `SC_RICHARDSON_48x48x48_TO_128x128x128_48x84_TO_128x224` | 1 | BCC 48x48x48 -> 128x128x128; CYL 48x84 -> 128x224 | second-order Richardson on both phases | 0.000821 | resolution_shift_exceeds_stored_half_width |
| `SC_RICHARDSON_48x48x48_TO_64x64x64_48x84_TO_128x224` | 10 | BCC 48x48x48 -> 64x64x64; CYL 48x84 -> 128x224 | second-order Richardson on both phases | 0.000399–0.00096 | resolution_shift_exceeds_stored_half_width;resolution_shift_within_stored_half_width_but_total_bound_missing |
| `SC_THETA_BCC48_CYL48x84_AUDIT4` | 10 | BCC 48x48x48; CYL 48x84 | direct; no Richardson correction | 0.0005 | direct_root_current_stress_gate_pass;legacy_accepted_endpoint_exceeds_current_stress_maximum |
| `SDIS_BCC48_ROOT_CHIN_REFINED` | 9 | BCC 48x48x48; DIS analytic homogeneous | direct grid; no fine-grid correction | 0.00375 | accepted_survival_root_resolution_bias_unquantified |
| `SDIS_BCC48_ROOT_FA` | 24 | BCC 48x48x48; DIS analytic homogeneous | direct grid; no fine-grid correction | 0.000312–0.00043 | accepted_survival_root_resolution_bias_unquantified |

## What the uncertainty column means

`coordinate_half_width = bracket_width/2`. For fixed-`chiN` roots it is in
`fA`; for the 9 near-critical fixed-`fA` S/DIS roots it is in `chiN`.
It is **not** automatically a total continuum uncertainty:

- S/DIS is a survival/bifurcation bracket. Its half-width excludes BCC `48^3`
  spatial bias.
- Direct S/C, G/C, and G/L rows use signed endpoint brackets, but their
  half-widths exclude remaining spatial-grid bias.
- The 19 legacy Richardson-corrected S/C and G/C rows retain the raw numerical
  bracket width while their plotted coordinate is shifted by a continuum
  extrapolation. No uncertainty of that extrapolation is tabulated.

Of the 19 corrected rows (11 S/C and
8 G/C), **14** have a resolution
shift larger than the stored half-width: 8 G/C and
6 S/C rows. These rows therefore do not possess a strict
total coordinate bound under the current evidence.

## Unified-gate audit

The canonical plot remains unchanged by this read-only audit. The original
campaigns used phase-specific acceptance criteria. Applying the later uniform
stress audit retrospectively identifies **10** direct
S/C/G/C rows with at least one endpoint above `2e-3`. This comparison does not
by itself invalidate their original phase-specific acceptance, but it prevents
those rows from being described as uniformly recertified under the later audit.
They are identified point-by-point in
`phase_boundary_resolution_uncertainty.csv` with
`audit_status=legacy_accepted_endpoint_exceeds_current_stress_maximum`.

All 24 G/L rows use direct factor-3 energies. The matched-spacing
subset uses the GYR and LAM grids recorded point by point; the remaining
direct rows use LAM `1024`. Their audited endpoints remain below `2e-3`.
Filtered G/C uses GYR `96^3` and CYL `128x224`; G/C at `chiN=54` uses the
direct GYR `112^3` / CYL `96x168` contract.

For reuse of these boundary coordinates, two numerical questions remain open:
total resolution uncertainty for the legacy corrected S/C/G/C rows, and
recertification of the flagged direct S/C/G/C rows under the uniform stress
criterion. The pointwise records identify which coordinates each issue affects.

## Files

- `phase_boundary_resolution_uncertainty.csv`: all 106 canonical points.
- `resolution_contracts.csv`: grouped resolution/solver contracts.
- `resolution_contract_map.svg`: pointwise and contract-range visualization.
- `acceptance_rejection_audit.svg`: accepted-contract and rejected-row audit.

All `canonical_source` and `evidence_sources` entries are repository-relative;
the generator resolves both these portable paths and legacy absolute paths.

To reassemble the audit from its source ledgers, run from the repository root:

```bash
uv run python scripts/write_bvk2_resolution_uncertainty_table.py
uv run python test/test_write_bvk2_resolution_uncertainty_table.py
```

Schema: `bvk2-resolution-uncertainty-table-v1`.
