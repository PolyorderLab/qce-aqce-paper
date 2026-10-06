# AQCE phase diagram (Figure 6)

[phase_diagram.svg](phase_diagram.svg) is Figure 6. It shows the
accepted AQCE phase-boundary curves and the dashed SCFT literature reference
from Li and Liu, *J. Chem. Phys.* **154**, 014903 (2021). The computed
`fA <= 0.5` half is reflected by A/B symmetry.

Solid curves show AQCE boundaries, dashed curves show SCFT references, and
gray disks mark the FCC/BCC/DIS, G/C/O70, and G/O70/L triple-point estimates.

## Data files

- [phase_diagram_plot_data.csv](phase_diagram_plot_data.csv) is the canonical
  stable-topology table consumed by the renderer. It contains accepted roots
  only and omits pairwise continuations outside displayed stable pockets.
- [accepted_phase_boundaries.csv](accepted_phase_boundaries.csv) contains all
  accepted boundary claims; [boundary_candidates.csv](boundary_candidates.csv)
  also includes provisional claims.
- [boundary_diagnostics.csv](boundary_diagnostics.csv) retains diagnostic and
  metastable pairwise crossings that are not plotted as stable boundaries.
- [scft_reference_boundaries.csv](scft_reference_boundaries.csv) contains the
  digitized SCFT comparison curves.
- [fcc_region_zoom_triple_point.csv](fcc_region_zoom_triple_point.csv),
  [estimated_g_c_o70_triple_point.csv](estimated_g_c_o70_triple_point.csv), and
  [estimated_g_o70_l_triple_point.csv](estimated_g_o70_l_triple_point.csv)
  provide the displayed topology junctions. These estimates connect accepted
  curves for display; they are not additional accepted free-energy roots.
- [rejected_points.csv](rejected_points.csv) records rejected rows and is never
  plotted.
- [rpa_stability_limit.csv](rpa_stability_limit.csv) gives the homogeneous-state
  RPA stability limit as additional reference data.

The SCFT curves are a literature overlay, not AQCE calculations. Display-only
smoothing is applied to selected SCFT pocket curves to reduce digitization
stair steps while retaining endpoints and source vertices. Use the tables
linked above to reproduce the plot.

## Reproduce

From the repository root, redraw Figure 6 from the supplied plot table:

```bash
uv run python scripts/render_bvk2_phase_diagram_snapshot.py \
  --output results/bvk2_publication_phase_diagram/phase_diagram.svg
```

To calculate new boundaries, you will need additional raw calculation files
and the numerical environment described in [REPRODUCING.md](../../REPRODUCING.md).
These are not needed to plot the supplied boundary data.
