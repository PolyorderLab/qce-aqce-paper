# Model-hierarchy and table-of-contents graphics

This directory contains Figure 1 and the table-of-contents graphic. Both can
be redrawn from the supplied tables without running simulations.

## Outputs

- [bvk2_model_hierarchy_weak_response.svg](bvk2_model_hierarchy_weak_response.svg)
  is Figure 1. It compares the homogeneous second-order composition vertices,
  spinodal errors, and weak-period errors of the model hierarchy.
- [bvk2_model_hierarchy_weak_response.png](bvk2_model_hierarchy_weak_response.png)
  is its inspection preview.
- [bvk2_macromolecules_toc.svg](bvk2_macromolecules_toc.svg) is the
  3.25 in × 1.75 in TOC graphic. Its left panel uses the Figure 6 AQCE and SCFT
  phase boundaries; its right panel uses the Figure 7 ABA period errors for UD,
  BURP, QCE, and AQCE.
- [bvk2_macromolecules_toc.png](bvk2_macromolecules_toc.png) is the 975 × 525
  inspection preview. Use the SVG for a scalable version.
- [source_manifest.csv](source_manifest.csv) lists the input tables.

Figure 1 reads
[../diblock_kernel_comparison/kernel_samples.csv](../diblock_kernel_comparison/kernel_samples.csv)
and
[../diblock_weak_response_map/weak_response_map.csv](../diblock_weak_response_map/weak_response_map.csv).
The TOC additionally reads
[../bvk2_aba_comprehensive_validation/benchmark_summary.csv](../bvk2_aba_comprehensive_validation/benchmark_summary.csv),
[../bvk2_aba_fixed_stiffness_validation/summary.csv](../bvk2_aba_fixed_stiffness_validation/summary.csv),
[../bvk2_publication_phase_diagram/phase_diagram_plot_data.csv](../bvk2_publication_phase_diagram/phase_diagram_plot_data.csv),
and
[../bvk2_publication_phase_diagram/scft_reference_boundaries.csv](../bvk2_publication_phase_diagram/scft_reference_boundaries.csv).

From the repository root, redraw the SVGs:

```bash
uv run python -c 'from scripts.render_macromolecules_figures import render_kernel_and_hierarchy; render_kernel_and_hierarchy()'
```

The PNG previews require the additional LibreOffice/Pillow path used by
`scripts/write_bvk2_figure1_and_toc.py`. Numerical regeneration of the source
tables is separate; see [REPRODUCING.md](../../REPRODUCING.md).
