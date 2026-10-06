# JCP 2021 AB phase-diagram digitization

This directory contains numerical curves extracted from the top AB-diblock panel of
Li and Liu, *J. Chem. Phys.* **154**, 014903 (2021),
[DOI 10.1063/5.0037979](https://doi.org/10.1063/5.0037979). The source contains
nine coexistence curves and preserves the close-packed-sphere and
\(O^{70}\) (Fddd) pockets, triple points, and critical point.

## Data files

- [jcp2021_ab_phase_boundaries.csv](jcp2021_ab_phase_boundaries.csv) contains
  the extracted source vertices, separated by transition and symmetry side;
  no curve fitting or smoothing was applied.
- [jcp2021_ab_topology_nodes.csv](jcp2021_ab_topology_nodes.csv) contains the
  critical and triple-point coordinates and their incident curves.
- [jcp2021_ab_extraction_audit.csv](jcp2021_ab_extraction_audit.csv) retains the
  numerical overlay and mirror-symmetry residuals.
- [extract_from_eps.py](extract_from_eps.py) is the extractor and validator.

Use the CSV tables directly to plot the reference phase boundaries. To repeat
the extraction, you need the EPS artwork from the cited article, which is not
included here.

With an authorized EPS copy, run from the repository root:

```bash
uv run python results/jcp2021_ab_phase_diagram_digitization/extract_from_eps.py \
  --source /path/to/AB_phase_diagram_fig1_jcp2021.eps
```

The vector-to-data mapping is \(f_A=(x-954)/1450\) and
\(\chi N=50(1270-y)/1172\). The extracted critical point is
`(fA, chiN)=(0.500000, 10.537543)`. See
[REPRODUCING.md](../../REPRODUCING.md) for the broader figure workflow.
