# QCE/AQCE: data, code, and figures

Companion materials for **Improved Density Functionals for Predicting Block
Copolymer Domain Structure and Phase Behavior**, by Yi-Xin Liu.

The paper introduces the quadratic-connectivity entropic (QCE) density
functional and its adaptive-stiffness extension, AQCE. This repository provides
the numerical results, model implementations, and plotting scripts behind the
AB-diblock and symmetric ABA-triblock comparisons.

You can inspect the data and reproduce the figures without Polyorder.jl or
running any new simulations.

## Start here

| I want to… | Where to start |
|---|---|
| Find a figure and its source tables | [Data guide](DATA_GUIDE.md#figures-and-source-data) |
| Understand model names, columns, and units | [Reading the tables](DATA_GUIDE.md#reading-the-tables) |
| Regenerate the figures | [Reproduction instructions](REPRODUCING.md#reproduce-the-figures) |
| Inspect the model implementation | [Paper-specific source](SOFTWARE.md#paper-specific-source) |
| Rerun numerical calculations | [Numerical calculations](REPRODUCING.md#numerical-calculations) |

The figures are supplied as SVG files, which you can open in a browser. After
installing the [figure-generation tools](SOFTWARE.md#tools-for-figure-generation),
run these commands from the repository root to verify the files and redraw them:

```sh
sha256sum --check checksums.sha256
./reproduce_figures.sh
```

## Repository contents

- [`results/`](results/): numerical tables, figures, and dataset-specific notes.
- [`src/`](src/): Julia implementations of the density functionals and numerical methods.
- [`scripts/`](scripts/): calculation, analysis, validation, and plotting scripts.
- [`test/`](test/): numerical and plotting tests. Some require additional
  campaign files; see the [reproduction guide](REPRODUCING.md).
- [`accepted_roots/`](accepted_roots/): phase-boundary coordinates and associated
  resolution, provenance, and crossing evidence.
- [`SOFTWARE.md`](SOFTWARE.md): dependencies, version records, and software availability.

In filenames and tables, `bvk2_fixed` identifies QCE and `bvk2` identifies AQCE.
See the [model-name key](DATA_GUIDE.md#model-names) for the comparison models.

## Software access and reuse

Paper-specific code is included. Public software dependencies are linked in
[SOFTWARE.md](SOFTWARE.md), rather than copied into this repository. SCFT
calculations used the private Polyorder.jl package, whose source is not included.
The supplied SCFT results are sufficient for figure reproduction. See
[REPRODUCING.md](REPRODUCING.md) for the requirements of new calculations.

For questions about the materials or permission to reuse them, contact
Yi-Xin Liu at [lyx@fudan.edu.cn](mailto:lyx@fudan.edu.cn).
