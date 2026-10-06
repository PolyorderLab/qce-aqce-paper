# Reproducing the results

There are two ways to use this repository: redraw the figures from the supplied
results, or run new numerical calculations. Figure reproduction is the simplest
starting point and requires only public software.

## Reproduce the figures

### 1. Prepare the tools

Install **uv**, **Poppler** (which provides `pdftocairo`), and **LibreOffice**
(for PNG previews). See
[SOFTWARE.md](SOFTWARE.md#tools-for-figure-generation) for upstream links.
uv installs the Python dependencies recorded in `uv.lock`.

Open a terminal in this repository's root directory, where `reproduce_figures.sh`
is located. Check that the tools are available:

```sh
uv --version
pdftocairo -v
libreoffice --version
```

The full workflow has been checked on Linux. It does not require Julia,
Polyorder.jl, or a LaTeX installation.

### 2. Verify the downloaded files

Before changing or regenerating files, run:

```sh
sha256sum --check checksums.sha256
```

Each listed file should report `OK`. A mismatch means that the file differs
from the supplied version. Check whether it was edited or regenerated before
interpreting the mismatch as a download problem.

### 3. Generate the figures

```sh
./reproduce_figures.sh
```

This script creates a local `.venv`, installs the locked Python dependencies,
and runs three plotting commands:

```sh
uv run python scripts/render_macromolecules_figures.py
uv run python scripts/render_macromolecules_figure8.py
uv run python scripts/write_bvk2_figure1_and_toc.py
```

Together they produce the eight manuscript figures, two Supporting Information
figures, and the table-of-contents graphic. `render_macromolecules_figure8.py`
renders both Figures 7 and 8. The renderers also
generate auxiliary diagnostic plots and previews. They read the supplied
numerical data without running field minimizations or SCFT calculations.

Open the resulting SVGs listed in the [figure index](DATA_GUIDE.md#figures-and-source-data).
They are written under `results/`, replacing the corresponding figure files.
The TOC renderer also writes a PDF under
`docs/manuscript/macromolecules/latex_review/figures/`.

Use a separate copy if you want to change the plots while preserving the
supplied figures.

## SVG-only workflow

If you do not need PNG previews, LibreOffice is unnecessary. Install uv and
Poppler, then run only:

```sh
uv sync --frozen
uv run python scripts/render_macromolecules_figures.py
uv run python scripts/render_macromolecules_figure8.py
```

These commands produce all manuscript, SI, and TOC SVGs. They omit the preview
conversion performed by `write_bvk2_figure1_and_toc.py`.

## Render selected figures

After `uv sync --frozen`, you can run individual renderers from the repository
root. For Figure 6:

```sh
uv run python -c 'from scripts.render_macromolecules_figures import render_phase_diagram; render_phase_diagram()'
```

For Figures 7 and 8:

```sh
uv run python scripts/render_macromolecules_figure8.py
```

For Figure 1 and the TOC, including PNG previews:

```sh
uv run python scripts/write_bvk2_figure1_and_toc.py
```

These commands read the supplied tables. To extract new literature curves or
recalculate phase boundaries, you will also need the reference artwork or
underlying calculation files.

## Numerical calculations

Run numerical experiments in a separate working copy because the calculation
scripts replace output tables. Choose an entry point below for the calculation
you want to repeat.

### Calculations using public dependencies

Use Julia **1.12.6** to match the supplied environment. Install the main
project's dependencies, then regenerate the quadratic-response tables:

```sh
julia --startup-file=no --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. scripts/write_diblock_kernel_comparison.jl
julia --project=. scripts/write_diblock_weak_response_map.jl
```

These scripts write the source tables for Figure 1. The model implementations
are described in [SOFTWARE.md](SOFTWARE.md#paper-specific-source).

### Calculations involving SCFT

The SCFT scripts require **Polyorder.jl 0.31.1** and its dependencies. If you
have access to these private packages, install the SCFT environment and add it
to Julia's package search path:

```sh
julia --startup-file=no --project=reproducibility/polyorder -e 'using Pkg; Pkg.instantiate()'
export JULIA_LOAD_PATH="@:${PWD}/reproducibility/polyorder:@stdlib"
```

Representative calculation entry points are:

```sh
julia --project=. scripts/write_liu2019_stress_free_period_map.jl
julia --project=. scripts/write_lamellar_cell_stress_mechanism.jl
julia --project=. scripts/write_bvk2_aba_comprehensive_validation.jl --recompute=true
```

Inspect each script's options and input paths before launching it. Defaults can
reuse existing results or cover a different subset from the final figure.
Repeating the full simulation campaigns requires additional raw case files and
intermediate solver states. Some tests also require those files; use the figure
workflow above to check reproduction with the files supplied here.

## If a command fails

- **`uv`, `pdftocairo`, or `libreoffice` is not found:** install the corresponding tool and make
  sure it is on your shell's `PATH`.
- **A plotting input is missing:** confirm that you have the complete repository
  and run the checksum check. Copying a plotting script alone is insufficient.
- **A renderer rejects a table or its hash:** the figure scripts check expected
  data and model inventories. Use the supplied tables for reproduction; when
  adapting the code, review the checks together with your intended data changes.
- **A private Julia package cannot be installed:** use the figure-only workflow
  to work with the supplied SCFT results. Installing the main project's public
  dependencies does not grant access to Polyorder.jl.
