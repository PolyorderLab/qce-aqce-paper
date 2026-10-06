# Software and dependencies

You need the Python plotting tools and Poppler to reproduce the SVG figures.
`reproduce_figures.sh` also uses LibreOffice to make PNG previews.
Julia is needed for numerical calculations, and Polyorder.jl is needed for new
SCFT calculations. See [REPRODUCING.md](REPRODUCING.md) for setup and commands.

## Tools for figure generation

- [uv](https://github.com/astral-sh/uv) installs and runs the locked Python
  environment. The project requires Python 3.10 or newer; figure reproduction
  was checked with Python 3.13.14 on Linux.
- [Poppler](https://poppler.freedesktop.org/) supplies `pdftocairo`, used to
  render the internally authored homopolymer source plots.
- [LibreOffice](https://www.libreoffice.org/) supplies the `libreoffice` command
  used for the Figure 1 and TOC PNG previews. It is required by
  `reproduce_figures.sh`, but can be omitted for the
  [SVG-only workflow](REPRODUCING.md#svg-only-workflow).

`uv sync --frozen` installs the packages listed in
[`pyproject.toml`](pyproject.toml) at the versions resolved in
[`uv.lock`](uv.lock). You do not need to install them individually:

| Package | Use in this repository |
|---|---|
| [Matplotlib](https://github.com/matplotlib/matplotlib) | Plotting and SVG export |
| [mpltex](https://github.com/liuyxpp/mpltex) | Journal-style plot formatting |
| [NumPy](https://github.com/numpy/numpy) | Numerical table processing |
| [Pillow](https://github.com/python-pillow/Pillow) | Image handling and previews |
| [pytest](https://github.com/pytest-dev/pytest) | Python tests |

## Paper-specific source

The Julia package in [`src/`](src/) is named `DFMMonteCarlo`. Useful entry
points for inspecting the implementation are:

| File | Contents |
|---|---|
| [`src/DFMMonteCarlo.jl`](src/DFMMonteCarlo.jl) | Main module and shared numerical routines |
| [`src/bvk2.jl`](src/bvk2.jl) | QCE/AQCE functional implementation, including the adaptive-stiffness option |
| [`src/bvk2_oversample.jl`](src/bvk2_oversample.jl) | Oversampled spectral discretization |
| [`src/bvk2_morphology.jl`](src/bvk2_morphology.jl) | Morphology-related routines |
| [`src/opf.jl`](src/opf.jl) | OPF comparison model |
| [`src/burp_ti_newton_gmres.jl`](src/burp_ti_newton_gmres.jl) | BURP numerical solver |

[`scripts/`](scripts/) contains the calculation drivers, table assemblers,
validators, and renderers. Use the [model-name key](DATA_GUIDE.md#model-names)
to identify QCE, AQCE, and the comparison models in filenames and tables.

## Julia dependencies

[`Project.toml`](Project.toml) lists the main project's dependencies, and
[`Manifest.toml`](Manifest.toml) records their exact versions. The manifest was
generated with [Julia 1.12.6](https://github.com/JuliaLang/julia).
Public dependencies are obtained from their upstream repositories rather than
bundled here:

- [Arianna.jl](https://github.com/TheDisorderedOrganization/Arianna.jl)
- [CSV.jl](https://github.com/JuliaData/CSV.jl)
- [DataFrames.jl](https://github.com/JuliaData/DataFrames.jl)
- [FFTW.jl](https://github.com/JuliaMath/FFTW.jl)
- [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl)
- [FourierTools.jl](https://github.com/bionanoimaging/FourierTools.jl)
- [MAT.jl](https://github.com/JuliaIO/MAT.jl)
- [OhMyThreads.jl](https://github.com/JuliaFolds2/OhMyThreads.jl)
- [Optim.jl](https://github.com/JuliaNLSolvers/Optim.jl)

The separate SCFT environment also records the public packages
[Polymer.jl](https://github.com/liuyxpp/Polymer.jl) and
[PolymerArchitecture.jl](https://github.com/liuyxpp/PolymerArchitecture.jl).
The manifests include additional transitive dependencies and Julia standard
libraries. For Python, install dependencies with `uv sync --frozen`, rather
than from `reproducibility/requirements.txt`.

## Polyorder.jl availability

The SCFT scripts use **Polyorder.jl 0.31.1** and the Scattering.jl dependency,
which are private and not included here. You need access to these packages to
run the SCFT scripts.

[`reproducibility/polyorder/Project.toml`](reproducibility/polyorder/Project.toml)
and its [manifest](reproducibility/polyorder/Manifest.toml) list the required
package versions. To work with the paper's SCFT results without installing
these packages, use the reference tables and profiles under `results/`.
