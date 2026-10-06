#!/bin/sh
set -eu

uv sync --frozen
uv run python scripts/render_macromolecules_figures.py
uv run python scripts/render_macromolecules_figure8.py
uv run python scripts/write_bvk2_figure1_and_toc.py
