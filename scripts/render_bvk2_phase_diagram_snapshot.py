#!/usr/bin/env python3
"""Render Figure 6 from the frozen canonical publication-plot snapshot."""

from __future__ import annotations

import argparse
from pathlib import Path

from render_macromolecules_figures import render_phase_diagram


PROJECT = Path(__file__).resolve().parents[1]
DEFAULT_PLOT_DATA = (
    PROJECT / "results/bvk2_publication_phase_diagram/phase_diagram_plot_data.csv"
)
DEFAULT_SCFT = (
    PROJECT / "results/bvk2_publication_phase_diagram/scft_reference_boundaries.csv"
)
DEFAULT_RPA = (
    PROJECT / "results/bvk2_publication_phase_diagram/rpa_stability_limit.csv"
)
DEFAULT_FCC_TRIPLE = (
    PROJECT / "results/bvk2_publication_phase_diagram/fcc_region_zoom_triple_point.csv"
)
DEFAULT_G_C_O70_TRIPLE = (
    PROJECT
    / "results/bvk2_publication_phase_diagram/estimated_g_c_o70_triple_point.csv"
)
DEFAULT_G_O70_L_TRIPLE = (
    PROJECT
    / "results/bvk2_publication_phase_diagram/estimated_g_o70_l_triple_point.csv"
)
DEFAULT_OUTPUT = (
    PROJECT / "results/bvk2_publication_phase_diagram/phase_diagram_snapshot.svg"
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plot-data", type=Path, default=DEFAULT_PLOT_DATA)
    parser.add_argument("--scft", type=Path, default=DEFAULT_SCFT)
    parser.add_argument("--rpa", type=Path, default=DEFAULT_RPA)
    parser.add_argument("--fcc-triple", type=Path, default=DEFAULT_FCC_TRIPLE)
    parser.add_argument(
        "--g-c-o70-triple", type=Path, default=DEFAULT_G_C_O70_TRIPLE
    )
    parser.add_argument(
        "--g-o70-l-triple", type=Path, default=DEFAULT_G_O70_L_TRIPLE
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    args.output.resolve().parent.mkdir(parents=True, exist_ok=True)
    render_phase_diagram(
        plot_data_path=args.plot_data.resolve(),
        scft_path=args.scft.resolve(),
        output_path=args.output.resolve(),
        rpa_path=args.rpa.resolve(),
        triple_point_sources=(
            (
                args.fcc_triple.resolve(),
                "fA_estimate",
                "chiN_estimate",
                ("S/DIS", "FCC/DIS", "FCC/BCC"),
            ),
            (
                args.g_c_o70_triple.resolve(),
                "fA",
                "chiN",
                ("G/C", "G/O70", "C/O70"),
            ),
            (
                args.g_o70_l_triple.resolve(),
                "fA",
                "chiN",
                ("G/L", "G/O70", "O70/L"),
            ),
        ),
    )
    with args.plot_data.resolve().open() as handle:
        accepted_count = sum(1 for _ in handle) - 1
    print(
        f"rendered {accepted_count} accepted roots to "
        f"{args.output.resolve()}"
    )


if __name__ == "__main__":
    main()
