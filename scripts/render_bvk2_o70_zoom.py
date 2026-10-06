#!/usr/bin/env python3
"""Render the low-chi BVK2 O70 pocket with a smooth JCP SCFT overlay.

The stored digitized SCFT vertices remain unchanged. Only the display path is
dequantized and densified under endpoint, monotonicity, and displacement gates.
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams.update(
    {
        "svg.fonttype": "none",
        "svg.hashsalt": "bvk2-o70-pocket-v1",
    }
)
import matplotlib.pyplot as plt
import mpltex
import mpltex.acs as mpltex_acs
from matplotlib.lines import Line2D

from assemble_bvk2_publication_phase_diagram import load_rpa_stability_limit
from scft_display import smooth_scft_display_curve


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
DEFAULT_BOUNDARIES = (
    RESULTS / "bvk2_publication_phase_diagram" / "phase_diagram_plot_data.csv"
)
DEFAULT_BRACKETS = (
    RESULTS / "bvk2_fcc_o70_boundary_campaign" / "boundary_brackets.csv"
)
DEFAULT_SCFT = (
    RESULTS / "bvk2_publication_phase_diagram" / "scft_reference_boundaries.csv"
)
DEFAULT_RPA_STABILITY = (
    RESULTS / "bvk2_publication_phase_diagram" / "rpa_stability_limit.csv"
)
DEFAULT_OUTPUT_PREFIX = (
    RESULTS / "bvk2_publication_phase_diagram" / "o70_region_zoom"
)
DEFAULT_TRIPLE_POINT = (
    RESULTS / "bvk2_publication_phase_diagram"
    / "estimated_g_c_o70_triple_point.csv"
)
DEFAULT_UPPER_TRIPLE_POINT = (
    RESULTS / "bvk2_publication_phase_diagram"
    / "estimated_g_o70_l_triple_point.csv"
)

SVG_METADATA = {"Date": None, "Creator": "Matplotlib + mpltex ACS via uv"}
TRANSITIONS = ("S/DIS", "S/C", "G/C", "G/L", "C/O70", "G/O70", "O70/L")
COLORS = {
    "S/C": "#d62728",
    "S/DIS": "#1f77b4",
    "G/C": "#9467bd",
    "G/L": "#ff7f0e",
    "C/O70": "#17becf",
    "G/O70": "#bcbd22",
    "O70/L": "#393b79",
}
MARKERS = {
    "S/C": "o",
    "S/DIS": "s",
    "G/C": "D",
    "G/L": "+",
    "C/O70": "^",
    "G/O70": "d",
    "O70/L": "h",
}
SCFT_COLOR = mpltex.almost_black


def read_rows(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def number(value: str | None) -> float | None:
    try:
        result = float(value or "")
    except ValueError:
        return None
    return result if result == result else None


def normalize_transition(value: str) -> str | None:
    if value in {"DIS/S", "DIS/S_cp", "S_cp/S"}:
        return "S/DIS"
    if value == "C/G":
        return "G/C"
    return value if value in TRANSITIONS else None


def accepted_bvk2(path: Path) -> dict[str, list[tuple[float, float]]]:
    grouped: dict[str, list[tuple[float, float]]] = defaultdict(list)
    for row in read_rows(path):
        transition = normalize_transition(row.get("transition", ""))
        f_a = number(row.get("fA"))
        chi_n = number(row.get("chiN"))
        if transition is None or f_a is None or chi_n is None:
            continue
        if row.get("status") not in {"accepted", "critical_endpoint"}:
            continue
        grouped[transition].append((f_a, chi_n))
    return {
        transition: sorted(set(points), key=lambda point: point[1])
        for transition, points in grouped.items()
    }


def provisional_o70(path: Path) -> list[tuple[str, float, float]]:
    points: list[tuple[str, float, float]] = []
    for row in read_rows(path):
        transition = normalize_transition(row.get("boundary", ""))
        f_a = number(row.get("fA_estimate"))
        chi_n = number(row.get("chiN"))
        if (
            row.get("status") == "provisional"
            and transition in {"C/O70", "G/O70", "O70/L"}
            and f_a is not None
            and chi_n is not None
        ):
            points.append((transition, f_a, chi_n))
    return points


def estimated_triple_point(path: Path) -> tuple[float, float]:
    rows = read_rows(path)
    if len(rows) != 1:
        raise RuntimeError(f"expected one estimated triple point in {path}")
    row = rows[0]
    if (
        row.get("claim_kind") != "estimated_triple_point"
        or row.get("status") != "estimated_display_only"
    ):
        raise RuntimeError(f"invalid estimated triple-point contract in {path}")
    f_a = number(row.get("fA"))
    chi_n = number(row.get("chiN"))
    if f_a is None or chi_n is None:
        raise RuntimeError(f"non-finite estimated triple point in {path}")
    return f_a, chi_n


def close_o70_pocket(
    curves: dict[str, list[tuple[float, float]]],
    point: tuple[float, float],
) -> None:
    for transition in ("G/C", "G/O70", "C/O70"):
        curves.setdefault(transition, []).append(point)
        curves[transition] = sorted(set(curves[transition]), key=lambda item: item[1])


def close_upper_o70_pocket(
    curves: dict[str, list[tuple[float, float]]],
    point: tuple[float, float],
) -> None:
    for transition in ("G/O70", "O70/L", "G/L"):
        curves.setdefault(transition, []).append(point)
        curves[transition] = sorted(set(curves[transition]), key=lambda item: item[1])


def curves_with_estimated_triple_points(
    accepted: dict[str, list[tuple[float, float]]],
    lower_triple: tuple[float, float],
    upper_triple: tuple[float, float],
) -> dict[str, list[tuple[float, float]]]:
    """Return curve geometry without relabeling junctions as accepted roots."""
    curves = {
        transition: list(points) for transition, points in accepted.items()
    }
    close_o70_pocket(curves, lower_triple)
    close_upper_o70_pocket(curves, upper_triple)
    return curves


def scft_curves(path: Path) -> list[tuple[str, list[tuple[float, float]]]]:
    grouped: dict[tuple[str, str, str], list[tuple[int, float, float]]] = defaultdict(list)
    for row in read_rows(path):
        transition = normalize_transition(row.get("transition", ""))
        if transition is None or row.get("side") != "left":
            continue
        f_a = number(row.get("fA"))
        chi_n = number(row.get("chiN"))
        if f_a is None or chi_n is None:
            continue
        grouped[(transition, row.get("side", ""), row.get("curve_id", ""))].append(
            (int(row.get("point_index", "0")), f_a, chi_n)
        )
    curves: list[tuple[str, list[tuple[float, float]]]] = []
    for (transition, _side, _curve), points in grouped.items():
        ordered = [(f_a, chi_n) for _, f_a, chi_n in sorted(points)]
        curves.append((transition, ordered))
    return curves




@mpltex.acs_decorator
def render(
    boundaries: Path,
    brackets: Path,
    scft: Path,
    triple_point: Path,
    upper_triple_point: Path,
    output_prefix: Path,
    rpa_stability: Path = DEFAULT_RPA_STABILITY,
    *,
    show_provisional: bool = False,
) -> None:
    accepted = accepted_bvk2(boundaries)
    lower_triple = estimated_triple_point(triple_point)
    upper_triple = estimated_triple_point(upper_triple_point)
    curves = curves_with_estimated_triple_points(
        accepted, lower_triple, upper_triple
    )
    # The publication zoom mirrors Figure 6 and therefore shows accepted
    # boundary claims only. Provisional campaign locators are available as an
    # explicit diagnostic overlay, never as part of the default figure.
    provisional = provisional_o70(brackets) if show_provisional else []
    fig, ax = plt.subplots(
        figsize=(mpltex_acs.width_double_column, mpltex_acs.width_double_column),
        constrained_layout=True,
    )

    rpa_curve = [
        point for point in load_rpa_stability_limit(rpa_stability)
        if 0.4 <= point[0] <= 0.5 and 10.0 <= point[1] <= 16.0
    ]
    ax.plot(
        [point[0] for point in rpa_curve],
        [point[1] for point in rpa_curve],
        color=SCFT_COLOR,
        linewidth=1.0,
        linestyle=":",
        alpha=0.75,
        zorder=2,
    )

    for _transition, raw_curve in scft_curves(scft):
        curve = smooth_scft_display_curve(raw_curve)
        ax.plot(
            [point[0] for point in curve],
            [point[1] for point in curve],
            color=SCFT_COLOR,
            linewidth=0.8,
            linestyle=(0, (4, 3)),
            alpha=0.75,
            zorder=1,
        )

    visible: list[str] = []
    for transition in TRANSITIONS:
        curve_points = curves.get(transition, [])
        accepted_points = accepted.get(transition, [])
        if not curve_points:
            continue
        if any(
            0.39 <= f_a <= 0.51 and 9.8 <= chi_n <= 16.2
            for f_a, chi_n in curve_points
        ):
            visible.append(transition)
        ax.plot(
            [point[0] for point in curve_points],
            [point[1] for point in curve_points],
            color=COLORS[transition],
            linewidth=1.2,
            zorder=3,
        )
        ax.plot(
            [point[0] for point in accepted_points],
            [point[1] for point in accepted_points],
            linestyle="none",
            color=COLORS[transition],
            marker=MARKERS[transition],
            markerfacecolor="white",
            markeredgewidth=0.9,
            markersize=4.2,
            zorder=4,
        )

    provisional_transitions: set[str] = set()
    for transition, f_a, chi_n in provisional:
        if 0.4 <= f_a <= 0.5 and 10.0 <= chi_n <= 16.0:
            provisional_transitions.add(transition)
            ax.plot(
                f_a,
                chi_n,
                linestyle="none",
                marker=MARKERS[transition],
                markerfacecolor="white",
                markeredgecolor=COLORS[transition],
                markeredgewidth=1.0,
                markersize=5.0,
                zorder=4,
            )

    for triple in (lower_triple, upper_triple):
        ax.plot(
            triple[0],
            triple[1],
            linestyle="none",
            marker="*",
            markerfacecolor=SCFT_COLOR,
            markeredgecolor=SCFT_COLOR,
            markeredgewidth=1.0,
            markersize=5.0,
            zorder=6,
        )

    # Retain the published 0.40 endpoint tick while leaving enough interior
    # padding to show an accepted point exactly at f_A=0.40 in full.
    ax.set_xlim(0.3995, 0.5)
    ax.set_ylim(10.0, 16.0)
    ax.set_xticks([0.40, 0.42, 0.44, 0.46, 0.48, 0.50])
    ax.set_yticks([10, 11, 12, 13, 14, 15, 16])
    ax.set_xlabel(r"$f_A$")
    ax.set_ylabel(r"$\chi N$")
    ax.text(
        0.441,
        11.72,
        "O70",
        color=COLORS["O70/L"],
        fontsize=8,
        ha="center",
        va="center",
        zorder=5,
    )
    ax.grid(False)
    for spine in ax.spines.values():
        spine.set_visible(True)

    handles = [
        Line2D(
            [0], [0],
            color=COLORS[transition],
            marker=MARKERS[transition],
            markerfacecolor="white",
            markeredgewidth=0.9,
            markersize=4.2,
            linewidth=1.2,
            label=transition,
        )
        for transition in visible
    ]
    if provisional_transitions:
        handles.append(
            Line2D(
                [0], [0],
                color=mpltex.almost_black,
                marker="o",
                markerfacecolor="white",
                linestyle="none",
                markersize=4.5,
                label="provisional BVK2",
            )
        )
    handles.append(
        Line2D(
            [0], [0],
            color=SCFT_COLOR,
            linewidth=1.0,
            linestyle=":",
            label="RPA stability limit",
        )
    )
    handles.append(
        Line2D(
            [0], [0],
            color=SCFT_COLOR,
            marker="*",
            markerfacecolor=SCFT_COLOR,
            linestyle="none",
            markersize=4.5,
            label="estimated triple point",
        )
    )
    handles.append(
        Line2D(
            [0], [0],
            color=SCFT_COLOR,
            linewidth=0.8,
            linestyle=(0, (4, 3)),
            label="SCFT",
        )
    )
    ax.legend(
        handles=handles,
        loc="upper right",
        frameon=False,
        borderaxespad=0.55,
        labelspacing=0.35,
        handlelength=2.4,
    )

    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    svg_path = output_prefix.with_suffix(".svg")
    fig.savefig(
        svg_path,
        format="svg",
        metadata=SVG_METADATA,
    )
    svg_path.write_text(
        "\n".join(
            line.rstrip()
            for line in svg_path.read_text(encoding="utf-8").splitlines()
        )
        + "\n",
        encoding="utf-8",
    )
    fig.savefig(
        output_prefix.with_suffix(".png"),
        format="png",
        dpi=240,
        metadata={"Software": "Matplotlib + mpltex ACS via uv"},
    )
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--boundaries", type=Path, default=DEFAULT_BOUNDARIES)
    parser.add_argument("--brackets", type=Path, default=DEFAULT_BRACKETS)
    parser.add_argument("--scft", type=Path, default=DEFAULT_SCFT)
    parser.add_argument(
        "--rpa-stability", type=Path, default=DEFAULT_RPA_STABILITY
    )
    parser.add_argument("--triple-point", type=Path, default=DEFAULT_TRIPLE_POINT)
    parser.add_argument(
        "--upper-triple-point", type=Path, default=DEFAULT_UPPER_TRIPLE_POINT
    )
    parser.add_argument("--output-prefix", type=Path, default=DEFAULT_OUTPUT_PREFIX)
    parser.add_argument(
        "--show-provisional",
        action="store_true",
        help="overlay provisional O70 campaign locators for diagnostics",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    render(
        args.boundaries.resolve(),
        args.brackets.resolve(),
        args.scft.resolve(),
        args.triple_point.resolve(),
        args.upper_triple_point.resolve(),
        args.output_prefix.resolve(),
        args.rpa_stability.resolve(),
        show_provisional=args.show_provisional,
    )
