#!/usr/bin/env python3
"""Compatibility entry point for Figure 3 and the Macromolecules TOC graphic.

The authoritative publication renderer is ``render_macromolecules_figures.py``.
This legacy entry point delegates to it so older reproduction commands cannot
overwrite Figure 3 with the retired four-panel layout.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import html
import shutil
import subprocess
import tempfile
from collections import Counter
from pathlib import Path

from PIL import Image


PROJECT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = PROJECT / "results/bvk2_model_hierarchy_graphics"

KERNEL = PROJECT / "results/diblock_kernel_comparison/kernel_samples.csv"
WEAK = PROJECT / "results/diblock_weak_response_map/weak_response_map.csv"
PROFILES = PROJECT / "results/unified_lamellar_benchmark/profiles.csv"
SUMMARY = PROJECT / "results/unified_lamellar_benchmark/summary.csv"
BOUNDARIES = (
    PROJECT / "results/bvk2_publication_phase_diagram/accepted_phase_boundaries.csv"
)
SCFT_BOUNDARIES = (
    PROJECT / "results/bvk2_publication_phase_diagram/scft_reference_boundaries.csv"
)

BLUE = "#1f5fbf"
GREEN = "#23864a"
PURPLE = "#7551a8"
RED = "#d94141"
ORANGE = "#e78520"
TEAL = "#1297a6"
INK = "#18212b"
MID = "#52606d"
GRID = "#d9e0e7"
PALE_BLUE = "#edf4fc"
PALE_GREEN = "#edf8f1"
PALE_ORANGE = "#fff4e7"


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def esc(value: object) -> str:
    return html.escape(str(value), quote=True)


def f3(value: float) -> str:
    return f"{value:.3f}"


class SVG:
    def __init__(
        self,
        width: int,
        height: int,
        *,
        physical_width: str | None = None,
        physical_height: str | None = None,
        title: str,
        description: str,
    ) -> None:
        width_attr = physical_width or str(width)
        height_attr = physical_height or str(height)
        self.parts = [
            (
                f'<svg xmlns="http://www.w3.org/2000/svg" '
                f'width="{width_attr}" height="{height_attr}" '
                f'viewBox="0 0 {width} {height}" role="img" '
                f'aria-labelledby="svg-title svg-desc">'
            ),
            f"<title id=\"svg-title\">{esc(title)}</title>",
            f"<desc id=\"svg-desc\">{esc(description)}</desc>",
            "<defs>",
            (
                '<marker id="arrow" markerWidth="10" markerHeight="8" '
                'refX="9" refY="4" orient="auto">'
                f'<path d="M0,0 L10,4 L0,8 Z" fill="{MID}"/>'
                "</marker>"
            ),
            (
                '<marker id="arrow-blue" markerWidth="10" markerHeight="8" '
                'refX="9" refY="4" orient="auto">'
                f'<path d="M0,0 L10,4 L0,8 Z" fill="{BLUE}"/>'
                "</marker>"
            ),
            (
                "<style>"
                "text{font-family:Arial,Helvetica,sans-serif}"
                ".axis{stroke:#18212b;stroke-width:1.5}"
                ".grid{stroke:#d9e0e7;stroke-width:1}"
                ".curve{fill:none;stroke-width:3;stroke-linecap:round;"
                "stroke-linejoin:round}"
                ".thin{fill:none;stroke-width:2.1;stroke-linecap:round;"
                "stroke-linejoin:round}"
                ".box{stroke:#b9c5d0;stroke-width:1.5}"
                "</style>"
            ),
            "</defs>",
            '<rect width="100%" height="100%" fill="#ffffff"/>',
        ]

    def add(self, markup: str) -> None:
        self.parts.append(markup)

    def text(
        self,
        x: float,
        y: float,
        value: object,
        *,
        size: float = 22,
        weight: int = 400,
        anchor: str = "start",
        fill: str = INK,
        extra: str = "",
    ) -> None:
        self.add(
            f'<text x="{f3(x)}" y="{f3(y)}" font-size="{size}" '
            f'font-weight="{weight}" text-anchor="{anchor}" '
            f'fill="{fill}" {extra}>{esc(value)}</text>'
        )

    def line(
        self,
        x1: float,
        y1: float,
        x2: float,
        y2: float,
        *,
        stroke: str = INK,
        width: float = 2,
        dash: str = "",
        marker: str = "",
    ) -> None:
        dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
        marker_attr = f' marker-end="url(#{marker})"' if marker else ""
        self.add(
            f'<line x1="{f3(x1)}" y1="{f3(y1)}" x2="{f3(x2)}" '
            f'y2="{f3(y2)}" stroke="{stroke}" stroke-width="{width}"'
            f"{dash_attr}{marker_attr}/>"
        )

    def rect(
        self,
        x: float,
        y: float,
        width: float,
        height: float,
        *,
        fill: str = "white",
        stroke: str = "#b9c5d0",
        radius: float = 12,
        stroke_width: float = 1.5,
    ) -> None:
        self.add(
            f'<rect x="{f3(x)}" y="{f3(y)}" width="{f3(width)}" '
            f'height="{f3(height)}" rx="{f3(radius)}" fill="{fill}" '
            f'stroke="{stroke}" stroke-width="{stroke_width}"/>'
        )

    def circle(
        self,
        x: float,
        y: float,
        radius: float,
        *,
        fill: str,
        stroke: str = "white",
        stroke_width: float = 1.5,
    ) -> None:
        self.add(
            f'<circle cx="{f3(x)}" cy="{f3(y)}" r="{f3(radius)}" '
            f'fill="{fill}" stroke="{stroke}" stroke-width="{stroke_width}"/>'
        )

    def polyline(
        self,
        points: list[tuple[float, float]],
        *,
        stroke: str,
        width: float = 3,
        dash: str = "",
        css_class: str = "curve",
        fill: str = "none",
    ) -> None:
        points_attr = " ".join(f"{f3(x)},{f3(y)}" for x, y in points)
        dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
        self.add(
            f'<polyline class="{css_class}" points="{points_attr}" '
            f'stroke="{stroke}" stroke-width="{width}" fill="{fill}"'
            f"{dash_attr}/>"
        )

    def path(
        self,
        d: str,
        *,
        stroke: str,
        width: float = 2,
        fill: str = "none",
        dash: str = "",
    ) -> None:
        dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
        self.add(
            f'<path d="{d}" stroke="{stroke}" stroke-width="{width}" '
            f'fill="{fill}" stroke-linecap="round" '
            f'stroke-linejoin="round"{dash_attr}/>'
        )

    def finish(self, path: Path) -> None:
        self.parts.append("</svg>")
        path.write_text("\n".join(self.parts) + "\n", encoding="utf-8")


def map_x(value: float, x0: float, width: float, lo: float, hi: float) -> float:
    return x0 + width * (value - lo) / (hi - lo)


def map_y(value: float, y0: float, height: float, lo: float, hi: float) -> float:
    return y0 + height * (hi - value) / (hi - lo)


def panel_label(svg: SVG, x: float, y: float, label: str) -> None:
    svg.circle(x, y, 19, fill=INK, stroke=INK)
    svg.text(x, y + 7, label, size=22, weight=700, anchor="middle", fill="white")


def draw_axes(
    svg: SVG,
    *,
    x0: float,
    y0: float,
    width: float,
    height: float,
    xlim: tuple[float, float],
    ylim: tuple[float, float],
    xticks: list[float],
    yticks: list[float],
    xlabel: str,
    ylabel: str,
    tick_size: int = 17,
) -> None:
    for tick in xticks:
        x = map_x(tick, x0, width, *xlim)
        svg.line(x, y0, x, y0 + height, stroke=GRID, width=1)
        svg.text(x, y0 + height + 26, f"{tick:g}", size=tick_size, anchor="middle")
    for tick in yticks:
        y = map_y(tick, y0, height, *ylim)
        svg.line(x0, y, x0 + width, y, stroke=GRID, width=1)
        svg.text(x0 - 12, y + 6, f"{tick:g}", size=tick_size, anchor="end")
    svg.line(x0, y0 + height, x0 + width, y0 + height, width=1.6)
    svg.line(x0, y0, x0, y0 + height, width=1.6)
    svg.text(x0 + width / 2, y0 + height + 57, xlabel, size=21, anchor="middle")
    svg.add(
        f'<text transform="translate({f3(x0 - 60)} '
        f'{f3(y0 + height / 2)}) rotate(-90)" font-size="21" '
        f'text-anchor="middle">{esc(ylabel)}</text>'
    )


def select_kernel(
    kernel_rows: list[dict[str, str]], f_value: float
) -> dict[str, list[tuple[float, float]]]:
    selected: dict[str, list[tuple[float, float]]] = {}
    for model in ("exact_rpa", "uneyama_doi", "bvk1_bvk2", "ok"):
        points = [
            (float(row["qRg"]), float(row["susceptibility_over_exact_smax"]))
            for row in kernel_rows
            if row["model"] == model
            and abs(float(row["f"]) - f_value) < 1e-12
            and 0.5 <= float(row["qRg"]) <= 3.6
        ]
        selected[model] = points[::3]
    return selected


def weak_curves(
    weak_rows: list[dict[str, str]], field: str
) -> dict[str, list[tuple[float, float]]]:
    curves: dict[str, list[tuple[float, float]]] = {}
    for model in ("exact_rpa", "uneyama_doi", "bvk1_bvk2", "ok"):
        curves[model] = sorted(
            (
                float(row["f"]),
                100.0 * float(row[field]),
            )
            for row in weak_rows
            if row["model"] == model
        )
    return curves


def draw_legend(
    svg: SVG,
    x: float,
    y: float,
    entries: list[tuple[str, str, str]],
    *,
    spacing: float = 31,
    font_size: int = 16,
) -> None:
    for index, (label, color, dash) in enumerate(entries):
        yy = y + index * spacing
        svg.line(x, yy, x + 38, yy, stroke=color, width=3, dash=dash)
        svg.text(x + 49, yy + 6, label, size=font_size)


def write_figure1(
    output: Path,
    kernel_rows: list[dict[str, str]],
    weak_rows: list[dict[str, str]],
) -> None:
    svg = SVG(
        1800,
        1050,
        title="QCE/AQCE model hierarchy and quantitative weak response",
        description=(
            "Four-panel publication figure: model lineage, normalized Gaussian "
            "susceptibility, spinodal error, and weak-period error."
        ),
    )
    svg.add(
        "<metadata>"
        "schema=bvk2-model-hierarchy-figure-v1;"
        "kernel_fA=0.35;weak_response_fA=0.10:0.50;"
        "AQCE_second_order_vertex=QCE"
        "</metadata>"
    )

    panel_label(svg, 38, 42, "a")
    svg.text(74, 43, "Model hierarchy", size=25, weight=700)
    svg.text(74, 70, "chain statistics → local mechanics", size=18, weight=600, fill=MID)

    x, w = 55.0, 520.0
    boxes = [
        (82.0, 122.0, PALE_BLUE, BLUE, "SCFT", "chain propagators; mean-field reference"),
        (82.0, 278.0, PALE_BLUE, BLUE, "Exact diblock RPA", "ideal-chain composition vertex"),
        (82.0, 466.0, "#f7f9fb", "#8b99a8", "Reduced UD Gaussian skeleton", "A/k² + C + Bk²"),
        (82.0, 647.0, PALE_GREEN, GREEN, "QCE / AQCE", "quadratic connectivity; adaptive gradient refinement"),
    ]
    for bx, by, fill, stroke, title, subtitle in boxes:
        svg.rect(bx, by, w - 54, 104, fill=fill, stroke=stroke, radius=14, stroke_width=2)
        svg.text(bx + 22, by + 39, title, size=23, weight=700, fill=stroke)
        svg.text(bx + 22, by + 73, subtitle, size=17, fill=MID)

    svg.line(315, 226, 315, 273, stroke=BLUE, width=2.5, marker="arrow-blue")
    svg.text(330, 253, "linearize", size=15, fill=BLUE)

    svg.rect(82, 389, 213, 68, fill=PALE_BLUE, stroke=BLUE, radius=11, stroke_width=2)
    svg.text(188.5, 418, "BURP-TI", size=21, weight=700, anchor="middle", fill=BLUE)
    svg.text(188.5, 443, "exact Gaussian order", size=14, anchor="middle", fill=MID)
    svg.path("M252 382 C245 365,225 358,202 358", stroke=BLUE, width=2.3)
    svg.line(202, 358, 202, 384, stroke=BLUE, width=2.3, marker="arrow-blue")
    svg.text(109, 364, "retain full vertex", size=14, fill=BLUE)

    svg.path("M373 382 C389 408,390 431,370 463", stroke=MID, width=2.2)
    svg.line(370, 449, 370, 462, stroke=MID, width=2.2, marker="arrow")
    svg.text(392, 414, "local reduction", size=14, fill=MID)

    svg.line(315, 570, 315, 642, stroke=GREEN, width=2.5, marker="arrow")
    svg.text(330, 603, "same Gaussian vertex", size=15, fill=GREEN)
    svg.text(330, 626, "adaptive term starts at O(δφ⁴)", size=15, fill=GREEN)

    svg.rect(82, 797, w - 54, 160, fill=PALE_ORANGE, stroke=ORANGE, radius=14, stroke_width=2)
    svg.text(103, 832, "The discriminating observable is nonlinear", size=21, weight=700, fill=ORANGE)
    svg.text(103, 869, "weak response", size=18, weight=700)
    svg.line(225, 863, 327, 863, stroke=ORANGE, width=2.3, marker="arrow")
    svg.text(344, 869, "cell stress", size=18, weight=700)
    svg.text(103, 907, "vertex / spinodal", size=16, fill=MID)
    svg.line(245, 901, 327, 901, stroke=ORANGE, width=2.3, marker="arrow")
    svg.text(344, 907, "stress-free period", size=16, fill=MID)
    svg.text(103, 938, "Matching a weak minimum does not fix finite-amplitude mechanics.", size=15, fill=MID)

    # Panel B: normalized susceptibility at fA = 0.35.
    panel_label(svg, 642, 42, "b")
    svg.text(680, 49, "Gaussian susceptibility at fA = 0.35", size=24, weight=700)
    x0, y0, pw, ph = 705.0, 94.0, 1000.0, 250.0
    xlim, ylim = (0.5, 3.6), (0.0, 1.05)
    draw_axes(
        svg,
        x0=x0,
        y0=y0,
        width=pw,
        height=ph,
        xlim=xlim,
        ylim=ylim,
        xticks=[0.5, 1, 1.5, 2, 2.5, 3, 3.5],
        yticks=[0, 0.25, 0.5, 0.75, 1],
        xlabel="q Rg",
        ylabel="S(q) / max S_RPA",
    )
    kernel = select_kernel(kernel_rows, 0.35)
    model_style = {
        "exact_rpa": (BLUE, ""),
        "uneyama_doi": (ORANGE, "3 5"),
        "bvk1_bvk2": (GREEN, "10 6"),
        "ok": (PURPLE, "3 6"),
    }
    for model in ("exact_rpa", "uneyama_doi", "bvk1_bvk2", "ok"):
        points = [
            (map_x(q, x0, pw, *xlim), map_y(value, y0, ph, *ylim))
            for q, value in kernel[model]
        ]
        svg.polyline(points, stroke=model_style[model][0], dash=model_style[model][1])
        q_peak, y_peak = max(kernel[model], key=lambda point: point[1])
        svg.circle(
            map_x(q_peak, x0, pw, *xlim),
            map_y(y_peak, y0, ph, *ylim),
            5.2,
            fill=model_style[model][0],
        )
    draw_legend(
        svg,
        1280,
        122,
        [
            ("Full RPA (BURP-TI)", BLUE, ""),
            ("Literal nonlinear UD Hessian", ORANGE, "3 5"),
            ("Reduced QCE/AQCE vertex", GREEN, "10 6"),
            ("OK", PURPLE, "3 6"),
        ],
        spacing=29,
        font_size=15,
    )

    # Panels C and D: weak minima across composition.
    panels = [
        (
            "c",
            642.0,
            "Spinodal error",
            "relative error (%)",
            "relative_spinodal_error",
            (-1.6, 0.8),
            [-1.5, -1, -0.5, 0, 0.5],
        ),
        (
            "d",
            1212.0,
            "Weak-period error",
            "relative error (%)",
            "relative_period_error",
            (-2.5, 5.0),
            [-2, 0, 2, 4],
        ),
    ]
    for label, panel_x, title, ylabel, field, panel_ylim, yticks in panels:
        panel_label(svg, panel_x, 426, label)
        svg.text(panel_x + 38, 433, title, size=25, weight=700)
        px0, py0, ppw, pph = panel_x + 64, 478.0, 470.0, 330.0
        draw_axes(
            svg,
            x0=px0,
            y0=py0,
            width=ppw,
            height=pph,
            xlim=(0.1, 0.5),
            ylim=panel_ylim,
            xticks=[0.1, 0.2, 0.3, 0.4, 0.5],
            yticks=yticks,
            xlabel="fA",
            ylabel=ylabel,
        )
        curves = weak_curves(weak_rows, field)
        for model in ("exact_rpa", "ok", "uneyama_doi", "bvk1_bvk2"):
            points = [
                (
                    map_x(f_value, px0, ppw, 0.1, 0.5),
                    map_y(error, py0, pph, *panel_ylim),
                )
                for f_value, error in curves[model]
            ]
            svg.polyline(
                points,
                stroke=model_style[model][0],
                dash=model_style[model][1],
                width=3,
            )
        svg.line(
            px0,
            map_y(0, py0, pph, *panel_ylim),
            px0 + ppw,
            map_y(0, py0, pph, *panel_ylim),
            stroke="#727f8c",
            width=1.2,
            dash="5 5",
        )

    draw_legend(
        svg,
        711,
        895,
        [
            ("Full RPA (BURP-TI): exact by construction", BLUE, ""),
            ("Literal nonlinear UD Hessian", ORANGE, "3 5"),
            ("Reduced QCE/AQCE vertex", GREEN, "10 6"),
            ("OK: matched RPA spinodal", PURPLE, "3 6"),
        ],
        spacing=28,
        font_size=15,
    )
    svg.rect(1215, 864, 490, 116, fill=PALE_ORANGE, stroke=ORANGE, radius=12, stroke_width=1.8)
    svg.text(1240, 900, "Weak-limit agreement is a necessary test,", size=19, weight=700, fill=ORANGE)
    svg.text(1240, 932, "not a stress-free-period prediction.", size=19, weight=700, fill=ORANGE)
    svg.text(1240, 962, "Finite-amplitude mechanics are tested in Figures 2–6.", size=15, fill=MID)

    svg.text(
        1758,
        1022,
        "Frozen vertex and weak-response tables; no fitted plotting data",
        size=13,
        anchor="end",
        fill=MID,
    )
    svg.finish(output)


def profile_points(
    profile_rows: list[dict[str, str]], model: str
) -> list[tuple[float, float]]:
    return sorted(
        (float(row["s"]), float(row["phi_a"]))
        for row in profile_rows
        if row["case_id"] == "f0.5_chiN20" and row["model"] == model
    )


def write_toc(
    output: Path,
    kernel_rows: list[dict[str, str]],
    profile_rows: list[dict[str, str]],
    summary_rows: list[dict[str, str]],
    boundary_rows: list[dict[str, str]],
) -> None:
    svg = SVG(
        975,
        525,
        physical_width="3.25in",
        physical_height="1.75in",
        title="AQCE from weak response to phase diagram topology",
        description=(
            "Macromolecules table-of-contents graphic using frozen Gaussian "
            "response, lamellar profile, and accepted phase-boundary data."
        ),
    )
    svg.add(
        "<metadata>"
        "schema=bvk2-macromolecules-toc-v1;"
        "physical_size=3.25in_x_1.75in;viewBox=975x525;"
        "kernel_fA=0.35;profile=f0.5_chiN20;accepted_roots=101"
        "</metadata>"
    )
    svg.add(
        '<rect x="12" y="12" width="951" height="501" rx="24" '
        'fill="#fbfcfe" stroke="#cad4de" stroke-width="2"/>'
    )
    svg.text(42, 60, "AQCE", size=36, weight=700, fill=GREEN)
    svg.text(160, 58, "response", size=21, weight=700, fill=BLUE)
    svg.line(250, 52, 306, 52, stroke=MID, width=2.5, marker="arrow")
    svg.text(329, 58, "period", size=21, weight=700, fill=ORANGE)
    svg.line(413, 52, 469, 52, stroke=MID, width=2.5, marker="arrow")
    svg.text(492, 58, "phase diagram", size=21, weight=700, fill=TEAL)

    # Left: actual normalized Gaussian response at fA = 0.35.
    x0, y0, pw, ph = 52.0, 115.0, 245.0, 285.0
    svg.rect(31, 87, 286, 347, fill="#ffffff", stroke="#d5dee7", radius=16)
    svg.text(52, 111, "weak response", size=17, weight=700, fill=BLUE)
    for y_value in (0.25, 0.5, 0.75, 1.0):
        yy = map_y(y_value, y0, ph, 0.0, 1.05)
        svg.line(x0, yy, x0 + pw, yy, stroke=GRID, width=1)
    svg.line(x0, y0 + ph, x0 + pw, y0 + ph, width=1.4)
    svg.line(x0, y0, x0, y0 + ph, width=1.4)
    kernel = select_kernel(kernel_rows, 0.35)
    for model, color, dash in (
        ("exact_rpa", BLUE, ""),
        ("bvk1_bvk2", GREEN, "8 5"),
    ):
        points = [
            (
                map_x(q, x0, pw, 0.5, 3.6),
                map_y(value, y0, ph, 0.0, 1.05),
            )
            for q, value in kernel[model]
        ]
        svg.polyline(points, stroke=color, width=3.2, dash=dash)
    svg.text(x0 + pw / 2, y0 + ph + 27, "q Rg", size=16, anchor="middle")
    svg.text(68, 461, "RPA", size=14, weight=700, fill=BLUE)
    svg.line(105, 456, 143, 456, stroke=BLUE, width=3)
    svg.text(160, 461, "QCE/AQCE quadratic", size=14, weight=700, fill=GREEN)

    svg.line(325, 260, 360, 260, stroke=ORANGE, width=3, marker="arrow")

    # Center: frozen SCFT and AQCE own-cell lamellar profiles.
    svg.rect(374, 87, 275, 347, fill="#ffffff", stroke="#d5dee7", radius=16)
    svg.text(395, 111, "stress-free lamella", size=17, weight=700, fill=ORANGE)
    px0, py0, ppw, pph = 396.0, 135.0, 230.0, 230.0
    svg.line(px0, py0 + pph, px0 + ppw, py0 + pph, width=1.4)
    svg.line(px0, py0, px0, py0 + pph, width=1.4)
    for model, color, dash in (("scft", INK, "7 5"), ("bvk2", GREEN, "")):
        points = [
            (
                map_x(s, px0, ppw, 0, 1),
                map_y(phi, py0, pph, 0, 1),
            )
            for s, phi in profile_points(profile_rows, model)[::2]
        ]
        svg.polyline(points, stroke=color, width=3, dash=dash)
    summary = {
        row["model"]: row
        for row in summary_rows
        if row["case_id"] == "f0.5_chiN20" and row["model"] in {"scft", "bvk2"}
    }
    scft_period = float(summary["scft"]["period_rg"])
    bvk2_period = float(summary["bvk2"]["period_rg"])
    svg.text(511, 392, f"L/Rg = {bvk2_period:.3f}", size=18, weight=700, anchor="middle", fill=GREEN)
    svg.text(511, 416, f"SCFT {scft_period:.3f}", size=14, anchor="middle", fill=MID)
    svg.text(395, 461, "AQCE", size=14, weight=700, fill=GREEN)
    svg.line(441, 456, 479, 456, stroke=GREEN, width=3)
    svg.text(495, 461, "SCFT", size=14, weight=700)
    svg.line(541, 456, 579, 456, stroke=INK, width=3, dash="7 5")

    svg.line(659, 260, 694, 260, stroke=TEAL, width=3, marker="arrow")

    # Right: accepted left-half phase-boundary roots.
    svg.rect(708, 87, 236, 347, fill="#ffffff", stroke="#d5dee7", radius=16)
    svg.text(728, 111, "phase boundaries", size=17, weight=700, fill=TEAL)
    bx0, by0, bpw, bph = 737.0, 135.0, 177.0, 240.0
    svg.line(bx0, by0 + bph, bx0 + bpw, by0 + bph, width=1.4)
    svg.line(bx0, by0, bx0, by0 + bph, width=1.4)
    boundary_style = {
        "S/DIS": (BLUE, ""),
        "S/C": (RED, ""),
        "G/C": (PURPLE, ""),
        "G/L": (ORANGE, ""),
    }
    for transition in ("S/DIS", "S/C", "G/C", "G/L"):
        rows = sorted(
            (
                float(row["fA"]),
                float(row["chiN"]),
            )
            for row in boundary_rows
            if row["transition"] == transition
        )
        points = [
            (
                map_x(f_value, bx0, bpw, 0.05, 0.5),
                map_y(chi_n, by0, bph, 10, 60),
            )
            for f_value, chi_n in rows
        ]
        svg.polyline(points, stroke=boundary_style[transition][0], width=2.6)
    svg.text(bx0 + bpw / 2, by0 + bph + 27, "fA", size=16, anchor="middle")
    svg.add(
        f'<text transform="translate({f3(bx0 - 22)} {f3(by0 + bph / 2)}) '
        'rotate(-90)" font-size="16" text-anchor="middle">χN</text>'
    )
    legend_positions = [
        ("S/DIS", BLUE, 718),
        ("S/C", RED, 774),
        ("G/C", PURPLE, 821),
        ("G/L", ORANGE, 868),
    ]
    for label, color, xx in legend_positions:
        svg.line(xx, 456, xx + 23, 456, stroke=color, width=3)
        svg.text(xx + 27, 461, label, size=12, weight=700, fill=color)

    svg.text(
        946,
        497,
        "Exact weak response • accurate period • full AB phase diagram",
        size=15,
        weight=700,
        anchor="end",
        fill=MID,
    )
    svg.finish(output)


def validate_inputs(
    kernel_rows: list[dict[str, str]],
    weak_rows: list[dict[str, str]],
    profile_rows: list[dict[str, str]],
    summary_rows: list[dict[str, str]],
    boundary_rows: list[dict[str, str]],
) -> None:
    assert len(kernel_rows) >= 3000
    assert len(weak_rows) == 324
    assert Counter(row["model"] for row in weak_rows) == {
        "exact_rpa": 81,
        "uneyama_doi": 81,
        "bvk1_bvk2": 81,
        "ok": 81,
    }
    for model in ("exact_rpa", "uneyama_doi", "bvk1_bvk2", "ok"):
        assert any(
            row["model"] == model and abs(float(row["f"]) - 0.35) < 1e-12
            for row in kernel_rows
        )
    assert len(profile_points(profile_rows, "scft")) == 256
    assert len(profile_points(profile_rows, "bvk2")) == 256
    assert any(
        row["case_id"] == "f0.5_chiN20" and row["model"] == "scft"
        for row in summary_rows
    )
    assert any(
        row["case_id"] == "f0.5_chiN20" and row["model"] == "bvk2"
        for row in summary_rows
    )
    assert len(boundary_rows) == 101
    assert all(row["status"] == "accepted" for row in boundary_rows)
    assert Counter(row["transition"] for row in boundary_rows) == {
        "S/C": 21,
        "S/DIS": 31,
        "G/C": 23,
        "G/L": 26,
    }


def write_manifest(output_dir: Path, sources: list[tuple[Path, str, int]]) -> None:
    path = output_dir / "source_manifest.csv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["source", "role", "data_rows", "sha256"])
        for source, role, count in sources:
            writer.writerow(
                [source.relative_to(PROJECT), role, count, sha256(source)]
            )


def write_readme(output_dir: Path) -> None:
    readme = """# QCE/AQCE model-hierarchy and TOC graphics

This directory contains deterministic publication graphics assembled only from
frozen CSV artifacts. The generator performs no SCFT, field, cell, or
phase-boundary solve.

## Outputs

- `bvk2_model_hierarchy_weak_response.svg`: final main Figure 1, 1800 × 1050
  SVG user units. Panel (a) compares the chi-independent parts of the
  homogeneous second-order composition vertices; panels (b) and (c) show the composition-dependent spinodal and
  weak-period errors.
- `bvk2_model_hierarchy_weak_response.png`: 1800 × 1050 inspection preview.
- `bvk2_macromolecules_toc.svg`: **3.25 in × 1.75 in**, with a
  **234 × 126 pt** viewBox. Paired color-segmented curved AB and symmetric ABA
  chains identify the architectures without a title or explanatory captions.
  The larger left panel retains
  the Figure 6 AQCE phase boundaries and SCFT reference curves from Li and Liu,
  *J. Chem. Phys.* **2021**, *154*, 014903. The right panel uses the Figure 7d
  signed ABA period errors relative to SCFT for UD, BURP, QCE, and AQCE at
  `fA=0.5` over the seven conditions `chiN=20,22,25,30,35,40,45`.
  It has no background shading and uses solid curves with small markers,
  including palette-blue QCE and green AQCE. A thin gray horizontal line
  denotes zero error relative to SCFT. Errors are converted to percent without
  further processing. No new calculations or boundary fitting are performed.
- `bvk2_macromolecules_toc.png`: exact 975 × 525, 300 dpi-tagged inspection
  preview. The SVG remains the authoritative publication artifact.
- `source_manifest.csv`: source roles, row counts, and SHA-256 hashes.

## Frozen sources

1. `results/diblock_kernel_comparison/kernel_samples.csv`
2. `results/diblock_weak_response_map/weak_response_map.csv`
3. `results/bvk2_aba_comprehensive_validation/benchmark_summary.csv`
4. `results/bvk2_aba_fixed_stiffness_validation/summary.csv`
5. `results/bvk2_publication_phase_diagram/phase_diagram_plot_data.csv`
6. `results/bvk2_publication_phase_diagram/scft_reference_boundaries.csv`

The Figure 1 labels identify the chi-independent composition vertex associated
with each finite-amplitude functional. The QCE and AQCE vertices both equal
`Gamma_psi` because the differentiable AQCE correction begins at fourth order in
a weak composition modulation. The UD and QCE/AQCE vertices have the same
nonlocal-to-gradient coefficient ratio and therefore share the weak-period
error, while their different local curvatures produce distinct spinodal errors.

## Reproduce

From the `MonteCarlo` project root:

```bash
uv run python scripts/write_bvk2_figure1_and_toc.py
uv run python test/test_write_bvk2_figure1_and_toc.py
```

To regenerate only the TOC (including the review PDF figure):

```bash
uv run python -c 'from scripts.render_macromolecules_figures import render_toc; render_toc()'
```

Both SVGs are byte-deterministic for unchanged input tables.
The checked-in PNG previews are also byte-deterministic in the recorded
LibreOffice/Pillow environment.
"""
    (output_dir / "README.md").write_text(readme, encoding="utf-8")


def render_png(svg_path: Path, png_path: Path, size: tuple[int, int]) -> None:
    """Rasterize one SVG, then normalize dimensions and PNG encoding."""
    libreoffice = shutil.which("libreoffice")
    if libreoffice is None:
        raise RuntimeError("LibreOffice is required to render PNG previews")
    width, height = size
    with tempfile.TemporaryDirectory() as temporary:
        temporary_path = Path(temporary)
        render_svg = temporary_path / svg_path.name
        svg_text = svg_path.read_text(encoding="utf-8")
        if 'width="3.25in" height="1.75in"' in svg_text:
            svg_text = svg_text.replace(
                'width="3.25in" height="1.75in"',
                f'width="{width}" height="{height}"',
                1,
            )
        render_svg.write_text(svg_text, encoding="utf-8")
        home = temporary_path / "home"
        home.mkdir()
        completed = subprocess.run(
            [
                libreoffice,
                "--headless",
                "--convert-to",
                "png",
                "--outdir",
                str(temporary_path),
                str(render_svg),
            ],
            env={"HOME": str(home), "PATH": str(Path(libreoffice).parent)},
            capture_output=True,
            text=True,
            timeout=60,
            check=False,
        )
        rendered = temporary_path / f"{render_svg.stem}.png"
        if completed.returncode != 0 or not rendered.exists():
            raise RuntimeError(
                "LibreOffice SVG rasterization failed: "
                + completed.stderr.strip()
            )
        with Image.open(rendered) as image:
            normalized = image.convert("RGB")
            if normalized.size != size:
                resampling = getattr(Image, "Resampling", Image)
                normalized = normalized.resize(size, resampling.LANCZOS)
            normalized.save(
                png_path,
                format="PNG",
                compress_level=9,
                optimize=False,
                dpi=(300, 300),
            )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    output_dir = args.output_dir.resolve()
    if output_dir != DEFAULT_OUTPUT.resolve():
        raise ValueError(
            "custom --output-dir is not supported by the authoritative renderer"
        )

    from render_macromolecules_figures import (
        _write_model_hierarchy_source_manifest,
        render_kernel_and_hierarchy,
    )

    render_kernel_and_hierarchy()
    _write_model_hierarchy_source_manifest()
    render_png(
        output_dir / "bvk2_model_hierarchy_weak_response.svg",
        output_dir / "bvk2_model_hierarchy_weak_response.png",
        (1800, 1050),
    )
    render_png(
        output_dir / "bvk2_macromolecules_toc.svg",
        output_dir / "bvk2_macromolecules_toc.png",
        (975, 525),
    )
    print("wrote authoritative three-panel Figure 1 and Macromolecules TOC")


if __name__ == "__main__":
    main()
