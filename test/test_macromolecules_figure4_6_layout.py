#!/usr/bin/env python3
"""Focused publication-layout regressions for Macromolecules Figures 4 and 6."""

from __future__ import annotations

import re
import sys
from pathlib import Path
from unittest.mock import patch

import numpy as np
import matplotlib.pyplot as plt


PROJECT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT))

import scripts.render_macromolecules_figures as figures  # noqa: E402


def test_main_and_aba_figures_share_model_colors() -> None:
    from matplotlib.colors import to_hex
    from scripts.render_macromolecules_figure8 import COLORS as aba_colors

    expected = {
        "scft": "#7f7f7f", "uneyama_doi": "#2ca02c",
        "burp_ti": "#9467bd", "bvk2_fixed": "#ff7f0e", "bvk2": "#d62728",
    }
    for model, color in expected.items():
        assert to_hex(figures.COLORS[model]) == to_hex(aba_colors[model]) == color
    assert to_hex(figures.COLORS["liu2019_opf"]) == "#1f77b4"
    model_colors = {*expected.values(), to_hex(figures.COLORS["liu2019_opf"])}
    assert len(model_colors) == 6
    assert model_colors <= {to_hex(color) for color in figures.ACS_COLORS[:10]}


def test_figure4_rightmost_profile_tick_stays_inside_canvas() -> None:
    svg = (
        PROJECT
        / "results/liu2019_stress_free_period_map"
        / "liu2019_stress_free_profile_comparison.svg"
    ).read_text(encoding="utf-8")
    root = re.search(r'<svg [^>]*width="([0-9.]+)pt"', svg)
    assert root is not None
    canvas_width = float(root.group(1))
    tick_positions = [
        float(value)
        for value in re.findall(
            r"<!-- 1\.0 -->\s*<g transform=\"translate\(([0-9.]+) ", svg
        )
    ]
    assert len(tick_positions) == 4
    # The rightmost label is about 11 pt wide at the rendered 8 pt font.
    assert max(tick_positions) + 11.0 < canvas_width


def test_figure6_labels_and_insets_reuse_identical_plot_layers() -> None:
    captured = {}

    def inspect_figure(fig, _output, **_kwargs) -> None:
        main_ax, fcc_ax, o70_ax = fig.axes
        reference_groups = figures._groups(
            figures.rows(figures.FIGURE6_SCFT_DATA),
            ("transition", "side", "curve_id"),
        )
        displayed = {}
        for line, ((transition, side, _), group) in zip(
            main_ax.lines, reference_groups.items()
        ):
            raw = np.array([
                (float(row["fA"]), float(row["chiN"]))
                for row in sorted(group, key=lambda row: int(row["point_index"]))
            ])
            curve = np.column_stack((line.get_xdata(), line.get_ydata()))
            if transition in {"C/O70", "G/O70", "O70/L", "G/L"}:
                assert len(curve) == 5 * (len(raw) - 1) + 1
                np.testing.assert_array_equal(curve[[0, -1]], raw[[0, -1]])
                displacement = np.abs(curve[::5] - raw)
                assert displacement[:, 0].max() <= 1e-3
                assert displacement[:, 1].max() <= 0.03
                if transition == "G/L":
                    np.testing.assert_array_equal(curve[::5, 1], raw[:, 1])
                assert np.all(np.diff(curve[:, 1]) >= 0)
                assert np.all(np.diff(curve[:, 0]) * (-1 if side == "left" else 1) >= 0)
                displayed[transition, side] = curve
            else:
                np.testing.assert_array_equal(curve, raw)
        for transition in ("C/O70", "G/O70", "O70/L", "G/L"):
            left, right = displayed[transition, "left"], displayed[transition, "right"]
            np.testing.assert_allclose(left[:, 0], 1 - right[:, 0], atol=2e-12, rtol=0)
            np.testing.assert_array_equal(left[:, 1], right[:, 1])
        # The two pocket edges retain their ordering, including near junctions.
        right_edge = displayed["O70/L", "left"]
        for transition in ("C/O70", "G/O70"):
            left_edge = displayed[transition, "left"]
            chi = np.linspace(left_edge[0, 1], left_edge[-1, 1], 501)[1:-1]
            assert np.all(
                np.interp(chi, right_edge[:, 1], right_edge[:, 0])
                > np.interp(chi, left_edge[:, 1], left_edge[:, 0])
            )
        # Exclude the main-only SCFT legend handle and f_A=0.5 guide.
        main_layers = main_ax.lines[:-2]
        # FCC/BCC must use only the original vertices, joined linearly.
        triple = figures.rows(figures.FIGURE6_TRIPLE_POINT_SOURCES[0][0])[0]
        knots = sorted(
            [(float(triple["fA_estimate"]), float(triple["chiN_estimate"]))]
            + [(float(row["fA"]), float(row["chiN"]))
               for row in figures.rows(figures.FIGURE6_PLOT_DATA)
               if row["transition"] == "FCC/BCC"],
            key=lambda point: point[1],
        )
        fcc_lines = [line for line in main_layers
                     if line.get_color() == figures.FIGURE6_COLORS["FCC/BCC"]]
        assert len(fcc_lines) == 2
        for line, expected in zip(fcc_lines, (knots[::-1], [(1-f, chi) for f, chi in knots])):
            np.testing.assert_array_equal(
                np.column_stack((line.get_xdata(), line.get_ydata())), expected
            )
            assert line.get_linestyle() == "-"
        assert len(main_layers) == len(fcc_ax.lines) == len(o70_ax.lines)
        for axis in fig.axes:
            scft_lines = [line for line in axis.lines if line.get_linestyle() == "--"]
            assert scft_lines
            assert all(line.get_color() == figures.COLORS["scft"] for line in scft_lines)
            markers = [line for line in axis.lines if line.get_marker() != "None"]
            assert len(markers) == len(figures.FIGURE6_TRIPLE_POINT_SOURCES)
            for marker, (path, f_key, chi_key, _) in zip(
                markers, figures.FIGURE6_TRIPLE_POINT_SOURCES
            ):
                row = figures.rows(path)[0]
                f_a, chi_n = float(row[f_key]), float(row[chi_key])
                np.testing.assert_array_equal(marker.get_xdata(), [f_a, 1 - f_a])
                np.testing.assert_array_equal(marker.get_ydata(), [chi_n, chi_n])
                assert marker.get_marker() == "o"
                assert marker.get_markerfacecolor() == "#808080"
                assert marker.get_markeredgewidth() == 0
        assert set(main_ax.get_legend_handles_labels()[1]) == {
            "SCFT", "triple point"
        }
        # Only the vertical composition-symmetry guide remains dotted.
        dotted = [line for line in main_ax.lines if line.get_linestyle() == ":"]
        assert len(dotted) == 1
        assert all(value == 0.5 for value in dotted[0].get_xdata())
        for axis in (fcc_ax, o70_ax):
            assert all(line.get_linestyle() != ":" for line in axis.lines)
        fig.canvas.draw()
        legend_box = main_ax.get_legend().get_window_extent()
        assert main_ax.get_window_extent().contains(legend_box.x0, legend_box.y0)
        assert main_ax.get_window_extent().contains(legend_box.x1, legend_box.y1)
        assert legend_box.x0 > main_ax.get_window_extent().x0 + main_ax.get_window_extent().width / 2
        assert legend_box.y1 < main_ax.get_window_extent().y0 + main_ax.get_window_extent().height / 2
        for main_line, fcc_line, o70_line in zip(
            main_layers, fcc_ax.lines, o70_ax.lines
        ):
            np.testing.assert_array_equal(main_line.get_xdata(), fcc_line.get_xdata())
            np.testing.assert_array_equal(main_line.get_ydata(), fcc_line.get_ydata())
            np.testing.assert_array_equal(main_line.get_xdata(), o70_line.get_xdata())
            np.testing.assert_array_equal(main_line.get_ydata(), o70_line.get_ydata())
        captured["labels"] = {
            text.get_text() for axis in fig.axes for text in axis.texts
        }
        captured["titles"] = {axis.get_title() for axis in fig.axes}
        captured["main_positions"] = {
            text.get_text(): text.get_position() for text in main_ax.texts
        }
        captured["main_anchors"] = {
            text.get_text(): text.xy for text in main_ax.texts if hasattr(text, "xy")
        }
        captured["fcc_anchor"] = next(
            text.xy for text in fcc_ax.texts if text.get_text() == "FCC"
        )
        captured["o70_position"] = next(
            text.get_position()
            for text in o70_ax.texts
            if text.get_text() == r"$O^{70}$"
        )
        captured["o70_anchor"] = next(
            text.xy for text in o70_ax.texts if text.get_text() == r"$O^{70}$"
        )
        plt.close(fig)

    original_geometry = figures.phase_boundary_geometry
    geometry_calls = []

    def record_geometry(*args, **kwargs):
        geometry_calls.append(args[1])
        return original_geometry(*args, **kwargs)

    with (
        patch.object(figures, "phase_boundary_geometry", side_effect=record_geometry),
        patch.object(figures, "save", side_effect=inspect_figure),
    ):
        figures.render_phase_diagram()

    assert geometry_calls == list(figures.FIGURE6_TRANSITIONS)
    assert {"FCC", "BCC", "CYL", "GYR", "LAM", "DIS"} <= captured["labels"]
    assert r"$O^{70}$" in captured["labels"]
    assert captured["titles"] == {"", "FCC pocket", r"$O^{70}$ pocket"}

    boundary_rows = figures.rows(figures.FIGURE6_PLOT_DATA)
    anchors = figures.figure6_region_anchors(boundary_rows)
    for label, chi_n, left_transition, right_transition in (
        ("FCC", 37.0, "FCC/DIS", "FCC/BCC"),
        ("BCC", 37.0, "FCC/BCC", "S/C"),
        ("CYL", 33.0, "S/C", "G/C"),
        ("GYR", 29.0, "G/C", "G/L"),
        ("O70", 13.5, "G/O70", "O70/L"),
    ):
        anchor_f_a, anchor_chi_n = anchors[label]
        left = figures._boundary_f_a_at_chi(
            boundary_rows, left_transition, chi_n
        )
        right = figures._boundary_f_a_at_chi(
            boundary_rows, right_transition, chi_n
        )
        assert anchor_chi_n == chi_n
        assert left < anchor_f_a < right
    lam_f_a, lam_chi_n = anchors["LAM"]
    lam_left = figures._boundary_f_a_at_chi(boundary_rows, "G/L", lam_chi_n)
    assert lam_left < lam_f_a < 1.0 - lam_left

    for label in ("BCC", "CYL", "LAM", "DIS"):
        assert captured["main_positions"][label] == anchors[label]
    for label in ("FCC", "GYR"):
        assert captured["main_anchors"][label] == anchors[label]
    assert captured["fcc_anchor"] == anchors["FCC"]
    assert captured["o70_anchor"] == anchors["O70"]
    label_f, label_chi = captured["o70_position"]
    assert label_f > figures._boundary_f_a_at_chi(boundary_rows, "O70/L", label_chi)
