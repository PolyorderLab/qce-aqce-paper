#!/usr/bin/env python3
"""Regression checks for Matplotlib Figure 1 and the TOC graphic."""

import csv
import hashlib
from pathlib import Path
import sys
import tempfile
from unittest.mock import patch
from xml.etree import ElementTree as ET

import numpy as np
from PIL import Image


PROJECT = Path(__file__).resolve().parents[1]
RESULT = PROJECT / "results/bvk2_model_hierarchy_graphics"
sys.path.insert(0, str(PROJECT))

import scripts.render_macromolecules_figures as figure_renderer

from scripts.render_macromolecules_figures import (  # noqa: E402
    FIGURE6_PLOT_DATA,
    FIGURE6_SCFT_DATA,
    FIGURE6_SOURCE_SVG,
    FIGURE6_TRANSITIONS,
    SCFT_CRITICAL_CHI,
    _boundary_f_a_at_chi,
    _groups,
    _scft_display_curve,
    _toc_aba_period_errors,
    figure6_toc_sources,
    render_kernel_and_hierarchy,
    render_phase_diagram,
    render_toc,
    rows,
)

for name in (
    "bvk2_model_hierarchy_weak_response.svg",
    "bvk2_macromolecules_toc.svg",
):
    path = RESULT / name
    assert path.is_file(), path
    root = ET.parse(path).getroot()
    assert root.tag.endswith("svg")
    text = path.read_text(encoding="utf-8")
    assert "Matplotlib + mpltex ACS via uv" in text
    assert "<title>" not in text

figure_text = (RESULT / "bvk2_model_hierarchy_weak_response.svg").read_text()
for phrase in (
    "AQCE",
    "spinodal error",
    "weak-period error",
    "Ideal-chain vertex",
    "UD vertex",
    "QCE/AQCE vertex",
    "UD; QCE/AQCE vertices",
    "f_A=0.10",
    "f_A=0.30",
    "f_A=0.50",
    "k^2R_g^2",
):
    assert phrase in figure_text, phrase
assert "BURP-TI" not in figure_text
assert "\\Gamma_\\psi(k)" in figure_text
assert "\\bar{\\Gamma}_\\psi(k)" not in figure_text
assert "Gaussian response" not in figure_text
assert "Liu-2019 OK" not in figure_text
assert "OK" not in figure_text
assert "stroke-dasharray" not in figure_text
for label in ("(a)", "(b)", "(c)"):
    assert label in figure_text, label
assert "(d)" not in figure_text
assert "UD/BVK1" not in figure_text
assert "OK/AQCE" not in figure_text
assert "(fixed; AQCE)" not in figure_text
assert "PRM" not in figure_text
assert "APRM" not in figure_text
kernel_rows = rows(PROJECT / "results/diblock_kernel_comparison/kernel_samples.csv")
assert {0.1, 0.3, 0.5} <= {float(row["f"]) for row in kernel_rows}
weak_rows = rows(PROJECT / "results/diblock_weak_response_map/weak_response_map.csv")
for f_a in sorted({float(row["f"]) for row in weak_rows}):
    shared_period_errors = [
        float(row["relative_period_error"])
        for row in weak_rows
        if float(row["f"]) == f_a
        and row["model"] in {"uneyama_doi", "bvk1_bvk2"}
    ]
    assert len(shared_period_errors) == 2
    np.testing.assert_allclose(
        shared_period_errors,
        shared_period_errors[0],
        rtol=0.0,
        atol=5.0e-8,
    )
renderer = (PROJECT / "scripts/render_macromolecules_figures.py").read_text()
assert 'models = ("uneyama_doi", "bvk1_bvk2", "exact_rpa")' in renderer
assert 'for model in ("uneyama_doi", "bvk1_bvk2")' in renderer
assert "ylim=(1.0e1, 1.0e3)" in renderer
kernel_panel_source = renderer[
    renderer.index("def _kernel_panel") : renderer.index("def _weak_panels")
]
assert 'loc="upper right"' in kernel_panel_source
assert "bbox_to_anchor" not in kernel_panel_source

manifest_path = RESULT / "source_manifest.csv"
with manifest_path.open(newline="", encoding="utf-8") as handle:
    manifest_rows = list(csv.DictReader(handle))
manifest_by_source = {row["source"]: row for row in manifest_rows}
expected_figure1_roles = {
    "results/diblock_kernel_comparison/kernel_samples.csv": "Figure 1 panel (a)",
    "results/diblock_weak_response_map/weak_response_map.csv": (
        "Figure 1 panels (b) and (c)"
    ),
}
for source, role in expected_figure1_roles.items():
    row = manifest_by_source[source]
    source_path = PROJECT / source
    with source_path.open(newline="", encoding="utf-8") as handle:
        data_rows = sum(1 for _row in csv.DictReader(handle))
    assert row["role"] == role
    assert int(row["data_rows"]) == data_rows
    assert row["sha256"] == hashlib.sha256(source_path.read_bytes()).hexdigest()
assert not any("Figure 3 panel" in row["role"] for row in manifest_rows)
expected_toc_roles = {
    "results/bvk2_aba_comprehensive_validation/benchmark_summary.csv": (
        "TOC ABA UD, BURP, and AQCE period errors"
    ),
    "results/bvk2_aba_fixed_stiffness_validation/summary.csv": (
        "TOC ABA QCE period errors"
    ),
}
for source, role in expected_toc_roles.items():
    row = manifest_by_source[source]
    source_path = PROJECT / source
    assert row["role"] == role
    assert row["sha256"] == hashlib.sha256(source_path.read_bytes()).hexdigest()

with tempfile.TemporaryDirectory() as first_dir, tempfile.TemporaryDirectory() as second_dir:
    first = Path(first_dir)
    second = Path(second_dir)
    render_kernel_and_hierarchy(weak_output_dir=first, include_toc=False)
    render_kernel_and_hierarchy(weak_output_dir=second, include_toc=False)
    for name in (
        "diblock_kernel_comparison.svg",
        "bvk2_model_hierarchy_weak_response.svg",
    ):
        first_svg = first / name
        second_svg = second / name
        assert first_svg.is_file()
        assert second_svg.is_file()
        assert first_svg.read_bytes() == second_svg.read_bytes(), name
        assert first_svg.read_bytes() == (
            RESULT / name
            if name.startswith("bvk2_")
            else PROJECT / "results/diblock_kernel_comparison" / name
        ).read_bytes(), name
legacy_writer = (PROJECT / "scripts/write_bvk2_figure1_and_toc.py").read_text()
legacy_main = legacy_writer.split("def main() -> None:", 1)[1]
assert "render_kernel_and_hierarchy()" in legacy_main
assert "write_figure1(" not in legacy_main
assert "phase " + "topology" not in legacy_writer
assert "phase diagram" in legacy_writer

toc_text = (RESULT / "bvk2_macromolecules_toc.svg").read_text()
for phrase in (
    "<!-- AB -->",
    "<!-- ABA -->",
    'id="chain-AB-0"',
    'id="chain-ABA-0.75"',
    "f_A=0.5",
    r"Period error (\%)",
    'id="scft-zero-reference"',
    "<!-- UD -->",
    "<!-- BURP -->",
    "<!-- FCC -->",
    "<!-- BCC -->",
    "<!-- G -->",
    "O^{70}",
    "<!-- $f_A$ -->",
    "<!-- $\\chi N$ -->",
    "<!-- AQCE -->",
    "<!-- SCFT -->",
):
    assert phrase in toc_text, phrase
assert "phase " + "topology" not in toc_text
assert "Density-only AQCE" not in toc_text
assert "captures SCFT structure" not in toc_text
assert "PRM" not in toc_text
assert "APRM" not in toc_text
for removed_annotation in (
    "QCE and AQCE density functionals", "Shared nonlinear construction",
    "Phase behavior", "Architectural transfer",
    "AQCE solid", "SCFT dashed",
):
    assert removed_annotation not in toc_text
toc_root = ET.parse(RESULT / "bvk2_macromolecules_toc.svg").getroot()
assert toc_root.attrib["width"] == "234pt"
assert toc_root.attrib["height"] == "126pt"
with Image.open(RESULT / "bvk2_macromolecules_toc.png") as toc_png:
    assert toc_png.size == (975, 525)

# Verify layout against the actual plotted data, not just label strings.
original_save = figure_renderer.save
def check_toc_layout(fig, *args, **kwargs):
    fig.canvas.draw()
    phase_ax = fig.axes[0]
    assert phase_ax.get_xlim() == (0.0, 1.0)
    np.testing.assert_array_equal(phase_ax.get_xticks(), [0.0, 0.5, 1.0])
    for transition in FIGURE6_TRANSITIONS:
        curves = [line for line in phase_ax.lines
                  if line.get_linestyle() == "-"
                  and line.get_color() == figure_renderer.FIGURE6_COLORS[transition]]
        assert len(curves) == 2
        np.testing.assert_allclose(curves[0].get_xdata(), 1 - curves[1].get_xdata(), atol=1e-15)
        np.testing.assert_array_equal(curves[0].get_ydata(), curves[1].get_ydata())
    scft_groups = _groups(rows(FIGURE6_SCFT_DATA), ("transition", "side", "curve_id"))
    scft_lines = [line for line in phase_ax.lines if line.get_linestyle() == "--"]
    assert len(scft_lines) == len(scft_groups)
    for line, group in zip(scft_lines, scft_groups.values()):
        ordered = sorted(group, key=lambda row: int(row["point_index"]))
        np.testing.assert_array_equal(line.get_xdata(), [float(row["fA"]) for row in ordered])
        np.testing.assert_array_equal(line.get_ydata(), [float(row["chiN"]) for row in ordered])
    period_ax = fig.axes[1]
    assert len(period_ax.collections) == 0  # No background shading.
    expected_colors = {"QCE": figure_renderer.COLORS["bvk2_fixed"],
                       "AQCE": figure_renderer.COLORS["bvk2"],
                       "UD": figure_renderer.COLORS["uneyama_doi"],
                       "BURP": figure_renderer.COLORS["burp_ti"]}
    reference = period_ax.lines[0]
    assert reference.get_gid() == "scft-zero-reference"
    assert np.all(np.asarray(reference.get_ydata()) == 0)
    assert reference.get_color() == figure_renderer.COLORS["scft"]
    assert set(line.get_label() for line in period_ax.lines[1:]) == set(expected_colors)
    for line in period_ax.lines[1:]:
        assert line.get_linestyle() == "-"
        assert line.get_color() == expected_colors[line.get_label()]
        assert len(line.get_xdata()) == 7
        assert line.get_marker() in ("o", "^", "s", "D")
    phase_text = {label.get_text(): label for label in fig.axes[0].texts}
    dis_box = phase_text["DIS"].get_window_extent()
    phase_legend = fig.axes[0].get_legend()
    assert [label.get_text() for label in phase_legend.get_texts()] == ["AQCE", "SCFT"]
    assert [line.get_linestyle() for line in phase_legend.get_lines()] == ["-", "--"]
    assert all(label.get_fontsize() == 5.0 for label in phase_legend.get_texts())
    key_box = phase_legend.get_window_extent()
    assert not dis_box.overlaps(key_box)
    for label in phase_text.values():
        if not hasattr(label, "xy"):
            assert not label.get_window_extent().overlaps(key_box), label.get_text()
    assert phase_ax.get_window_extent().contains(key_box.x0, key_box.y0)
    assert phase_ax.get_window_extent().contains(key_box.x1, key_box.y1)
    chains = [artist for artist in fig.artists
              if (artist.get_gid() or "").startswith("chain-")]
    assert len(chains) == 5
    assert all(np.ptp(chain.get_ydata()) > 0.01 for chain in chains)
    legend_box = period_ax.get_legend().get_window_extent().transformed(
        period_ax.transData.inverted()
    )
    for line in period_ax.lines[1:]:
        # Check the connecting segments as well as the seven data markers.
        x = np.linspace(20, 45, 501)
        y = np.interp(x, *line.get_data())
        covered = ((x >= legend_box.x0) & (x <= legend_box.x1)
                   & (y >= legend_box.y0) & (y <= legend_box.y1))
        assert not np.any(covered), f"TOC legend covers {line.get_label()}"
    for ax in fig.axes:
        assert all(spine.get_visible() for spine in ax.spines.values())
    for label in fig.texts:
        box = label.get_window_extent()
        assert fig.bbox.contains(box.x0, box.y0), label.get_text()
        assert fig.bbox.contains(box.x1, box.y1), label.get_text()
    original_save(fig, *args, **kwargs)

with patch.object(figure_renderer, "save", side_effect=check_toc_layout):
    render_toc()
assert (RESULT / "bvk2_macromolecules_toc.svg").read_text() == toc_text

figure6_rows, figure6_scft, figure6_triples = figure6_toc_sources()
manuscript = (
    PROJECT / "docs/manuscript/macromolecules/manuscript.md"
).read_text(encoding="utf-8")
assert (
    "../../../results/bvk2_publication_phase_diagram/phase_diagram.svg"
    in manuscript
)
assert FIGURE6_SOURCE_SVG == PROJECT / "results/bvk2_publication_phase_diagram/phase_diagram.svg"
assert FIGURE6_PLOT_DATA.name == "phase_diagram_plot_data.csv"
assert FIGURE6_SCFT_DATA.name == "scft_reference_boundaries.csv"
assert len(figure6_rows) == len(rows(FIGURE6_PLOT_DATA))
assert len(figure6_scft) == len(rows(FIGURE6_SCFT_DATA))
assert len(figure6_triples) == 3
assert {row["transition"] for row in figure6_rows} == set(FIGURE6_TRANSITIONS)
figure6_svg = FIGURE6_SOURCE_SVG.read_text(encoding="utf-8")
display_transition = {
    "S/C": "BCC/CYL",
    "S/DIS": "BCC/DIS",
    "G/C": "GYR/CYL",
    "G/L": "GYR/LAM",
    "FCC/BCC": "FCC/BCC",
    "FCC/DIS": "FCC/DIS",
    "C/O70": "CYL/$O^{70}$",
    "G/O70": "GYR/$O^{70}$",
    "O70/L": "$O^{70}$/LAM",
}
for transition in FIGURE6_TRANSITIONS:
    assert display_transition.get(transition, transition) not in figure6_svg, transition
for label in ("FCC", "BCC", "CYL", "GYR", "LAM", "DIS", "SCFT"):
    assert label in figure6_svg, label
assert "triple point" in figure6_svg
assert "RPA" not in figure6_svg
with tempfile.TemporaryDirectory() as directory:
    rerendered = Path(directory) / "phase_diagram.svg"
    render_phase_diagram(output_path=rerendered)
    assert hashlib.sha256(rerendered.read_bytes()).digest() == hashlib.sha256(
        FIGURE6_SOURCE_SVG.read_bytes()
    ).digest()
for chi_n, left_transition, right_transition in (
    (40.0, "FCC/DIS", "FCC/BCC"),
    (13.0, "G/O70", "O70/L"),
):
    left = _boundary_f_a_at_chi(figure6_rows, left_transition, chi_n)
    right = _boundary_f_a_at_chi(figure6_rows, right_transition, chi_n)
    anchor = 0.5 * (left + right)
    assert left < anchor < right
assert "accepted_phase_boundaries.csv" not in renderer[
    renderer.index("figure6_rows, scft, triple_points") :
    renderer.index("toc_png =")
]
toc_phase_renderer = renderer[
    renderer.index("figure6_rows, scft, triple_points") :
    renderer.index("toc_png =")
]
toc_renderer = renderer[renderer.index("def render_toc()") : renderer.index("def _profile_axes")]
assert "width_ratios=(1.5, 1.0)" in toc_renderer
assert '("A", 0.25), ("B", 0.5), ("A", 0.25)' in toc_renderer
assert "D_{\\rm QCE}/D_{\\rm SCFT}" not in toc_text
# Plot the Figure 7d signed relative errors with only the percent conversion.
toc_errors = _toc_aba_period_errors()
assert set(toc_errors) == {"uneyama_doi", "burp_ti", "bvk2", "bvk2_fixed"}
for model, source in (
    ("uneyama_doi", "results/bvk2_aba_comprehensive_validation/benchmark_summary.csv"),
    ("burp_ti", "results/bvk2_aba_comprehensive_validation/benchmark_summary.csv"),
    ("bvk2", "results/bvk2_aba_comprehensive_validation/benchmark_summary.csv"),
    ("bvk2_fixed", "results/bvk2_aba_fixed_stiffness_validation/summary.csv"),
):
    subset = sorted(
        (r for r in rows(source) if r["model"] == model),
        key=lambda r: float(r["chiN"]),
    )
    expected = np.array([(float(r["chiN"]), 100 * float(r["period_relative_error"])) for r in subset])
    np.testing.assert_array_equal(toc_errors[model], expected)
assert "marker=" not in toc_phase_renderer
assert "phase_ax.scatter" not in toc_phase_renderer
assert "chiN_spinodal" not in toc_phase_renderer
assert "xticks=(0.0, 0.5, 1.0)" in toc_phase_renderer
assert "yticks=(10, 50)" in toc_phase_renderer

scft = rows("results/bvk2_publication_phase_diagram/scft_reference_boundaries.csv")
left_curves = {}
for (transition, side, _curve), group in _groups(
    scft, ("transition", "side", "curve_id")
).items():
    if side != "left":
        continue
    f_a, chi_n = _scft_display_curve(group, transition)
    assert np.all(np.diff(chi_n) > 0)
    left_curves[transition] = (f_a, chi_n)

for transition in ("DIS/S", "S/C", "C/O70", "O70/L"):
    f_a, chi_n = left_curves[transition]
    assert np.isclose(chi_n[0], SCFT_CRITICAL_CHI)
    assert np.isclose(f_a[0], 0.5)

ordered_transitions = ("DIS/S_cp", "S_cp/S", "S/C", "C/G", "G/L")
for chi_n in (20.0, 30.0, 40.0):
    values = [
        float(np.interp(chi_n, left_curves[name][1], left_curves[name][0]))
        for name in ordered_transitions
    ]
    assert np.all(np.diff(values) > 0.0005), (chi_n, values)

digitization = PROJECT / "results/jcp2021_ab_phase_diagram_digitization"
with Image.open(digitization / "jcp2021_ab_overlay_check.png") as overlay_check:
    assert overlay_check.width >= 2000
with (digitization / "jcp2021_ab_extraction_audit.csv").open(newline="") as handle:
    audit_rows = list(csv.DictReader(handle))
assert len(audit_rows) == 18
assert max(float(row["raster_residual_p90_px"]) for row in audit_rows) <= 1.0
assert max(float(row["mirror_fA_equivalent_max"]) for row in audit_rows) <= 1e-12

with (digitization / "jcp2021_ab_topology_nodes.csv").open(newline="") as handle:
    topology_nodes = list(csv.DictReader(handle))
assert len(topology_nodes) == 7
critical = next(row for row in topology_nodes if row["node"] == "critical")
assert np.isclose(float(critical["f_A"]), 0.5)
assert np.isclose(float(critical["chiN"]), SCFT_CRITICAL_CHI)
assert set(left_curves) == {
    "DIS/S_cp",
    "S_cp/S",
    "DIS/S",
    "S/C",
    "C/G",
    "C/O70",
    "G/O70",
    "G/L",
    "O70/L",
}

print("Matplotlib Figure 1 and TOC graphic passed")
