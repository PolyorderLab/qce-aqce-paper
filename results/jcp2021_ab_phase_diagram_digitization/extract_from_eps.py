#!/usr/bin/env python3
"""Extract the AB SCFT phase boundaries from the vector paths in JCP 154, 014903.

The supplied EPS is a DOS EPS binary containing an OriginLab PostScript stream and
an embedded preview.  The phase boundaries are recovered from the PostScript
paths, whereas the verification overlay is calibrated independently from a
Ghostscript rasterization of the complete source file.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import subprocess
import tempfile
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFont


MONTE_CARLO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_SOURCE = MONTE_CARLO_ROOT.parent / "references" / "AB_phase_diagram_fig1_jcp2021.eps"
DEFAULT_OUTPUT = Path(__file__).resolve().parent

# Native Origin coordinates of the top-panel axes.
SOURCE_X_LEFT = 954.0
SOURCE_X_RIGHT = 2404.0
SOURCE_Y_TOP = 98.0
SOURCE_Y_BOTTOM = 1270.0
SOURCE_X_CENTER = 1679.0
SOURCE_Y_CRITICAL = 1023.0

THEORY = "SCFT"
SOURCE_LABEL = "Li and Liu, J. Chem. Phys. 154, 014903 (2021), Fig. 1(a)"


@dataclass(frozen=True)
class Segment:
    transition: str
    side: str
    source_path_id: int
    points: tuple[tuple[float, float, bool], ...]


PATH_ASSIGNMENTS = {
    0: ("DIS/S_cp", "left", "reverse"),
    1: ("S_cp/S", "left", "reverse"),
    3: ("S_cp/S", "right", "reverse"),
    4: ("DIS/S_cp", "right", "reverse"),
    6: ("C/G", "left", "forward"),
    8: ("G/O70", "left", "reverse"),
    9: ("G/L", "left", "reverse"),
    11: ("G/L", "right", "reverse"),
    12: ("G/O70", "right", "reverse"),
    13: ("C/G", "right", "forward"),
}

SPLIT_PATH_ASSIGNMENTS = {
    2: "DIS/S",
    5: "S/C",
    7: "C/O70",
    10: "O70/L",
}

COLORS = {
    "DIS/S_cp": "#d62728",
    "S_cp/S": "#ff7f0e",
    "DIS/S": "#8c564b",
    "S/C": "#2ca02c",
    "C/G": "#1f77b4",
    "C/O70": "#17becf",
    "G/O70": "#9467bd",
    "G/L": "#e377c2",
    "O70/L": "#bcbd22",
}


def extract_postscript(eps_path: Path) -> str:
    raw = eps_path.read_bytes()
    start = raw.find(b"%!PS")
    if start < 0:
        raise ValueError(f"No PostScript stream found in {eps_path}")
    end = raw.find(b"%%EOF", start)
    if end < 0:
        raise ValueError(f"No PostScript EOF marker found in {eps_path}")
    return raw[start : end + len(b"%%EOF")].decode("latin-1").replace("\r\n", "\n")


def parse_top_panel_paths(postscript: str) -> list[list[tuple[int, int]]]:
    pattern = re.compile(
        r"pathproc 954 98 1451 1172 np rectpath\n"
        r"/eocl cland\nnp exec\n"
        r"(.*?)"
        r"(?=\ngr\ngs\npathproc 954 (?:98|1561))",
        re.DOTALL,
    )
    blocks = pattern.findall(postscript)
    if len(blocks) != 16:
        raise ValueError(f"Expected 16 top-panel clipped paths; found {len(blocks)}")

    paths: list[list[tuple[int, int]]] = []
    for block in blocks:
        raw_points = [
            (int(x), int(y))
            for x, y, _operator in re.findall(r"(-?\d+) (-?\d+) ([ml])(?:\n|$)", block)
        ]
        points: list[tuple[int, int]] = []
        for point in raw_points:
            if not points or point != points[-1]:
                points.append(point)
        paths.append(points)

    expected_lengths = [100, 100, 100, 100, 100, 100, 100, 100, 63, 100, 100, 100, 63, 100, 8, 0]
    actual_lengths = [len(path) for path in paths]
    if actual_lengths != expected_lengths:
        raise ValueError(f"Unexpected path structure: {actual_lengths}")
    return paths


def _with_insert_flag(points: list[tuple[int, int]]) -> list[tuple[float, float, bool]]:
    return [(float(x), float(y), False) for x, y in points]


def split_at_critical(path: list[tuple[int, int]]) -> tuple[list[tuple[float, float, bool]], list[tuple[float, float, bool]]]:
    left = _with_insert_flag([point for point in path if point[0] < SOURCE_X_CENTER])
    right = _with_insert_flag([point for point in path if point[0] > SOURCE_X_CENTER])
    center = (SOURCE_X_CENTER, SOURCE_Y_CRITICAL, True)

    # Each U-shaped boundary is stored from its left junction to its right
    # junction.  Standardize both sides from the low-chi node toward high chi.
    left = [center, *reversed(left)]
    right = [center, *right]
    return left, right


def build_segments(paths: list[list[tuple[int, int]]]) -> list[Segment]:
    segments: list[Segment] = []
    for path_id, (transition, side, orientation) in PATH_ASSIGNMENTS.items():
        points = _with_insert_flag(paths[path_id])
        if orientation == "reverse":
            points.reverse()
        segments.append(Segment(transition, side, path_id, tuple(points)))

    for path_id, transition in SPLIT_PATH_ASSIGNMENTS.items():
        left, right = split_at_critical(paths[path_id])
        segments.append(Segment(transition, "left", path_id, tuple(left)))
        segments.append(Segment(transition, "right", path_id, tuple(right)))

    transition_order = {name: index for index, name in enumerate(COLORS)}
    return sorted(segments, key=lambda segment: (transition_order[segment.transition], segment.side))


def source_to_data(x: float, y: float) -> tuple[float, float]:
    f_a = (x - SOURCE_X_LEFT) / (SOURCE_X_RIGHT - SOURCE_X_LEFT)
    chi_n = 50.0 * (SOURCE_Y_BOTTOM - y) / (SOURCE_Y_BOTTOM - SOURCE_Y_TOP)
    return f_a, chi_n


def write_boundaries_csv(segments: list[Segment], output_path: Path) -> None:
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "theory",
                "source",
                "transition",
                "side",
                "source_path_id",
                "point_index",
                "f_A",
                "chiN",
                "source_x",
                "source_y",
                "inserted_critical_point",
            ]
        )
        for segment in segments:
            for point_index, (x, y, inserted) in enumerate(segment.points):
                f_a, chi_n = source_to_data(x, y)
                writer.writerow(
                    [
                        THEORY,
                        SOURCE_LABEL,
                        segment.transition,
                        segment.side,
                        segment.source_path_id,
                        point_index,
                        f"{f_a:.12f}",
                        f"{chi_n:.12f}",
                        f"{x:.3f}",
                        f"{y:.3f}",
                        str(inserted).lower(),
                    ]
                )


def topology_nodes() -> list[dict[str, object]]:
    nodes = [
        ("critical", "center", 1679, 1023, "DIS/S;S/C;C/O70;O70/L"),
        ("sphere_triple", "left", 1298, 854, "DIS/S_cp;S_cp/S;DIS/S"),
        ("sphere_triple", "right", 2060, 854, "DIS/S_cp;S_cp/S;DIS/S"),
        ("O70_G_L_triple", "left", 1557, 953, "G/L;G/O70;O70/L"),
        ("O70_G_L_triple", "right", 1801, 953, "G/L;G/O70;O70/L"),
        ("O70_C_G_triple", "left", 1585, 999, "C/G;G/O70;C/O70"),
        ("O70_C_G_triple", "right", 1773, 999, "C/G;G/O70;C/O70"),
    ]
    rows = []
    for node, side, x, y, transitions in nodes:
        f_a, chi_n = source_to_data(x, y)
        rows.append(
            {
                "node": node,
                "side": side,
                "f_A": f_a,
                "chiN": chi_n,
                "source_x": x,
                "source_y": y,
                "connected_transitions": transitions,
            }
        )
    return rows


def write_topology_csv(rows: list[dict[str, object]], output_path: Path) -> None:
    fieldnames = ["node", "side", "f_A", "chiN", "source_x", "source_y", "connected_transitions"]
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            rendered = dict(row)
            rendered["f_A"] = f"{float(row['f_A']):.12f}"
            rendered["chiN"] = f"{float(row['chiN']):.12f}"
            writer.writerow(rendered)


def rasterize_eps(source: Path, target: Path, dpi: int) -> None:
    subprocess.run(
        [
            "gs",
            "-q",
            "-dSAFER",
            "-dBATCH",
            "-dNOPAUSE",
            "-dEPSCrop",
            "-sDEVICE=pnggray",
            f"-r{dpi}",
            f"-sOutputFile={target}",
            str(source),
        ],
        check=True,
    )


def _clusters(indices: np.ndarray) -> list[np.ndarray]:
    if len(indices) == 0:
        return []
    split_locations = np.where(np.diff(indices) > 1)[0] + 1
    return list(np.split(indices, split_locations))


def detect_top_panel_frame(gray: np.ndarray) -> tuple[float, float, float, float, float | None]:
    height, width = gray.shape
    dark = gray < 64
    search_height = int(0.36 * height)
    row_counts = dark[:search_height].sum(axis=1)
    row_candidates = np.flatnonzero(row_counts > 0.70 * width)
    row_clusters = _clusters(row_candidates)
    if len(row_clusters) < 2:
        raise ValueError("Could not independently locate the top-panel horizontal axes")
    y_top = float(row_clusters[0].mean())
    y_bottom = float(row_clusters[1].mean())
    if y_bottom - y_top < 0.20 * height:
        raise ValueError("Detected top-panel axes are implausibly close")

    y0 = max(0, int(math.floor(y_top)))
    y1 = min(height, int(math.ceil(y_bottom)) + 1)
    col_counts = dark[y0:y1].sum(axis=0)
    col_candidates = np.flatnonzero(col_counts > 0.80 * (y1 - y0))
    col_clusters = _clusters(col_candidates)
    if len(col_clusters) != 2:
        raise ValueError(f"Expected two vertical axis clusters; found {len(col_clusters)}")
    x_left = float(col_clusters[0].mean())
    x_right = float(col_clusters[1].mean())

    next_panel_top = float(row_clusters[2].mean()) if len(row_clusters) > 2 else None
    return x_left, x_right, y_top, y_bottom, next_panel_top


def data_to_raster(
    f_a: float,
    chi_n: float,
    frame: tuple[float, float, float, float, float | None],
) -> tuple[float, float]:
    x_left, x_right, y_top, y_bottom, _ = frame
    x = x_left + f_a * (x_right - x_left)
    y = y_bottom - (chi_n / 50.0) * (y_bottom - y_top)
    return x, y


def nearest_dark_distance(gray: np.ndarray, x: float, y: float, radius: int) -> float:
    x0 = max(0, int(math.floor(x)) - radius)
    x1 = min(gray.shape[1], int(math.ceil(x)) + radius + 1)
    y0 = max(0, int(math.floor(y)) - radius)
    y1 = min(gray.shape[0], int(math.ceil(y)) + radius + 1)
    yy, xx = np.nonzero(gray[y0:y1, x0:x1] < 96)
    if len(xx) == 0:
        return float("inf")
    dx = xx + x0 - x
    dy = yy + y0 - y
    return float(np.sqrt(dx * dx + dy * dy).min())


def mirror_source_distance(segment: Segment, counterpart: Segment) -> float:
    counterpart_xy = np.array([(x, y) for x, y, _ in counterpart.points], dtype=float)
    distances = []
    for x, y, _ in segment.points:
        reflected = np.array([2.0 * SOURCE_X_CENTER - x, y])
        distances.append(float(np.linalg.norm(counterpart_xy - reflected, axis=1).min()))
    return max(distances)


def write_audit_csv(
    segments: list[Segment],
    gray: np.ndarray,
    frame: tuple[float, float, float, float, float | None],
    output_path: Path,
    dpi: int,
) -> list[dict[str, object]]:
    by_key = {(segment.transition, segment.side): segment for segment in segments}
    radius = max(8, int(round(dpi / 25)))
    rows: list[dict[str, object]] = []
    for segment in segments:
        residuals = []
        for x, y, _ in segment.points:
            f_a, chi_n = source_to_data(x, y)
            raster_x, raster_y = data_to_raster(f_a, chi_n, frame)
            residuals.append(nearest_dark_distance(gray, raster_x, raster_y, radius))
        other_side = "right" if segment.side == "left" else "left"
        counterpart = by_key[(segment.transition, other_side)]
        rows.append(
            {
                "transition": segment.transition,
                "side": segment.side,
                "source_path_id": segment.source_path_id,
                "n_points": len(segment.points),
                "raster_residual_median_px": float(np.median(residuals)),
                "raster_residual_p90_px": float(np.percentile(residuals, 90)),
                "raster_residual_max_px": max(residuals),
                "mirror_source_distance_max": mirror_source_distance(segment, counterpart),
                "mirror_fA_equivalent_max": mirror_source_distance(segment, counterpart)
                / (SOURCE_X_RIGHT - SOURCE_X_LEFT),
            }
        )

    fieldnames = list(rows[0])
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            rendered = dict(row)
            for key in (
                "raster_residual_median_px",
                "raster_residual_p90_px",
                "raster_residual_max_px",
                "mirror_source_distance_max",
                "mirror_fA_equivalent_max",
            ):
                rendered[key] = f"{float(row[key]):.8f}"
            writer.writerow(rendered)
    return rows


def draw_overlay(
    gray: np.ndarray,
    segments: list[Segment],
    nodes: list[dict[str, object]],
    frame: tuple[float, float, float, float, float | None],
    output_path: Path,
    dpi: int,
) -> None:
    base = Image.fromarray(gray).convert("RGB")
    draw = ImageDraw.Draw(base)
    line_width = max(3, int(round(dpi / 120)))

    for segment in segments:
        raster_points = []
        for x, y, _ in segment.points:
            f_a, chi_n = source_to_data(x, y)
            raster_points.append(data_to_raster(f_a, chi_n, frame))
        draw.line(raster_points, fill=COLORS[segment.transition], width=line_width, joint="curve")

    node_radius = max(5, int(round(dpi / 60)))
    for row in nodes:
        x, y = data_to_raster(float(row["f_A"]), float(row["chiN"]), frame)
        draw.ellipse(
            (x - node_radius, y - node_radius, x + node_radius, y + node_radius),
            outline="#000000",
            fill="#ffffff",
            width=max(2, line_width // 2),
        )

    _, _, _, y_bottom, next_panel_top = frame
    natural_panel_bottom = int(y_bottom + 0.20 * (y_bottom - frame[2]))
    crop_bottom = min(int(next_panel_top - 10), natural_panel_bottom) if next_panel_top is not None else natural_panel_bottom
    crop_bottom = min(base.height, crop_bottom)
    cropped = base.crop((0, 0, base.width, crop_bottom))

    legend_height = max(150, int(round(dpi * 0.55)))
    canvas = Image.new("RGB", (cropped.width, cropped.height + legend_height), "white")
    canvas.paste(cropped, (0, 0))
    legend = ImageDraw.Draw(canvas)
    font = ImageFont.load_default(size=max(12, int(round(dpi / 22))))
    columns = 3
    column_width = canvas.width / columns
    row_height = legend_height / 3
    for index, (transition, color) in enumerate(COLORS.items()):
        column = index % columns
        row = index // columns
        x0 = int(column * column_width + 0.08 * column_width)
        y0 = int(cropped.height + (row + 0.5) * row_height)
        sample_length = int(0.16 * column_width)
        legend.line((x0, y0, x0 + sample_length, y0), fill=color, width=line_width + 1)
        legend.text((x0 + sample_length + 12, y0), transition, fill="black", font=font, anchor="lm")
    canvas.save(output_path, dpi=(dpi, dpi))


def validate_extraction(segments: list[Segment], audit_rows: list[dict[str, object]]) -> None:
    expected_transitions = set(COLORS)
    actual_transitions = {segment.transition for segment in segments}
    if actual_transitions != expected_transitions:
        raise ValueError(f"Boundary set mismatch: {actual_transitions}")
    counts = defaultdict(set)
    for segment in segments:
        counts[segment.transition].add(segment.side)
    if any(sides != {"left", "right"} for sides in counts.values()):
        raise ValueError(f"Every transition must have left and right branches: {dict(counts)}")

    segment_lookup = {(segment.transition, segment.side): segment for segment in segments}
    for node in topology_nodes():
        node_xy = (float(node["source_x"]), float(node["source_y"]))
        node_sides = ("left", "right") if node["side"] == "center" else (str(node["side"]),)
        for transition in str(node["connected_transitions"]).split(";"):
            for side in node_sides:
                endpoints = {
                    segment_lookup[(transition, side)].points[0][:2],
                    segment_lookup[(transition, side)].points[-1][:2],
                }
                if node_xy not in endpoints:
                    raise ValueError(
                        f"Topology node {node['node']} ({side}) is not an endpoint of {transition}: {endpoints}"
                    )

    max_p90 = max(float(row["raster_residual_p90_px"]) for row in audit_rows)
    if max_p90 > 3.0:
        raise ValueError(f"Raster overlay audit failed: maximum segment p90 residual is {max_p90:.3f} px")

    # The high-chi phase sequence on the A-poor side is an independent topology
    # check on the source-path assignment.
    high_chi = {}
    for transition in ("DIS/S_cp", "S_cp/S", "S/C", "C/G", "G/L"):
        segment = next(item for item in segments if item.transition == transition and item.side == "left")
        x, y, _ = min(segment.points, key=lambda point: abs(source_to_data(point[0], point[1])[1] - 50.0))
        high_chi[transition] = source_to_data(x, y)[0]
    ordered = [high_chi[name] for name in ("DIS/S_cp", "S_cp/S", "S/C", "C/G", "G/L")]
    if ordered != sorted(ordered):
        raise ValueError(f"Incorrect high-chi phase ordering: {high_chi}")


def write_readme(
    output_path: Path,
    source: Path,
    frame: tuple[float, float, float, float, float | None],
    audit_rows: list[dict[str, object]],
) -> None:
    critical_f, critical_chi = source_to_data(SOURCE_X_CENTER, SOURCE_Y_CRITICAL)
    max_p90 = max(float(row["raster_residual_p90_px"]) for row in audit_rows)
    max_residual = max(float(row["raster_residual_max_px"]) for row in audit_rows)
    max_mirror_f = max(float(row["mirror_fA_equivalent_max"]) for row in audit_rows)
    x_left, x_right, y_top, y_bottom, _ = frame
    try:
        displayed_source = source.relative_to(MONTE_CARLO_ROOT.parent)
    except ValueError:
        displayed_source = source
    text = rf"""# JCP 2021 AB SCFT phase-diagram extraction

This directory contains a vector-path extraction of the **top AB-diblock panel**
of Fig. 1 in W. Li and Y.-X. Liu, *J. Chem. Phys.* **154**, 014903
(2021), DOI: 10.1063/5.0037979.  The source file is
`{displayed_source}`.

## What was extracted

The nine coexistence curves are `DIS/S_cp`, `S_cp/S`, `DIS/S`, `S/C`, `C/G`,
`C/O70`, `G/O70`, `G/L`, and `O70/L`.  Together they retain the complete
topology of the small orthorhombic \(O^{{70}}\) (Fddd) pocket, the close-packed
sphere pocket, both outer triple points, the four \(O^{{70}}\)-related triple
points, and the critical point.

The vector axes are

- \(f_A=(x-954)/1450\),
- \(\chi N=50(1270-y)/1172\).

The critical point read from the source geometry is
\((f_A,\chi N)=({critical_f:.6f},{critical_chi:.6f})\).

## Files

- `jcp2021_ab_phase_boundaries.csv`: all source vertices, separated by
  transition and symmetry side.  No curve fitting or smoothing was applied.
- `jcp2021_ab_topology_nodes.csv`: critical and triple-point coordinates with
  their incident coexistence curves.
- `jcp2021_ab_overlay_check.png`: colored vector paths over an independently
  calibrated Ghostscript rasterization of the original EPS.
- `jcp2021_ab_extraction_audit.csv`: per-branch overlay and mirror-symmetry
  residuals.
- `extract_from_eps.py`: complete reproducible extractor and validator.

## Verification

The raster frame was detected from the image rather than inherited from the
PostScript mapping: left/right = {x_left:.2f}/{x_right:.2f} px and top/bottom =
{y_top:.2f}/{y_bottom:.2f} px.  Across all 18 left/right branch segments, the
largest 90th-percentile distance from an extracted point to the corresponding
black source line is {max_p90:.3f} px; the largest individual distance is
{max_residual:.3f} px.  The maximum source-side mirror discrepancy is equivalent
to \(\Delta f_A={max_mirror_f:.6g}\).  The validator also checks the expected
high-\(\chi N\) phase order and the presence of both symmetry branches for every
transition.

Regenerate from the project root with:

```bash
uv run python results/jcp2021_ab_phase_diagram_digitization/extract_from_eps.py
```
"""
    output_path.write_text(text, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--dpi", type=int, default=600)
    args = parser.parse_args()

    source = args.source.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    postscript = extract_postscript(source)
    paths = parse_top_panel_paths(postscript)
    segments = build_segments(paths)
    nodes = topology_nodes()

    write_boundaries_csv(segments, output_dir / "jcp2021_ab_phase_boundaries.csv")
    write_topology_csv(nodes, output_dir / "jcp2021_ab_topology_nodes.csv")

    with tempfile.TemporaryDirectory(prefix="jcp2021_ab_eps_") as temporary_directory:
        raster_path = Path(temporary_directory) / "source.png"
        rasterize_eps(source, raster_path, args.dpi)
        gray = np.asarray(Image.open(raster_path).convert("L"))
    frame = detect_top_panel_frame(gray)
    audit_rows = write_audit_csv(
        segments,
        gray,
        frame,
        output_dir / "jcp2021_ab_extraction_audit.csv",
        args.dpi,
    )
    validate_extraction(segments, audit_rows)
    draw_overlay(
        gray,
        segments,
        nodes,
        frame,
        output_dir / "jcp2021_ab_overlay_check.png",
        args.dpi,
    )
    write_readme(output_dir / "README.md", source, frame, audit_rows)

    print(f"Extracted {len(segments)} symmetry-resolved segments for {len(COLORS)} transitions")
    print(f"Independent raster frame: x=({frame[0]:.2f}, {frame[1]:.2f}), y=({frame[2]:.2f}, {frame[3]:.2f})")
    print(f"Wrote extraction and overlay audit to {output_dir}")


if __name__ == "__main__":
    main()
