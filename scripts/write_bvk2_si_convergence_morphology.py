#!/usr/bin/env python3
"""Curate accepted BVK2 convergence and morphology evidence for the SI.

The script is intentionally solver-free.  It selects four phase endpoints
from the canonical accepted crossing packet, validates their provenance, and
renders diagnostic convergence traces plus representative density sections.
"""

from __future__ import annotations

import argparse
import base64
import csv
import hashlib
import html
import math
import struct
import zlib
from pathlib import Path
from typing import Any, Iterable

from PIL import Image, ImageDraw, ImageFont


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCHEMA = "bvk2-si-convergence-morphology-v1"
CANONICAL_ENDPOINTS = (
    PROJECT_ROOT / "results/bvk2_boundary_crossing_panels/crossing_endpoints.csv"
)
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "results/bvk2_si_convergence_morphology"

# Each state is a canonical accepted endpoint, not a convenient standalone run.
# The shared chiN=54 point gives GYR, CYL, and LAM at the same thermodynamic
# coordinate; BCC comes from the all-target-stress S/C representative bracket.
SELECTIONS = (
    ("BCC", "S/C", 47.5, 0.2085, "phase_a"),
    ("CYL", "G/C", 54.0, 0.363, "phase_b"),
    ("GYR", "G/C", 54.0, 0.363, "phase_a"),
    ("LAM", "G/L", 54.0, 0.363, "phase_b"),
)

PHASE_COLORS = {
    "LAM": "#0072B2",
    "CYL": "#E69F00",
    "BCC": "#009E73",
    "GYR": "#CC79A7",
}


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def read_one(path: Path) -> dict[str, str]:
    rows = read_rows(path)
    if len(rows) != 1:
        raise ValueError(f"expected one row in {path}, found {len(rows)}")
    return rows[0]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(PROJECT_ROOT.resolve()).as_posix()


def resolve_artifact(value: str) -> Path:
    path = Path(value)
    if not path.is_absolute():
        path = PROJECT_ROOT / path
    path = path.resolve()
    if PROJECT_ROOT.resolve() not in (path, *path.parents):
        raise ValueError(f"artifact lies outside project: {path}")
    return path


def as_float(row: dict[str, str], key: str) -> float:
    value = float(row[key])
    if not math.isfinite(value):
        raise ValueError(f"{key} is not finite in {row}")
    return value


def parse_dims(value: str) -> tuple[int, ...]:
    return tuple(int(part) for part in value.split("x"))


def parse_lengths(summary: dict[str, str]) -> tuple[float, ...]:
    value = summary.get("reference_lengths") or summary.get("reference_length")
    if not value:
        raise ValueError(f"missing reference length in summary: {summary}")
    return tuple(float(part) for part in value.split(";"))


def product(values: Iterable[float]) -> float:
    answer = 1.0
    for value in values:
        answer *= value
    return answer


def _find_selection(
    rows: list[dict[str, str]],
    phase: str,
    transition: str,
    chi_n: float,
    coordinate: float,
    side: str,
) -> tuple[dict[str, str], str]:
    matches = [
        row
        for row in rows
        if row["transition"] == transition
        and math.isclose(float(row["chiN"]), chi_n, abs_tol=1.0e-12)
        and math.isclose(float(row["coordinate"]), coordinate, abs_tol=1.0e-12)
        and row[side] == phase
    ]
    if len(matches) != 1:
        raise ValueError(
            f"expected one canonical endpoint for {phase} at "
            f"{transition}, chiN={chi_n}, fA={coordinate}; found {len(matches)}"
        )
    return matches[0], side


def select_states() -> list[dict[str, Any]]:
    endpoint_rows = read_rows(CANONICAL_ENDPOINTS)
    states: list[dict[str, Any]] = []
    for phase, transition, chi_n, coordinate, side in SELECTIONS:
        endpoint, prefix = _find_selection(
            endpoint_rows, phase, transition, chi_n, coordinate, side
        )
        if endpoint["endpoint_status"] != "accepted":
            raise ValueError(f"canonical endpoint is not accepted: {endpoint}")
        if endpoint[f"{prefix}_status"] != "accepted":
            raise ValueError(f"selected phase endpoint is not accepted: {endpoint}")

        source = resolve_artifact(endpoint["source"])
        if sha256(source) != endpoint["source_sha256"]:
            raise ValueError(f"canonical endpoint source hash changed: {source}")

        artifact = resolve_artifact(endpoint[f"{prefix}_artifact"])
        summary_path = artifact / "summary.csv"
        trace_path = artifact / "convergence_trace.csv"
        density_path = artifact / "density.csv"
        comparison_path = artifact / "quadrature_comparison.csv"
        for required in (summary_path, trace_path, density_path):
            if not required.is_file():
                raise FileNotFoundError(required)

        summary = read_one(summary_path)
        if summary["phase"] != phase:
            raise ValueError(f"phase mismatch in {summary_path}")
        if not math.isclose(float(summary["chiN"]), chi_n, abs_tol=1.0e-12):
            raise ValueError(f"chiN mismatch in {summary_path}")
        if not math.isclose(float(summary["fA"]), coordinate, abs_tol=2.0e-10):
            raise ValueError(f"fA mismatch in {summary_path}")
        if summary["shell_pass"].lower() != "true":
            raise ValueError(f"morphology shell audit did not pass: {summary_path}")

        lengths = parse_lengths(summary)
        dims = parse_dims(summary["dims"])
        if len(lengths) != len(dims):
            raise ValueError(f"dimension/length mismatch in {summary_path}")
        volume = product(lengths)
        final_energy = as_float(summary, "final_energy")
        final_energy_density = final_energy / volume
        audit_energy_density = as_float(endpoint, f"{prefix}_energy_density")
        trace_rows = read_rows(trace_path)
        if not trace_rows:
            raise ValueError(f"empty convergence trace: {trace_path}")
        trace_last_energy = float(trace_rows[-1]["energy"])

        state = {
            "schema": SCHEMA,
            "phase": phase,
            "transition": transition,
            "chiN": chi_n,
            "fA": coordinate,
            "endpoint_role": endpoint["bracket_role"],
            "claim_kind": "phase_endpoint",
            "endpoint_status": endpoint["endpoint_status"],
            "phase_status": endpoint[f"{prefix}_status"],
            "grid": endpoint[f"{prefix}_grid"],
            "dims": summary["dims"],
            "reference_lengths": ";".join(f"{value:.15g}" for value in lengths),
            "energy_audit_factor": endpoint["energy_audit_factor"],
            "accepted_energy_density": audit_energy_density,
            "trace_objective_energy_density": final_energy_density,
            "audit_minus_trace_energy_density": (
                audit_energy_density - final_energy_density
            ),
            "cell_stress": as_float(endpoint, f"{prefix}_cell_stress"),
            "shell_label": endpoint[f"{prefix}_shell"],
            "expected_shell": summary["expected_shell"],
            "shell_pass": summary["shell_pass"].lower(),
            "mean_phi": as_float(summary, "mean_phi_fine"),
            "mean_error": as_float(summary, "mean_error"),
            "min_phi": as_float(summary, "min_phi"),
            "max_phi": as_float(summary, "max_phi"),
            "contrast": as_float(summary, "contrast"),
            "solver": summary["solver"],
            "iterations": int(summary["iterations"]),
            "solver_converged": summary["converged"].lower(),
            "stop_reason": summary["stop_reason"],
            "selected_plateau_energy_gap_density": (
                as_float(summary, "selected_plateau_energy_gap") / volume
            ),
            "trace_terminal_gap_density": abs(
                trace_last_energy - final_energy
            )
            / volume,
            "final_projected_theta_r2": as_float(
                summary, "final_projected_theta_r2"
            ),
            "final_projected_theta_rinf": as_float(
                summary, "final_projected_theta_rinf"
            ),
            "artifact": relative(artifact),
            "canonical_endpoint_source": relative(source),
            "canonical_endpoint_source_sha256": endpoint["source_sha256"],
            "summary_sha256": sha256(summary_path),
            "convergence_trace_sha256": sha256(trace_path),
            "density_sha256": sha256(density_path),
            "quadrature_comparison_sha256": (
                sha256(comparison_path) if comparison_path.is_file() else ""
            ),
            "selection_reason": (
                "accepted S/C endpoint with all endpoint stresses below 1e-3"
                if phase == "BCC"
                else "accepted shared chiN=54,fA=0.363 G-pocket endpoint"
            ),
            "_artifact_path": artifact,
            "_summary": summary,
            "_volume": volume,
        }
        states.append(state)
    return states


def build_convergence_rows(states: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for state in states:
        trace_path = state["_artifact_path"] / "convergence_trace.csv"
        selected_energy = float(state["_summary"]["final_energy"])
        volume = state["_volume"]
        for source_row in read_rows(trace_path):
            energy = float(source_row["energy"])
            rows.append(
                {
                    "schema": SCHEMA,
                    "phase": state["phase"],
                    "chiN": state["chiN"],
                    "fA": state["fA"],
                    "iteration": int(source_row["iteration"]),
                    "record_kind": "trace",
                    "energy": energy,
                    "energy_density": energy / volume,
                    "abs_gap_to_selected_energy_density": abs(
                        energy - selected_energy
                    )
                    / volume,
                    "projected_theta_r2": float(
                        source_row["projected_theta_r2"]
                    ),
                    "projected_theta_rinf": float(
                        source_row["projected_theta_rinf"]
                    ),
                    "residual_role": "diagnostic_only",
                }
            )
        rows.append(
            {
                "schema": SCHEMA,
                "phase": state["phase"],
                "chiN": state["chiN"],
                "fA": state["fA"],
                "iteration": int(state["_summary"]["selected_gradient_call"]),
                "record_kind": "selected_plateau",
                "energy": selected_energy,
                "energy_density": selected_energy / volume,
                "abs_gap_to_selected_energy_density": 0.0,
                "projected_theta_r2": state["final_projected_theta_r2"],
                "projected_theta_rinf": state["final_projected_theta_rinf"],
                "residual_role": "diagnostic_only",
            }
        )
    return rows


def _read_density_header(path: Path) -> list[str]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle)
        return next(reader)


def _max_variance_k(path: Path, dims: tuple[int, int, int]) -> int:
    nx, ny, nz = dims
    count = [0] * nz
    total = [0.0] * nz
    total2 = [0.0] * nz
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            k = int(row["k"]) - 1
            phi = float(row["phi"])
            count[k] += 1
            total[k] += phi
            total2[k] += phi * phi
    if any(value != nx * ny for value in count):
        raise ValueError(f"incomplete 3-D density in {path}")
    variances = [
        total2[k] / count[k] - (total[k] / count[k]) ** 2 for k in range(nz)
    ]
    return max(range(nz), key=lambda k: (variances[k], -k)) + 1


def build_morphology_rows(
    states: list[dict[str, Any]],
) -> tuple[list[dict[str, Any]], dict[str, dict[str, Any]]]:
    rows: list[dict[str, Any]] = []
    metadata: dict[str, dict[str, Any]] = {}
    for state in states:
        phase = state["phase"]
        path = state["_artifact_path"] / "density.csv"
        dims = parse_dims(state["dims"])
        header = _read_density_header(path)
        if len(dims) == 1 and header != ["i", "phi"]:
            raise ValueError(f"unexpected 1-D density schema: {path}")
        if len(dims) == 2 and header != ["i", "j", "phi"]:
            raise ValueError(f"unexpected 2-D density schema: {path}")
        if len(dims) == 3 and header != ["i", "j", "k", "phi"]:
            raise ValueError(f"unexpected 3-D density schema: {path}")

        plane_k = _max_variance_k(path, dims) if len(dims) == 3 else None
        nx = dims[0]
        ny = dims[1] if len(dims) >= 2 else 1
        selected = 0
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            for source in reader:
                if plane_k is not None and int(source["k"]) != plane_k:
                    continue
                i = int(source["i"])
                j = int(source.get("j", "1"))
                phi = float(source["phi"])
                rows.append(
                    {
                        "schema": SCHEMA,
                        "phase": phase,
                        "chiN": state["chiN"],
                        "fA": state["fA"],
                        "section_kind": (
                            "periodic_profile"
                            if len(dims) == 1
                            else "full_periodic_cell"
                            if len(dims) == 2
                            else "max_variance_xy_plane"
                        ),
                        "fixed_axis": "none" if plane_k is None else "k",
                        "fixed_index": "" if plane_k is None else plane_k,
                        "i": i,
                        "j": "" if len(dims) == 1 else j,
                        "u": (i - 0.5) / nx,
                        "v": "" if len(dims) == 1 else (j - 0.5) / ny,
                        "phi": phi,
                    }
                )
                selected += 1
        expected = nx if len(dims) == 1 else nx * ny
        if selected != expected:
            raise ValueError(
                f"selected {selected} density samples for {phase}, expected {expected}"
            )
        metadata[phase] = {
            "dims": dims,
            "section_kind": rows[-1]["section_kind"],
            "fixed_index": plane_k,
            "samples": selected,
        }
    return rows, metadata


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty table: {path}")
    public_rows = [
        {key: value for key, value in row.items() if not key.startswith("_")}
        for row in rows
    ]
    fieldnames = list(public_rows[0])
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(public_rows)


def _svg_text(
    x: float,
    y: float,
    value: str,
    *,
    size: int = 14,
    anchor: str = "start",
    weight: str = "normal",
    fill: str = "#222222",
) -> str:
    return (
        f'<text x="{x:.2f}" y="{y:.2f}" font-family="Arial,Helvetica,sans-serif" '
        f'font-size="{size}" text-anchor="{anchor}" font-weight="{weight}" '
        f'fill="{fill}">{html.escape(value)}</text>'
    )


def _interpolate_color(phi: float) -> tuple[int, int, int]:
    # Compact viridis-like, monotonically lightening color map.
    anchors = (
        (0.00, (68, 1, 84)),
        (0.25, (59, 82, 139)),
        (0.50, (33, 145, 140)),
        (0.75, (94, 201, 98)),
        (1.00, (253, 231, 37)),
    )
    value = min(1.0, max(0.0, phi))
    for (x0, c0), (x1, c1) in zip(anchors, anchors[1:]):
        if value <= x1:
            t = (value - x0) / (x1 - x0)
            return tuple(
                round(c0[index] + t * (c1[index] - c0[index]))
                for index in range(3)
            )
    return anchors[-1][1]


def _png_data_uri(width: int, height: int, values: list[float]) -> str:
    if len(values) != width * height:
        raise ValueError("PNG value count does not match dimensions")
    raw = bytearray()
    for row in range(height):
        raw.append(0)
        for column in range(width):
            raw.extend(_interpolate_color(values[row * width + column]))

    def chunk(kind: bytes, data: bytes) -> bytes:
        payload = kind + data
        return (
            struct.pack(">I", len(data))
            + payload
            + struct.pack(">I", zlib.crc32(payload) & 0xFFFFFFFF)
        )

    png = b"\x89PNG\r\n\x1a\n"
    png += chunk(b"IHDR", struct.pack(">IIBBBBB", width, height, 8, 2, 0, 0, 0))
    png += chunk(b"IDAT", zlib.compress(bytes(raw), level=9))
    png += chunk(b"IEND", b"")
    return "data:image/png;base64," + base64.b64encode(png).decode("ascii")


def _log_limits(values: list[float], floor: float) -> tuple[float, float]:
    positive = [max(value, floor) for value in values if math.isfinite(value)]
    lower = 10 ** math.floor(math.log10(min(positive)))
    upper = 10 ** math.ceil(math.log10(max(positive)))
    if upper <= lower:
        upper = lower * 10.0
    return lower, upper


def write_svg(
    path: Path,
    states: list[dict[str, Any]],
    convergence: list[dict[str, Any]],
    morphology: list[dict[str, Any]],
    metadata: dict[str, dict[str, Any]],
) -> None:
    width, height = 1500, 1030
    lines = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" '
        f'height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        _svg_text(
            750,
            36,
            "BVK2 accepted-state convergence and morphology diagnostics",
            size=22,
            anchor="middle",
            weight="bold",
        ),
        _svg_text(
            750,
            60,
            "Canonical accepted boundary endpoints; no field solve rerun",
            size=13,
            anchor="middle",
            fill="#555555",
        ),
    ]

    # Top convergence panels.
    panels = (
        (70, 100, 650, 350, "Energy plateau", "energy"),
        (780, 100, 650, 350, "Projected θ residual (diagnostic)", "residual"),
    )
    trace_only = [row for row in convergence if row["record_kind"] == "trace"]
    max_iteration = max(int(row["iteration"]) for row in trace_only)
    energy_floor = 1.0e-14
    energy_limits = _log_limits(
        [
            float(row["abs_gap_to_selected_energy_density"])
            for row in trace_only
        ],
        energy_floor,
    )
    residual_limits = _log_limits(
        [float(row["projected_theta_r2"]) for row in trace_only], 1.0e-8
    )

    for x0, y0, panel_w, panel_h, title, kind in panels:
        px0, px1 = x0 + 76, x0 + panel_w - 24
        py0, py1 = y0 + 44, y0 + panel_h - 55
        limits = energy_limits if kind == "energy" else residual_limits
        log0, log1 = math.log10(limits[0]), math.log10(limits[1])

        def sx(value: float) -> float:
            return px0 + value / max_iteration * (px1 - px0)

        def sy(value: float) -> float:
            safe = max(
                value, energy_floor if kind == "energy" else residual_limits[0]
            )
            return py1 - (math.log10(safe) - log0) / (log1 - log0) * (py1 - py0)

        lines.extend(
            [
                _svg_text(x0 + 10, y0 + 22, title, size=17, weight="bold"),
                f'<rect x="{px0}" y="{py0}" width="{px1-px0}" '
                f'height="{py1-py0}" fill="#fafafa" stroke="#444444"/>',
            ]
        )
        for exponent in range(math.floor(log0), math.ceil(log1) + 1):
            value = 10.0**exponent
            yy = sy(value)
            lines.append(
                f'<line x1="{px0}" y1="{yy:.2f}" x2="{px1}" y2="{yy:.2f}" '
                'stroke="#dddddd" stroke-width="1"/>'
            )
            lines.append(
                _svg_text(
                    px0 - 9,
                    yy + 4,
                    f"10^{exponent}",
                    size=11,
                    anchor="end",
                    fill="#555555",
                )
            )
        for tick in range(0, max_iteration + 1, 200):
            xx = sx(tick)
            lines.append(
                f'<line x1="{xx:.2f}" y1="{py1:.2f}" x2="{xx:.2f}" '
                f'y2="{py1+5:.2f}" stroke="#444444"/>'
            )
            lines.append(
                _svg_text(
                    xx, py1 + 20, str(tick), size=11, anchor="middle", fill="#555555"
                )
            )
        for phase in ("LAM", "CYL", "BCC", "GYR"):
            phase_rows = sorted(
                (row for row in trace_only if row["phase"] == phase),
                key=lambda row: int(row["iteration"]),
            )
            key = (
                "abs_gap_to_selected_energy_density"
                if kind == "energy"
                else "projected_theta_r2"
            )
            points = " ".join(
                f"{sx(int(row['iteration'])):.2f},{sy(float(row[key])):.2f}"
                for row in phase_rows
            )
            color = PHASE_COLORS[phase]
            lines.append(
                f'<polyline points="{points}" fill="none" stroke="{color}" '
                'stroke-width="2.2"/>'
            )
            last = phase_rows[-1]
            lines.append(
                f'<circle cx="{sx(int(last["iteration"])):.2f}" '
                f'cy="{sy(float(last[key])):.2f}" r="3.5" fill="{color}"/>'
            )
        lines.extend(
            [
                _svg_text(
                    (px0 + px1) / 2,
                    y0 + panel_h - 13,
                    "iteration",
                    size=13,
                    anchor="middle",
                ),
                _svg_text(
                    x0 + 16,
                    (py0 + py1) / 2,
                    (
                        "|e − e_selected|"
                        if kind == "energy"
                        else "projected θ R2"
                    ),
                    size=12,
                    anchor="middle",
                ),
            ]
        )

    legend_x = 1010
    for index, phase in enumerate(("LAM", "CYL", "BCC", "GYR")):
        xx = legend_x + index * 95
        lines.append(
            f'<line x1="{xx}" y1="84" x2="{xx+24}" y2="84" '
            f'stroke="{PHASE_COLORS[phase]}" stroke-width="3"/>'
        )
        lines.append(_svg_text(xx + 30, 89, phase, size=12))

    # Bottom morphology cards.
    morphology_by_phase: dict[str, list[dict[str, Any]]] = {
        phase: [row for row in morphology if row["phase"] == phase]
        for phase in ("LAM", "CYL", "BCC", "GYR")
    }
    state_by_phase = {state["phase"]: state for state in states}
    card_y, card_h, card_w, gap = 505, 440, 325, 25
    for index, phase in enumerate(("LAM", "CYL", "BCC", "GYR")):
        state = state_by_phase[phase]
        x0 = 50 + index * (card_w + gap)
        lines.append(
            f'<rect x="{x0}" y="{card_y}" width="{card_w}" height="{card_h}" '
            'rx="5" fill="#fbfbfb" stroke="#b0b0b0"/>'
        )
        lines.append(
            _svg_text(
                x0 + 16,
                card_y + 28,
                phase,
                size=18,
                weight="bold",
                fill=PHASE_COLORS[phase],
            )
        )
        lines.append(
            _svg_text(
                x0 + card_w - 15,
                card_y + 27,
                f"χN={state['chiN']:g}, f_A={state['fA']:.4f}",
                size=11,
                anchor="end",
                fill="#555555",
            )
        )
        phase_rows = morphology_by_phase[phase]
        plot_x, plot_y, plot_w, plot_h = x0 + 28, card_y + 54, 269, 268
        if phase == "LAM":
            lines.append(
                f'<rect x="{plot_x}" y="{plot_y}" width="{plot_w}" '
                f'height="{plot_h}" fill="white" stroke="#444444"/>'
            )
            points = " ".join(
                f"{plot_x + float(row['u']) * plot_w:.2f},"
                f"{plot_y + (1.0-float(row['phi'])) * plot_h:.2f}"
                for row in phase_rows
            )
            lines.append(
                f'<polyline points="{points}" fill="none" '
                f'stroke="{PHASE_COLORS[phase]}" stroke-width="2"/>'
            )
            lines.extend(
                [
                    _svg_text(
                        plot_x + plot_w / 2,
                        plot_y + plot_h + 18,
                        "x / D",
                        size=11,
                        anchor="middle",
                    ),
                    _svg_text(plot_x - 8, plot_y + 4, "1", size=10, anchor="end"),
                    _svg_text(
                        plot_x - 8, plot_y + plot_h + 4, "0", size=10, anchor="end"
                    ),
                ]
            )
        else:
            nx, ny = metadata[phase]["dims"][:2]
            image_values = [0.0] * (nx * ny)
            for row in phase_rows:
                i, j = int(row["i"]) - 1, int(row["j"]) - 1
                # SVG/PNG rows run from top to bottom; flip j for Cartesian view.
                image_values[(ny - 1 - j) * nx + i] = float(row["phi"])
            uri = _png_data_uri(nx, ny, image_values)
            lines.append(
                f'<image x="{plot_x}" y="{plot_y}" width="{plot_w}" '
                f'height="{plot_h}" preserveAspectRatio="none" href="{uri}"/>'
            )
            lines.append(
                f'<rect x="{plot_x}" y="{plot_y}" width="{plot_w}" '
                f'height="{plot_h}" fill="none" stroke="#444444"/>'
            )
        section = metadata[phase]["section_kind"]
        if metadata[phase]["fixed_index"] is not None:
            section += f" (k={metadata[phase]['fixed_index']})"
        lines.extend(
            [
                _svg_text(
                    x0 + card_w / 2,
                    card_y + 341,
                    section.replace("_", " "),
                    size=11,
                    anchor="middle",
                    fill="#555555",
                ),
                _svg_text(
                    x0 + 16,
                    card_y + 371,
                    f"grid {state['grid']}   shell {state['shell_label']}",
                    size=11,
                ),
                _svg_text(
                    x0 + 16,
                    card_y + 393,
                    f"stress {state['cell_stress']:.3g}   Δφ {state['contrast']:.5f}",
                    size=11,
                ),
                _svg_text(
                    x0 + 16,
                    card_y + 415,
                    (
                        f"plateau gap {state['selected_plateau_energy_gap_density']:.2e}"
                    ),
                    size=11,
                ),
            ]
        )

    # Shared color bar and caveat.
    bar_x, bar_y, bar_w, bar_h = 1240, 968, 180, 12
    gradient_id = "phi_gradient"
    lines.extend(
        [
            f'<defs><linearGradient id="{gradient_id}">',
            '<stop offset="0%" stop-color="#440154"/>',
            '<stop offset="25%" stop-color="#3b528b"/>',
            '<stop offset="50%" stop-color="#21918c"/>',
            '<stop offset="75%" stop-color="#5ec962"/>',
            '<stop offset="100%" stop-color="#fde725"/>',
            "</linearGradient></defs>",
            f'<rect x="{bar_x}" y="{bar_y}" width="{bar_w}" height="{bar_h}" '
            f'fill="url(#{gradient_id})" stroke="#444444"/>',
            _svg_text(bar_x - 8, bar_y + 11, "0", size=10, anchor="end"),
            _svg_text(bar_x + bar_w + 8, bar_y + 11, "1", size=10),
            _svg_text(
                bar_x + bar_w / 2,
                bar_y - 5,
                "φ_A",
                size=11,
                anchor="middle",
            ),
            _svg_text(
                60,
                981,
                (
                    "Projected θ norms are diagnostics, not universal acceptance "
                    "gates; phase identity is certified by the shell audit."
                ),
                size=12,
                fill="#555555",
            ),
            _svg_text(
                60,
                1004,
                (
                    "BCC/GYR show the deterministic maximum-variance xy section; "
                    "sections illustrate morphology but do not replace 3-D shell tests."
                ),
                size=12,
                fill="#555555",
            ),
            "</svg>",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _preview_font(size: int, *, bold: bool = False) -> ImageFont.FreeTypeFont:
    name = "DejaVuSans-Bold.ttf" if bold else "DejaVuSans.ttf"
    candidates = (
        Path("/usr/share/fonts/truetype/dejavu") / name,
        Path("/usr/share/fonts/dejavu") / name,
    )
    for path in candidates:
        if path.is_file():
            return ImageFont.truetype(str(path), size=size)
    return ImageFont.load_default()


def _draw_centered(
    draw: ImageDraw.ImageDraw,
    xy: tuple[float, float],
    value: str,
    *,
    font: ImageFont.ImageFont,
    fill: str = "#222222",
) -> None:
    box = draw.textbbox((0, 0), value, font=font)
    width = box[2] - box[0]
    draw.text((xy[0] - width / 2, xy[1]), value, font=font, fill=fill)


def write_png_preview(
    path: Path,
    states: list[dict[str, Any]],
    convergence: list[dict[str, Any]],
    morphology: list[dict[str, Any]],
    metadata: dict[str, dict[str, Any]],
) -> None:
    """Render a lightweight inspection preview from the same tabulated data."""
    canvas = Image.new("RGB", (1500, 1030), "white")
    draw = ImageDraw.Draw(canvas)
    title_font = _preview_font(22, bold=True)
    heading_font = _preview_font(17, bold=True)
    body_font = _preview_font(12)
    small_font = _preview_font(10)
    _draw_centered(
        draw,
        (750, 18),
        "BVK2 accepted-state convergence and morphology diagnostics",
        font=title_font,
    )
    _draw_centered(
        draw,
        (750, 48),
        "PNG inspection preview; SVG is authoritative",
        font=body_font,
        fill="#555555",
    )

    trace_only = [row for row in convergence if row["record_kind"] == "trace"]
    max_iteration = max(int(row["iteration"]) for row in trace_only)
    energy_floor = 1.0e-14
    energy_limits = _log_limits(
        [
            float(row["abs_gap_to_selected_energy_density"])
            for row in trace_only
        ],
        energy_floor,
    )
    residual_limits = _log_limits(
        [float(row["projected_theta_r2"]) for row in trace_only], 1.0e-8
    )
    for x0, y0, panel_w, panel_h, title, kind in (
        (70, 100, 650, 350, "Energy plateau", "energy"),
        (780, 100, 650, 350, "Projected theta R2 (diagnostic)", "residual"),
    ):
        px0, px1 = x0 + 76, x0 + panel_w - 24
        py0, py1 = y0 + 44, y0 + panel_h - 55
        limits = energy_limits if kind == "energy" else residual_limits
        log0, log1 = math.log10(limits[0]), math.log10(limits[1])
        draw.text((x0 + 10, y0 + 7), title, font=heading_font, fill="#222222")
        draw.rectangle((px0, py0, px1, py1), fill="#fafafa", outline="#444444")

        def sx(value: float) -> float:
            return px0 + value / max_iteration * (px1 - px0)

        def sy(value: float) -> float:
            safe = max(value, energy_floor if kind == "energy" else limits[0])
            return py1 - (math.log10(safe) - log0) / (log1 - log0) * (py1 - py0)

        for exponent in range(math.floor(log0), math.ceil(log1) + 1):
            value = 10.0**exponent
            yy = sy(value)
            draw.line((px0, yy, px1, yy), fill="#dddddd")
            label = f"10^{exponent}"
            box = draw.textbbox((0, 0), label, font=small_font)
            draw.text(
                (px0 - (box[2] - box[0]) - 8, yy - 6),
                label,
                font=small_font,
                fill="#555555",
            )
        for tick in range(0, max_iteration + 1, 200):
            xx = sx(tick)
            draw.line((xx, py1, xx, py1 + 5), fill="#444444")
            _draw_centered(
                draw, (xx, py1 + 8), str(tick), font=small_font, fill="#555555"
            )
        for phase in ("LAM", "CYL", "BCC", "GYR"):
            rows = sorted(
                (row for row in trace_only if row["phase"] == phase),
                key=lambda row: int(row["iteration"]),
            )
            key = (
                "abs_gap_to_selected_energy_density"
                if kind == "energy"
                else "projected_theta_r2"
            )
            points = [
                (sx(int(row["iteration"])), sy(float(row[key]))) for row in rows
            ]
            draw.line(points, fill=PHASE_COLORS[phase], width=3)
            xx, yy = points[-1]
            draw.ellipse((xx - 3, yy - 3, xx + 3, yy + 3), fill=PHASE_COLORS[phase])
        _draw_centered(
            draw, ((px0 + px1) / 2, y0 + panel_h - 26), "iteration", font=body_font
        )

    for index, phase in enumerate(("LAM", "CYL", "BCC", "GYR")):
        x = 1010 + index * 95
        draw.line((x, 84, x + 24, 84), fill=PHASE_COLORS[phase], width=3)
        draw.text((x + 30, 77), phase, font=body_font, fill="#222222")

    morphology_by_phase = {
        phase: [row for row in morphology if row["phase"] == phase]
        for phase in ("LAM", "CYL", "BCC", "GYR")
    }
    state_by_phase = {state["phase"]: state for state in states}
    card_y, card_h, card_w, gap = 505, 440, 325, 25
    for index, phase in enumerate(("LAM", "CYL", "BCC", "GYR")):
        state = state_by_phase[phase]
        x0 = 50 + index * (card_w + gap)
        draw.rounded_rectangle(
            (x0, card_y, x0 + card_w, card_y + card_h),
            radius=5,
            fill="#fbfbfb",
            outline="#b0b0b0",
        )
        draw.text(
            (x0 + 16, card_y + 9),
            phase,
            font=heading_font,
            fill=PHASE_COLORS[phase],
        )
        draw.text(
            (x0 + 158, card_y + 12),
            f"chiN={state['chiN']:g}, fA={state['fA']:.4f}",
            font=small_font,
            fill="#555555",
        )
        rows = morphology_by_phase[phase]
        plot_x, plot_y, plot_w, plot_h = x0 + 28, card_y + 54, 269, 268
        if phase == "LAM":
            draw.rectangle(
                (plot_x, plot_y, plot_x + plot_w, plot_y + plot_h),
                fill="white",
                outline="#444444",
            )
            points = [
                (
                    plot_x + float(row["u"]) * plot_w,
                    plot_y + (1.0 - float(row["phi"])) * plot_h,
                )
                for row in rows
            ]
            draw.line(points, fill=PHASE_COLORS[phase], width=2)
        else:
            nx, ny = metadata[phase]["dims"][:2]
            native = Image.new("RGB", (nx, ny))
            pixels = native.load()
            for row in rows:
                i, j = int(row["i"]) - 1, int(row["j"]) - 1
                pixels[i, ny - 1 - j] = _interpolate_color(float(row["phi"]))
            resampling = getattr(Image, "Resampling", Image)
            native = native.resize((plot_w, plot_h), resampling.BILINEAR)
            canvas.paste(native, (plot_x, plot_y))
            draw.rectangle(
                (plot_x, plot_y, plot_x + plot_w, plot_y + plot_h),
                outline="#444444",
            )
        section = metadata[phase]["section_kind"].replace("_", " ")
        if metadata[phase]["fixed_index"] is not None:
            section += f" (k={metadata[phase]['fixed_index']})"
        _draw_centered(
            draw,
            (x0 + card_w / 2, card_y + 329),
            section,
            font=small_font,
            fill="#555555",
        )
        draw.text(
            (x0 + 16, card_y + 366),
            f"grid {state['grid']}   shell {state['shell_label']}",
            font=small_font,
            fill="#222222",
        )
        draw.text(
            (x0 + 16, card_y + 388),
            f"stress {state['cell_stress']:.3g}   contrast {state['contrast']:.5f}",
            font=small_font,
            fill="#222222",
        )
        draw.text(
            (x0 + 16, card_y + 410),
            f"plateau gap {state['selected_plateau_energy_gap_density']:.2e}",
            font=small_font,
            fill="#222222",
        )
    draw.text(
        (60, 975),
        (
            "Projected theta norms are diagnostics, not universal acceptance "
            "gates; shell audits certify identity."
        ),
        font=body_font,
        fill="#555555",
    )
    draw.text(
        (60, 998),
        (
            "BCC/GYR use deterministic maximum-variance xy sections; "
            "the accepted reciprocal-shell audit remains the topology evidence."
        ),
        font=body_font,
        fill="#555555",
    )
    canvas.save(path, format="PNG", optimize=False, compress_level=9)


def write_readme(
    path: Path,
    states: list[dict[str, Any]],
    metadata: dict[str, dict[str, Any]],
) -> None:
    table = [
        "| Phase | Accepted endpoint | Grid | Cell stress | Shell | Stop reason |",
        "|---|---|---:|---:|---:|---|",
    ]
    for state in states:
        table.append(
            f"| {state['phase']} | {state['transition']} at "
            f"`chiN={state['chiN']:g}, fA={state['fA']:.4f}` | "
            f"{state['grid']} | {state['cell_stress']:.6g} | "
            f"{state['shell_label']} | {state['stop_reason']} |"
        )
    text = f"""# BVK2 SI convergence and morphology packet

This deterministic, solver-free packet curates one accepted representative of
LAM, CYL, BCC, and GYR from the canonical publication crossing endpoints:
`results/bvk2_boundary_crossing_panels/crossing_endpoints.csv`.

{chr(10).join(table)}

## Selection contract

- GYR, CYL, and LAM are the accepted endpoint states at the shared
  `chiN=54, fA=0.363` coordinate. The common-state comparison is already used
  to certify the nonempty GYR pocket.
- BCC is the accepted `chiN=47.5, fA=0.2085` S/C endpoint from the
  representative bracket where all endpoint stresses satisfy the `1e-3`
  target.
- Every selected row has `endpoint_status=accepted`,
  `phase_status=accepted`, and a passing morphology-shell audit. The script
  verifies the canonical endpoint source SHA-256 before reading the artifact.

## Files

- `representative_states.csv`: endpoint energy, stress, shell, plateau,
  residual, optimizer-stop, and content-hash provenance.
- `convergence_traces.csv`: raw objective trace normalized by cell volume plus
  the selected plateau record. Projected-theta norms are explicitly marked
  `diagnostic_only`.
- `morphology_sections.csv`: the complete 1-D LAM profile, complete 2-D CYL
  cell, and deterministic maximum-variance `xy` sections for BCC and GYR.
- `bvk2_si_convergence_morphology.svg`: publication-oriented convergence and
  morphology overview.
- `bvk2_si_convergence_morphology.png`: inspection preview rendered from the
  same tables. The SVG remains the authoritative publication asset.

## Interpretation and limitations

Acceptance is a boundary/phase-endpoint statement, not the value of the
optimizer's `converged` flag. Some exact-theta states stop at a bounded
iteration cap after an accepted energy plateau. The projected-theta R2/Rinf
values are therefore shown as diagnostics and are not silently converted into
a universal gate.

The canonical energy can use a higher post-solve quadrature factor than the
optimization trace. `accepted_energy_density` is the canonical endpoint value;
`trace_objective_energy_density` is the objective represented by the
convergence trace. Their difference is retained rather than hidden.

The BCC and GYR panels are scalar-field sections selected by maximum in-plane
variance (BCC `k={metadata['BCC']['fixed_index']}`, GYR
`k={metadata['GYR']['fixed_index']}`). They illustrate the accepted density
fields but do not by themselves prove 3-D topology; the reciprocal-shell
fingerprints (`110` and `211`) provide the phase-identity audit.

## Reproduction

```bash
python3 scripts/write_bvk2_si_convergence_morphology.py
python3 test/test_write_bvk2_si_convergence_morphology.py
```

No simulation is launched. Schema: `{SCHEMA}`.
"""
    path.write_text(text, encoding="utf-8")


def write_outputs(
    output_dir: Path,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    states = select_states()
    convergence = build_convergence_rows(states)
    morphology, metadata = build_morphology_rows(states)
    output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(output_dir / "representative_states.csv", states)
    write_csv(output_dir / "convergence_traces.csv", convergence)
    write_csv(output_dir / "morphology_sections.csv", morphology)
    write_svg(
        output_dir / "bvk2_si_convergence_morphology.svg",
        states,
        convergence,
        morphology,
        metadata,
    )
    write_png_preview(
        output_dir / "bvk2_si_convergence_morphology.png",
        states,
        convergence,
        morphology,
        metadata,
    )
    write_readme(output_dir / "README.md", states, metadata)
    return states, convergence, morphology


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=(
            "destination directory "
            "(default: results/bvk2_si_convergence_morphology)"
        ),
    )
    args = parser.parse_args()
    states, convergence, morphology = write_outputs(args.output_dir.resolve())
    print(
        f"wrote {len(states)} accepted states, {len(convergence)} convergence "
        f"records, and {len(morphology)} morphology samples to "
        f"{args.output_dir.resolve()}"
    )


if __name__ == "__main__":
    main()
