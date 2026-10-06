#!/usr/bin/env python3
"""Assemble publication-ready BVK2 free-energy crossing panels.

The script reads only accepted campaign artifacts and the canonical publication
ledger.  It does not rerun field minimizations or modify phase-boundary data.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import html
import math
from pathlib import Path
from typing import Any


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCHEMA = "bvk2-boundary-crossing-panels-v1"

PUBLICATION_LEDGER = (
    PROJECT_ROOT
    / "results/bvk2_publication_phase_diagram/accepted_phase_boundaries.csv"
)
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "results/bvk2_boundary_crossing_panels"


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def read_one(path: Path) -> dict[str, str]:
    rows = read_rows(path)
    if len(rows) != 1:
        raise ValueError(f"expected one row in {path}, found {len(rows)}")
    return rows[0]


def as_float(row: dict[str, str], key: str) -> float:
    value = float(row[key])
    if not math.isfinite(value):
        raise ValueError(f"{key} is not finite: {row[key]}")
    return value


def relative(path: Path) -> str:
    return path.resolve().relative_to(PROJECT_ROOT.resolve()).as_posix()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def secant_root(x0: float, y0: float, x1: float, y1: float) -> float:
    if y0 * y1 >= 0.0:
        raise ValueError("secant endpoints must have strict opposite signs")
    return x0 - y0 * (x1 - x0) / (y1 - y0)


def _accepted_publication_row(
    transition: str, chi_n: float, root_fa: float, boundary_source: Path
) -> dict[str, str]:
    matches = [
        row
        for row in read_rows(PUBLICATION_LEDGER)
        if row["transition"] == transition
        and math.isclose(float(row["chiN"]), chi_n, abs_tol=1.0e-10)
    ]
    if len(matches) != 1:
        raise ValueError(
            f"expected one publication row for {transition} at chiN={chi_n}"
        )
    row = matches[0]
    if row["status"] != "accepted" or row["claim_kind"] != "boundary_bracket":
        raise ValueError(f"publication row is not an accepted bracket: {row}")
    if not math.isclose(float(row["fA"]), root_fa, abs_tol=2.0e-12):
        raise ValueError(f"publication root disagrees with source: {row}")

    ledger_source = Path(row["source"])
    if not ledger_source.is_absolute():
        ledger_source = PROJECT_ROOT / ledger_source
    ledger_source = ledger_source.resolve()
    boundary = boundary_source.resolve()
    if ledger_source != boundary and ledger_source not in boundary.parents:
        raise ValueError(
            f"publication source {ledger_source} does not cover {boundary}"
        )
    return row


def _endpoint(
    *,
    transition: str,
    chi_n: float,
    coordinate: float,
    bracket_role: str,
    phase_a: str,
    phase_b: str,
    energy_a: float,
    energy_b: float,
    delta_f: float,
    grid_a: str,
    grid_b: str,
    stress_a: float,
    stress_b: float,
    shell_a: str,
    shell_b: str,
    status_a: str,
    status_b: str,
    audit_factor: str,
    artifact_a: str,
    artifact_b: str,
    source: Path,
) -> dict[str, Any]:
    if not math.isclose(energy_a - energy_b, delta_f, abs_tol=2.0e-12):
        raise ValueError(
            f"{transition} endpoint energy difference is inconsistent at {coordinate}"
        )
    if status_a != "accepted" or status_b not in {"accepted", "analytic_reference"}:
        raise ValueError(f"{transition} endpoint is not accepted at {coordinate}")
    return {
        "schema": SCHEMA,
        "transition": transition,
        "chiN": chi_n,
        "coordinate_name": "fA",
        "coordinate": coordinate,
        "bracket_role": bracket_role,
        "phase_a": phase_a,
        "phase_b": phase_b,
        "phase_a_energy_density": energy_a,
        "phase_b_energy_density": energy_b,
        "deltaF_phase_a_minus_phase_b": delta_f,
        "phase_a_grid": grid_a,
        "phase_b_grid": grid_b,
        "phase_a_cell_stress": stress_a,
        "phase_b_cell_stress": stress_b,
        "phase_a_shell": shell_a,
        "phase_b_shell": shell_b,
        "phase_a_status": status_a,
        "phase_b_status": status_b,
        "energy_audit_factor": audit_factor,
        "phase_a_artifact": artifact_a,
        "phase_b_artifact": artifact_b,
        "endpoint_status": "accepted",
        "source": relative(source),
        "source_sha256": sha256(source),
    }


def _summary(
    *,
    transition: str,
    chi_n: float,
    phase_a: str,
    phase_b: str,
    lower: float,
    upper: float,
    delta_lower: float,
    delta_upper: float,
    root_fa: float,
    root_uncertainty: float,
    root_method: str,
    crossing_kind: str,
    strict_signed_endpoint_root: bool,
    boundary_source: Path,
    endpoint_sources: list[Path],
) -> dict[str, Any]:
    width = upper - lower
    if width <= 0.0 or width > 0.0010000001:
        raise ValueError(f"{transition} bracket width is invalid: {width}")
    if delta_lower * delta_upper >= 0.0:
        raise ValueError(f"{transition} endpoints do not have strict opposite signs")
    diagnostic_root = secant_root(lower, delta_lower, upper, delta_upper)
    if not lower <= root_fa <= upper:
        raise ValueError(f"{transition} root is outside its bracket")
    if not math.isclose(root_uncertainty, width / 2.0, abs_tol=2.0e-12):
        raise ValueError(f"{transition} uncertainty is not half the bracket width")

    publication = _accepted_publication_row(
        transition, chi_n, root_fa, boundary_source
    )
    return {
        "schema": SCHEMA,
        "transition": transition,
        "chiN": chi_n,
        "phase_a": phase_a,
        "phase_b": phase_b,
        "deltaF_definition": f"F_{phase_a}-F_{phase_b}",
        "root_fA": root_fa,
        "root_uncertainty": root_uncertainty,
        "bracket_lower": lower,
        "bracket_upper": upper,
        "bracket_width": width,
        "deltaF_lower": delta_lower,
        "deltaF_upper": delta_upper,
        "diagnostic_secant_fA": diagnostic_root,
        "secant_slope": (delta_upper - delta_lower) / width,
        "root_method": root_method,
        "crossing_kind": crossing_kind,
        "strict_signed_endpoint_root": str(strict_signed_endpoint_root).lower(),
        "publication_status": publication["status"],
        "publication_source": publication["source"],
        "boundary_source": relative(boundary_source),
        "boundary_source_sha256": sha256(boundary_source),
        "endpoint_sources": ";".join(relative(path) for path in endpoint_sources),
        "endpoint_source_sha256": ";".join(sha256(path) for path in endpoint_sources),
    }


def _load_sdis() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    boundary = (
        PROJECT_ROOT
        / "results/bvk2_sdis_boundary_campaign/points/"
        "S_DIS_target_chiN14p000000.csv"
    )
    lower_source = (
        PROJECT_ROOT
        / "results/bvk2_sdis_boundary_campaign/evaluations/S_DIS/"
        "BCC_f0p301782_chiN14p000000.csv"
    )
    upper_source = (
        PROJECT_ROOT
        / "results/bvk2_sdis_boundary_campaign/evaluations/S_DIS/"
        "BCC_f0p302407_chiN14p000000.csv"
    )
    point = read_one(boundary)
    lower_row = read_one(lower_source)
    upper_row = read_one(upper_source)
    lower = as_float(point, "bracket_lower")
    upper = as_float(point, "bracket_upper")
    root = as_float(point, "root_fA")
    delta_lower = as_float(lower_row, "energy_density")
    delta_upper = as_float(upper_row, "energy_density")
    if not math.isclose(root, (lower + upper) / 2.0, abs_tol=2.0e-12):
        raise ValueError("S/DIS canonical root is not the survival-bracket midpoint")

    endpoints = []
    for role, row, source in (
        ("DIS side", lower_row, lower_source),
        ("BCC side", upper_row, upper_source),
    ):
        if row["accepted"].lower() != "true" or row["shell_pass"].lower() != "true":
            raise ValueError(f"S/DIS BCC endpoint failed its numerical gates: {row}")
        endpoints.append(
            _endpoint(
                transition="S/DIS",
                chi_n=14.0,
                coordinate=as_float(row, "fA"),
                bracket_role=role,
                phase_a="BCC",
                phase_b="DIS",
                energy_a=as_float(row, "energy_density"),
                energy_b=0.0,
                delta_f=as_float(row, "energy_density"),
                grid_a="48x48x48",
                grid_b="analytic homogeneous",
                stress_a=as_float(row, "cell_gradient_density"),
                stress_b=0.0,
                shell_a=row["shell_label"],
                shell_b="0",
                status_a="accepted",
                status_b="analytic_reference",
                audit_factor="1",
                artifact_a=row["density_path"],
                artifact_b="analytic DIS reference",
                source=source,
            )
        )
    summary = _summary(
        transition="S/DIS",
        chi_n=14.0,
        phase_a="BCC",
        phase_b="DIS",
        lower=lower,
        upper=upper,
        delta_lower=delta_lower,
        delta_upper=delta_upper,
        root_fa=root,
        root_uncertainty=(upper - lower) / 2.0,
        root_method="survival_bracket_midpoint",
        crossing_kind="ordered_branch_survival_edge",
        strict_signed_endpoint_root=False,
        boundary_source=boundary,
        endpoint_sources=[lower_source, upper_source],
    )
    return summary, endpoints


def _load_sc() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    directory = PROJECT_ROOT / "results/bvk2_ud_theta_sc_highchi/chi47p5"
    boundary = directory / "boundary_point.csv"
    endpoint_source = directory / "endpoint_energies.csv"
    point = read_one(boundary)
    rows = read_rows(endpoint_source)
    if len(rows) != 2:
        raise ValueError("S/C requires exactly two endpoint rows")
    endpoints = []
    for row in rows:
        endpoints.append(
            _endpoint(
                transition="S/C",
                chi_n=47.5,
                coordinate=as_float(row, "fA"),
                bracket_role=row["side"],
                phase_a="BCC",
                phase_b="CYL",
                energy_a=as_float(row, "bcc_energy_density_4x"),
                energy_b=as_float(row, "cyl_energy_density_4x"),
                delta_f=as_float(row, "deltaF_bcc_minus_cyl_4x"),
                grid_a="48x48x48",
                grid_b="48x84",
                stress_a=as_float(row, "bcc_cell_stress"),
                stress_b=as_float(row, "cyl_cell_stress"),
                shell_a=row["bcc_shell"],
                shell_b=row["cyl_shell"],
                status_a=row["status"],
                status_b=row["status"],
                audit_factor=row["audit_factor"],
                artifact_a=row["bcc_artifact"],
                artifact_b=row["cyl_artifact"],
                source=endpoint_source,
            )
        )
    summary = _summary(
        transition="S/C",
        chi_n=47.5,
        phase_a="BCC",
        phase_b="CYL",
        lower=as_float(point, "fA_lower"),
        upper=as_float(point, "fA_upper"),
        delta_lower=as_float(point, "deltaF_lower_4x"),
        delta_upper=as_float(point, "deltaF_upper_4x"),
        root_fa=as_float(point, "fA_secant"),
        root_uncertainty=as_float(point, "coordinate_half_width"),
        root_method="signed_endpoint_secant",
        crossing_kind="ordinary_free_energy_root",
        strict_signed_endpoint_root=True,
        boundary_source=boundary,
        endpoint_sources=[endpoint_source],
    )
    return summary, endpoints


def _load_gc() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    directory = PROJECT_ROOT / "results/bvk2_ud_theta_gc_chi54_direct"
    boundary = directory / "boundary_point.csv"
    endpoint_source = directory / "endpoint_energies.csv"
    point = read_one(boundary)
    rows = read_rows(endpoint_source)
    if len(rows) != 2:
        raise ValueError("G/C requires exactly two endpoint rows")
    endpoints = []
    for row in rows:
        endpoints.append(
            _endpoint(
                transition="G/C",
                chi_n=54.0,
                coordinate=as_float(row, "fA"),
                bracket_role=row["bracket_role"],
                phase_a="GYR",
                phase_b="CYL",
                energy_a=as_float(row, "gyr_energy_density_3x"),
                energy_b=as_float(row, "cyl_energy_density_3x"),
                delta_f=as_float(row, "deltaF_gyr_minus_cyl_3x"),
                grid_a=row["gyr_grid"],
                grid_b=row["cyl_grid"],
                stress_a=as_float(row, "gyr_stress"),
                stress_b=as_float(row, "cyl_stress"),
                shell_a=row["gyr_shell"],
                shell_b=row["cyl_shell"],
                status_a=row["gyr_status"],
                status_b=row["cyl_status"],
                audit_factor=row["audit_factor"],
                artifact_a=row["gyr_artifact"],
                artifact_b=row["cyl_artifact"],
                source=endpoint_source,
            )
        )
    summary = _summary(
        transition="G/C",
        chi_n=54.0,
        phase_a="GYR",
        phase_b="CYL",
        lower=as_float(point, "bracket_low"),
        upper=as_float(point, "bracket_high"),
        delta_lower=as_float(point, "delta_low"),
        delta_upper=as_float(point, "delta_high"),
        root_fa=as_float(point, "fA"),
        root_uncertainty=as_float(point, "fA_uncertainty"),
        root_method="signed_endpoint_secant",
        crossing_kind="ordinary_free_energy_root",
        strict_signed_endpoint_root=True,
        boundary_source=boundary,
        endpoint_sources=[endpoint_source],
    )
    return summary, endpoints


def _load_gl() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    directory = PROJECT_ROOT / "results/bvk2_ud_theta_gl_direct_fine/chi54"
    boundary = directory / "boundary_point.csv"
    endpoint_source = directory / "endpoint_energies.csv"
    point = read_one(boundary)
    rows = read_rows(endpoint_source)
    if len(rows) != 2:
        raise ValueError("G/L requires exactly two endpoint rows")
    endpoints = []
    for row in rows:
        endpoints.append(
            _endpoint(
                transition="G/L",
                chi_n=54.0,
                coordinate=as_float(row, "fA"),
                bracket_role=row["bracket_role"],
                phase_a="GYR",
                phase_b="LAM",
                energy_a=as_float(row, "gyr_energy_density_3x"),
                energy_b=as_float(row, "lam_energy_density_3x"),
                delta_f=as_float(row, "deltaF_gyr_minus_lam"),
                grid_a=row["gyr_grid"],
                grid_b=row["lam_grid"],
                stress_a=as_float(row, "gyr_stress"),
                stress_b=as_float(row, "lam_stress"),
                shell_a=row["gyr_shell"],
                shell_b=row["lam_shell"],
                status_a=row["gyr_status"],
                status_b=row["lam_status"],
                audit_factor=row["audit_factor"],
                artifact_a=row["gyr_artifact"],
                artifact_b=row["lam_artifact"],
                source=endpoint_source,
            )
        )
    summary = _summary(
        transition="G/L",
        chi_n=54.0,
        phase_a="GYR",
        phase_b="LAM",
        lower=as_float(point, "bracket_low"),
        upper=as_float(point, "bracket_high"),
        delta_lower=as_float(point, "delta_low"),
        delta_upper=as_float(point, "delta_high"),
        root_fa=as_float(point, "fA"),
        root_uncertainty=as_float(point, "fA_uncertainty"),
        root_method=point["root_method"],
        crossing_kind="ordinary_free_energy_root",
        strict_signed_endpoint_root=True,
        boundary_source=boundary,
        endpoint_sources=[endpoint_source],
    )
    return summary, endpoints


def build_dataset() -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    summaries: list[dict[str, Any]] = []
    endpoints: list[dict[str, Any]] = []
    for loader in (_load_sdis, _load_sc, _load_gc, _load_gl):
        summary, rows = loader()
        summaries.append(summary)
        endpoints.extend(rows)

    shared_gc = next(
        row
        for row in endpoints
        if row["transition"] == "G/C" and math.isclose(row["coordinate"], 0.363)
    )
    shared_gl = next(
        row
        for row in endpoints
        if row["transition"] == "G/L" and math.isclose(row["coordinate"], 0.363)
    )
    if not math.isclose(
        shared_gc["phase_a_energy_density"],
        shared_gl["phase_a_energy_density"],
        abs_tol=1.0e-12,
    ):
        raise ValueError("G/C and G/L do not share the same GYR state at fA=0.363")
    if not (
        shared_gc["deltaF_phase_a_minus_phase_b"] < 0.0
        and shared_gl["deltaF_phase_a_minus_phase_b"] < 0.0
    ):
        raise ValueError("the shared chiN=54 GYR state is not below both competitors")
    return summaries, endpoints


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty CSV: {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def _svg_text(
    x: float,
    y: float,
    text: str,
    *,
    size: int = 14,
    anchor: str = "start",
    weight: str = "normal",
    fill: str = "#222222",
) -> str:
    return (
        f'<text x="{x:.2f}" y="{y:.2f}" font-family="DejaVu Sans,Arial,sans-serif" '
        f'font-size="{size}" text-anchor="{anchor}" font-weight="{weight}" '
        f'fill="{fill}">{html.escape(text)}</text>'
    )


def write_svg(
    path: Path,
    summaries: list[dict[str, Any]],
    endpoints: list[dict[str, Any]],
) -> None:
    colors = {
        "S/DIS": "#0072b2",
        "S/C": "#d55e00",
        "G/C": "#7b4ab5",
        "G/L": "#238b45",
    }
    positions = [(70, 95), (665, 95), (70, 465), (665, 465)]
    panel_w, panel_h = 535, 285
    plot_left, plot_right, plot_top, plot_bottom = 78, 22, 52, 60
    lines = [
        '<svg xmlns="http://www.w3.org/2000/svg" width="1280" height="850" '
        'viewBox="0 0 1280 850">',
        '<rect width="1280" height="850" fill="white"/>',
        _svg_text(
            640,
            38,
            "BVK2 representative phase-boundary free-energy crossings",
            size=22,
            anchor="middle",
            weight="bold",
        ),
        _svg_text(
            640,
            63,
            "Accepted fixed-χN brackets; shaded width gives conservative coordinate uncertainty",
            size=13,
            anchor="middle",
            fill="#555555",
        ),
    ]

    for summary, (x0, y0) in zip(summaries, positions):
        transition = summary["transition"]
        rows = [row for row in endpoints if row["transition"] == transition]
        rows.sort(key=lambda row: row["coordinate"])
        color = colors[transition]
        lower = summary["bracket_lower"]
        upper = summary["bracket_upper"]
        width = upper - lower
        xmin, xmax = lower - 0.24 * width, upper + 0.24 * width
        values = [row["deltaF_phase_a_minus_phase_b"] * 1000.0 for row in rows]
        ymin, ymax = min(values + [0.0]), max(values + [0.0])
        span = max(ymax - ymin, 0.2)
        ymin -= 0.24 * span
        ymax += 0.24 * span

        px0 = x0 + plot_left
        px1 = x0 + panel_w - plot_right
        py0 = y0 + plot_top
        py1 = y0 + panel_h - plot_bottom

        def sx(value: float) -> float:
            return px0 + (value - xmin) / (xmax - xmin) * (px1 - px0)

        def sy(value: float) -> float:
            return py1 - (value - ymin) / (ymax - ymin) * (py1 - py0)

        lines.extend(
            [
                f'<rect x="{x0}" y="{y0}" width="{panel_w}" height="{panel_h}" '
                'rx="4" fill="#ffffff" stroke="#c7ccd1"/>',
                _svg_text(
                    x0 + 18,
                    y0 + 28,
                    f"{transition} at χN={summary['chiN']:g}",
                    size=17,
                    weight="bold",
                    fill=color,
                ),
                _svg_text(
                    x0 + panel_w - 18,
                    y0 + 27,
                    (
                        "accepted survival bracket"
                        if transition == "S/DIS"
                        else "accepted signed bracket"
                    ),
                    size=11,
                    anchor="end",
                    fill="#555555",
                ),
                f'<rect x="{sx(lower):.2f}" y="{py0:.2f}" '
                f'width="{sx(upper) - sx(lower):.2f}" height="{py1 - py0:.2f}" '
                f'fill="{color}" fill-opacity="0.08"/>',
                f'<line x1="{px0:.2f}" y1="{sy(0):.2f}" x2="{px1:.2f}" '
                f'y2="{sy(0):.2f}" stroke="#444444" stroke-width="1.2"/>',
                f'<line x1="{px0:.2f}" y1="{py0:.2f}" x2="{px0:.2f}" '
                f'y2="{py1:.2f}" stroke="#222222"/>',
                f'<line x1="{px0:.2f}" y1="{py1:.2f}" x2="{px1:.2f}" '
                f'y2="{py1:.2f}" stroke="#222222"/>',
            ]
        )

        root_x = sx(summary["root_fA"])
        lines.append(
            f'<line x1="{root_x:.2f}" y1="{py0:.2f}" x2="{root_x:.2f}" '
            f'y2="{py1:.2f}" stroke="{color}" stroke-width="1.2" '
            'stroke-dasharray="3,4"/>'
        )
        dash = ' stroke-dasharray="7,5"' if transition == "S/DIS" else ""
        lines.append(
            f'<line x1="{sx(rows[0]["coordinate"]):.2f}" '
            f'y1="{sy(values[0]):.2f}" x2="{sx(rows[1]["coordinate"]):.2f}" '
            f'y2="{sy(values[1]):.2f}" stroke="{color}" stroke-width="2.4"{dash}/>'
        )
        for row, value in zip(rows, values):
            lines.append(
                f'<circle cx="{sx(row["coordinate"]):.2f}" cy="{sy(value):.2f}" '
                f'r="5.5" fill="white" stroke="{color}" stroke-width="2.5"/>'
            )
            label_y = sy(value) - 10 if value >= 0 else sy(value) + 19
            lines.append(
                _svg_text(
                    sx(row["coordinate"]),
                    label_y,
                    f"{value:+.3f}",
                    size=11,
                    anchor="middle",
                    fill=color,
                )
            )

        zero_y = sy(0.0)
        diamond = (
            f"{root_x:.2f},{zero_y - 6:.2f} "
            f"{root_x + 6:.2f},{zero_y:.2f} "
            f"{root_x:.2f},{zero_y + 6:.2f} "
            f"{root_x - 6:.2f},{zero_y:.2f}"
        )
        lines.append(
            f'<polygon points="{diamond}" fill="{color}" stroke="white" '
            'stroke-width="1"/>'
        )
        lines.extend(
            [
                _svg_text(
                    x0 + 13,
                    (py0 + py1) / 2,
                    "ΔF × 10⁻³",
                    size=12,
                    anchor="middle",
                    fill="#333333",
                ),
                _svg_text(
                    (px0 + px1) / 2,
                    y0 + panel_h - 15,
                    "f_A",
                    size=13,
                    anchor="middle",
                ),
                _svg_text(
                    sx(lower),
                    py1 + 18,
                    f"{lower:.6f}",
                    size=10,
                    anchor="middle",
                    fill="#555555",
                ),
                _svg_text(
                    sx(upper),
                    py1 + 18,
                    f"{upper:.6f}",
                    size=10,
                    anchor="middle",
                    fill="#555555",
                ),
                _svg_text(
                    root_x,
                    py0 + 13,
                    f"root {summary['root_fA']:.6f}",
                    size=10,
                    anchor="middle",
                    fill=color,
                ),
                _svg_text(
                    x0 + panel_w - 18,
                    y0 + panel_h - 16,
                    f"ΔF = F_{summary['phase_a']} − F_{summary['phase_b']}",
                    size=11,
                    anchor="end",
                    fill="#555555",
                ),
            ]
        )

    lines.extend(
        [
            _svg_text(
                640,
                823,
                (
                    "S/DIS uses the ordered-branch survival midpoint; its secant is "
                    "diagnostic. Other panels use strict signed endpoint secants."
                ),
                size=12,
                anchor="middle",
                fill="#444444",
            ),
            "</svg>",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_readme(
    path: Path,
    summaries: list[dict[str, Any]],
    endpoints: list[dict[str, Any]],
) -> None:
    table = [
        "| Boundary | χN | Bracket | Root | Conservative uncertainty | Method |",
        "|---|---:|---:|---:|---:|---|",
    ]
    for row in summaries:
        table.append(
            f"| {row['transition']} | {row['chiN']:g} | "
            f"[{row['bracket_lower']:.9f}, {row['bracket_upper']:.9f}] | "
            f"{row['root_fA']:.9f} | ±{row['root_uncertainty']:.7f} | "
            f"{row['root_method']} |"
        )

    gc_shared = next(
        row
        for row in endpoints
        if row["transition"] == "G/C" and math.isclose(row["coordinate"], 0.363)
    )
    gl_shared = next(
        row
        for row in endpoints
        if row["transition"] == "G/L" and math.isclose(row["coordinate"], 0.363)
    )
    text = f"""# BVK2 representative boundary free-energy crossings

This directory is a deterministic publication packet assembled from accepted
fixed-`chiN` campaign endpoints. No field solve is rerun here.

{chr(10).join(table)}

`crossing_endpoints.csv` records the two accepted phase energies, their
difference, stress, morphology shell, grids, audit factor, source artifact, and
SHA-256 provenance for every panel. `crossing_summary.csv` records the bracket,
root, uncertainty, secant diagnostic, and canonical publication-ledger source.

## Interpretation

- S/C, G/C, and G/L are ordinary strict signed free-energy crossings. Their
  plotted roots are endpoint secants and the reported uncertainty is half the
  accepted bracket width.
- S/DIS is an order-disorder survival/bifurcation edge. Its canonical value is
  the accepted bracket midpoint, not the nearby secant diagnostic. The dashed
  connector in the figure makes this different semantics explicit.
- At `chiN=54`, `fA=0.363`, the shared GYR energy is
  `{gc_shared['phase_a_energy_density']:.15g}`. It is below CYL by
  `{abs(gc_shared['deltaF_phase_a_minus_phase_b']):.6g}` and below LAM by
  `{abs(gl_shared['deltaF_phase_a_minus_phase_b']):.6g}` in energy-density
  units. This common-state comparison directly establishes a nonempty GYR
  pocket without over-interpreting the separation between two interpolated
  roots.

The selected S/C slice is `chiN=47.5`: all four endpoint cell stresses are
inside the `1e-3` target and the independent factor-4 energy evaluation gives
strict opposite signs. The G/C and G/L panels use direct GYR `112^3`, CYL
`96x168`, and LAM `1024` states with factor-3 energy evaluation.

## Reproduction

From the `MonteCarlo` directory:

```bash
python3 scripts/write_bvk2_boundary_crossing_panels.py
python3 test/test_write_bvk2_boundary_crossing_panels.py
```

Generated files:

- `bvk2_boundary_crossings.svg`
- `crossing_endpoints.csv`
- `crossing_summary.csv`

Schema: `{SCHEMA}`.
"""
    path.write_text(text, encoding="utf-8")


def write_outputs(output_dir: Path) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    summaries, endpoints = build_dataset()
    output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(output_dir / "crossing_summary.csv", summaries)
    write_csv(output_dir / "crossing_endpoints.csv", endpoints)
    write_svg(output_dir / "bvk2_boundary_crossings.svg", summaries, endpoints)
    write_readme(output_dir / "README.md", summaries, endpoints)
    return summaries, endpoints


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="destination directory (default: results/bvk2_boundary_crossing_panels)",
    )
    args = parser.parse_args()
    summaries, endpoints = write_outputs(args.output_dir.resolve())
    print(
        f"wrote {len(summaries)} boundary panels and {len(endpoints)} endpoint rows "
        f"to {args.output_dir.resolve()}"
    )


if __name__ == "__main__":
    main()
