#!/usr/bin/env python3
"""Assemble the accepted BURP-TI boundaries into a publication phase diagram."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from math import cos, hypot, pi, sin
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent
RESULTS = PROJECT / "results"
DEFAULT_BASE = RESULTS / "burp_ti_gc_gl_kkt_boundary_campaign"
DEFAULT_SC = RESULTS / "burp_ti_sc_kkt_boundary_campaign" / "accepted_sc_boundaries.csv"
DEFAULT_SDIS = RESULTS / "burp_bvk1_sdis_boundary_trace" / "sdis_trace_roots.csv"
DEFAULT_SDIS_REFINEMENT = RESULTS / "burp_ti_sdis_refined_boundary" / "sdis_refined_roots.csv"
DEFAULT_SDIS_FIXED_FA = RESULTS / "burp_ti_sdis_fixed_fa_extension" / "sdis_fixed_fa_roots.csv"
DEFAULT_ODT = RESULTS / "burp_bvk1_ud_continuous_phase_diagram" / "burp_ti_oneperiod_boundaries.csv"
# SCFT overlay: Matsen revised diblock diagram (Fig. 2 + Table 1).
# Supersedes the Uneyama-Doi digitization, whose high-chiN gyroid region
# is inaccurate (it closes the Q230 pocket near chiN~38; Table 1 has Q230
# stable through chiN=100). See matsen_scft_phase_diagram_digitization/.
DEFAULT_SCFT_REFERENCE = (
    RESULTS / "matsen_scft_phase_diagram_digitization"
    / "matsen_scft_boundaries_overlay.csv"
)
DEFAULT_SCFT_SOURCE_SVG = (
    RESULTS / "matsen_scft_phase_diagram_digitization" / "extract_lowchi_from_zoom.py"
)
DEFAULT_OUTDIR = RESULTS / "burp_ti_publication_phase_diagram"
EXPECTED_SC_CHI = (60, 55, 50, 47.5, 45, 40, 35, 30, 25, 22, 20, 18, 16, 15)
EXPECTED_SDIS_CHI = (60, 55, 50, 45, 40, 35, 30, 25, 20, 15)
COLORS = {
    "S/DIS": "#1f77b4",
    "S/C": "#d62728",
    "G/C": "#9467bd",
    "G/L": "#ff7f0e",
    "C/L": "#2ca02c",
    "BCC/FCC": "#8c564b",
    "FCC/BCC": "#8c564b",
    "FCC/DIS": "#e377c2",
    "C/O70": "#17becf",
    "G/O70": "#bcbd22",
    "O70/L": "#393b79",
}
TRANSITIONS = ("S/C", "S/DIS", "C/L", "G/C", "G/L")
MARKERS = {
    "S/C": "circle",
    "S/DIS": "square",
    "C/L": "triangle",
    "G/C": "diamond",
    "G/L": "plus",
    "BCC/FCC": "diamond",
    "FCC/BCC": "diamond",
    "FCC/DIS": "square",
    "C/O70": "triangle",
    "G/O70": "diamond",
    "O70/L": "circle",
}
CRITICAL_POINT = (0.5, 10.508416112144529)
ODT_COLOR = "#17becf"
SCFT_COLOR = "#8b8f97"
RPA_COLOR = SCFT_COLOR
ESTIMATED_POINT_COLOR = "#2f2f2f"
EstimatedTriplePoint = tuple[float, float, tuple[str, ...]]


def truth(value: object) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes"}


def finite_float(value: object) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if number == number and abs(number) != float("inf") else None


def read_rows(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def interpolate(rows: list[dict[str, object]], chi: float) -> float | None:
    points = sorted(
        (float(row["chiN"]), float(row["fA"]))
        for row in rows
        if finite_float(row.get("chiN")) is not None
        and finite_float(row.get("fA")) is not None
    )
    if not points or chi < points[0][0] or chi > points[-1][0]:
        return None
    for x, y in points:
        if abs(x - chi) < 1.0e-9:
            return y
    for (x0, y0), (x1, y1) in zip(points, points[1:]):
        if x0 <= chi <= x1:
            weight = (chi - x0) / (x1 - x0)
            return y0 + weight * (y1 - y0)
    return None


def accepted_base(path: Path) -> list[dict[str, object]]:
    accepted: list[dict[str, object]] = []
    for row in read_rows(path):
        accepted.append(
            {
                "transition": row["transition"],
                "chiN": float(row["chiN"]),
                "fA": float(row["fA"]),
                "bracket_width": row["bracket_width"],
                "status": "accepted",
                "acceptance_basis": row["provenance"],
                "source": row["point_path"],
            }
        )
    return accepted


def cylinder_envelope(base: list[dict[str, object]], chi: float) -> float | None:
    closure_chi = min((float(row["chiN"]) for row in base if row["transition"] == "C/L"), default=35.0)
    transition = "G/C" if chi < closure_chi else "C/L"
    return interpolate([row for row in base if row["transition"] == transition], chi)


def accepted_sc(path: Path, base: list[dict[str, object]]) -> list[dict[str, object]]:
    accepted: list[dict[str, object]] = []
    for row in read_rows(path):
        kkt_schema = "root_fA" in row
        chi = finite_float(row.get("root_chiN" if kkt_schema else "chiN"))
        root = finite_float(row.get("root_fA" if kkt_schema else "fA_root"))
        envelope = cylinder_envelope(base, chi) if chi is not None else None
        if kkt_schema:
            pairwise_pass = row.get("status") == "accepted" and truth(row.get("phases_accepted"))
            dis_pass = True
            basis = "bounded-KKT BCC/CYL root; stress-free cells; accepted cylinder-side envelope"
        else:
            pairwise_pass = all(
                truth(row.get(field))
                for field in ("bracketed", "resolution_reached", "gates_pass", "class_ok")
            ) and row.get("blocker", "none") == "none"
            dis_pass = (
                finite_float(row.get("e_bcc_root")) is not None
                and finite_float(row.get("e_cyl_root")) is not None
                and float(row["e_bcc_root"]) < 0.0
                and float(row["e_cyl_root"]) < 0.0
            )
            basis = "BCC/CYL root; DIS energies negative; accepted cylinder-side envelope"
        topology_pass = root is not None and envelope is not None and root < envelope
        if chi is None or root is None or not (pairwise_pass and dis_pass and topology_pass):
            continue
        accepted.append(
            {
                "transition": "S/C",
                "chiN": chi,
                "fA": root,
                "bracket_width": row["bracket_width"],
                "status": "accepted",
                "acceptance_basis": basis,
                "source": str(path),
            }
        )
    return accepted


def merged_sdis_rows(path: Path, refinement_path: Path) -> list[dict[str, str]]:
    merged: dict[float, dict[str, str]] = {}
    for source in (path, refinement_path):
        for row in read_rows(source):
            merged[float(row["chiN"])] = {**row, "_source": str(source)}
    return [merged[chi] for chi in sorted(merged)]


def accepted_sdis(path: Path, refinement_path: Path) -> list[dict[str, object]]:
    accepted: list[dict[str, object]] = []
    for row in merged_sdis_rows(path, refinement_path):
        chi = finite_float(row.get("chiN"))
        root = finite_float(row.get("fA_root"))
        numerical_pass = all(
            truth(row.get(field))
            for field in ("bracketed", "resolution_reached", "orientation_ok", "gates_pass")
        ) and row.get("blocker", "none") == "none"
        if chi is None or root is None or not numerical_pass:
            continue
        accepted.append(
            {
                "transition": "S/DIS",
                "chiN": chi,
                "fA": root,
                "bracket_width": row["bracket_width"],
                "status": "accepted",
                "acceptance_basis": "ODT from validated BCC/DIS bracket; ordered-side gates; BCC is first stable ordered phase",
                "source": row["_source"],
            }
        )
    return accepted


def accepted_sdis_fixed_fa(path: Path) -> list[dict[str, object]]:
    accepted: list[dict[str, object]] = []
    for row in read_rows(path):
        f_a = finite_float(row.get("root_fA"))
        chi = finite_float(row.get("root_chiN"))
        numerical_pass = row.get("status") == "accepted" and all(
            truth(row.get(field))
            for field in ("bracketed", "resolution_reached", "orientation_ok", "gates_pass")
        ) and row.get("blocker", "none") == "none"
        if f_a is None or chi is None or not numerical_pass:
            continue
        accepted.append(
            {
                "transition": "S/DIS",
                "chiN": chi,
                "fA": f_a,
                "bracket_width": row["bracket_width"],
                "status": "accepted",
                "acceptance_basis": "fixed-fA ODT from validated BCC/DIS predicate bracket; production-grid ordered endpoint gates",
                "source": str(path),
            }
        )
    return accepted


def complete_fixed_fa_extension(rows: list[dict[str, object]]) -> bool:
    targets = {0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49}
    computed = {round(float(row["fA"]), 2) for row in rows}
    return targets <= computed


def critical_sdis_endpoint(source: Path) -> dict[str, object]:
    return {
        "transition": "S/DIS",
        "chiN": CRITICAL_POINT[1],
        "fA": CRITICAL_POINT[0],
        "bracket_width": "",
        "status": "critical_endpoint",
        "acceptance_basis": "A/B-symmetry critical endpoint shared by the completed ODT",
        "source": str(source),
    }


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    fields = (
        "transition",
        "chiN",
        "fA",
        "bracket_width",
        "claim_kind",
        "status",
        "acceptance_basis",
        "source",
    )
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows({
            **row,
            "claim_kind": row.get("claim_kind", "boundary_bracket"),
        } for row in rows)


def read_odt_approximation(path: Path) -> list[tuple[float, float]]:
    return sorted(
        (
            (float(row["fA"]), float(row["chiN"]))
            for row in read_rows(path)
            if row["transition"] == "L/DIS" and float(row["fA"]) >= 0.3
        ),
        key=lambda point: point[1],
        reverse=True,
    )


def remove_odt_overlap(
    points: list[tuple[float, float]], sdis: list[dict[str, object]]
) -> list[tuple[float, float]]:
    if not sdis:
        return points
    computed_limit = max(float(row["fA"]) for row in sdis)
    return [(f_a, chi_n) for f_a, chi_n in points if f_a > computed_limit]


def write_odt_approximation(path: Path, points: list[tuple[float, float]], source: Path) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=("transition", "fA", "chiN", "role", "source"), lineterminator="\n")
        writer.writeheader()
        for f_a, chi_n in points:
            writer.writerow(
                {
                    "transition": "L/DIS",
                    "fA": f_a,
                    "chiN": chi_n,
                    "role": "display approximation beyond the computed S/DIS composition range",
                    "source": str(source),
                }
            )


def load_scft_reference(path: Path) -> list[dict[str, object]]:
    reference: list[dict[str, object]] = []
    for row in read_rows(path):
        if row.get("theory") != "SCFT":
            continue
        f_a = finite_float(row.get("f_A"))
        chi_n = finite_float(row.get("chiN"))
        point_index = finite_float(row.get("point_index"))
        if f_a is None or chi_n is None or point_index is None:
            continue
        reference.append(
            {
                "theory": "SCFT",
                "transition": row.get("transition", ""),
                "side": row.get("side", ""),
                "curve_id": row.get("curve_id") or row.get("source_path_id", ""),
                "point_index": int(point_index),
                "fA": f_a,
                "chiN": chi_n,
            }
        )
    return reference


def scft_reference_curves(
    rows: list[dict[str, object]],
) -> list[tuple[str, list[tuple[float, float]]]]:
    grouped: dict[tuple[str, str, str], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        key = (str(row["transition"]), str(row["side"]), str(row["curve_id"]))
        grouped[key].append(row)

    curves: list[tuple[str, list[tuple[float, float]]]] = []
    for (transition, _side, _curve_id), group in grouped.items():
        run: list[tuple[float, float]] = []
        for row in sorted(group, key=lambda item: int(item["point_index"])):
            point = (float(row["fA"]), float(row["chiN"]))
            if 0.0 <= point[0] <= 1.0 and 8.0 <= point[1] <= 64.0:
                run.append(point)
            else:
                if len(run) >= 2:
                    curves.append((transition, run))
                run = []
        if len(run) >= 2:
            curves.append((transition, run))
    return curves


def write_scft_reference_csv(
    path: Path, rows: list[dict[str, object]], source: Path, source_svg: Path,
    model_label: str = "BURP-TI",
) -> None:
    fields = (
        "theory",
        "transition",
        "side",
        "curve_id",
        "point_index",
        "fA",
        "chiN",
        "role",
        "source",
        "source_svg",
    )
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            if not (0.0 <= float(row["fA"]) <= 1.0 and 8.0 <= float(row["chiN"]) <= 64.0):
                continue
            writer.writerow(
                {
                    **row,
                    "role": (
                        "literature reference overlay; not a "
                        f"{model_label} accepted root"
                    ),
                    "source": str(source),
                    "source_svg": str(source_svg),
                }
            )


def write_sc_identity_audit(path: Path, sc_path: Path) -> None:
    campaign = sc_path.parent
    audit: list[dict[str, object]] = []
    for root in read_rows(sc_path):
        if root.get("status") != "accepted":
            continue
        f = float(root["root_fA"])
        chi = float(root["root_chiN"])
        ft = f"{f:.6f}".replace(".", "p")
        ct = f"{chi:.6f}".replace(".", "p")
        for phase in ("BCC", "CYL"):
            density = campaign / "densities" / "S_C" / f"{phase}_f{ft}_chiN{ct}_density.csv"
            evaluation = campaign / "evaluations" / "S_C" / f"{phase}_f{ft}_chiN{ct}.csv"
            if not density.is_file() or not evaluation.is_file():
                if path.is_file():
                    return
                raise FileNotFoundError(
                    "S/C identity evidence is unavailable and no existing audit can be retained: "
                    f"{density}, {evaluation}"
                )
            values: list[float] = []
            with density.open(newline="") as handle:
                reader = csv.DictReader(handle)
                field = "phi" if "phi" in (reader.fieldnames or ()) else (reader.fieldnames or [""])[-1]
                values = [float(row[field]) for row in reader]
            evidence = read_rows(evaluation)[0]
            contrast = max(values) - min(values)
            audit.append(
                {
                    "chiN": chi,
                    "fA": f,
                    "phase": phase,
                    "density_min": min(values),
                    "density_max": max(values),
                    "density_contrast": contrast,
                    "contrast_tolerance": 0.02,
                    "contrast_pass": str(contrast > 0.02).lower(),
                    "kkt_force_norm": evidence["kkt_force_norm"],
                    "cell_gradient_density": evidence["cell_gradient_density"],
                    "field_accepted": evidence["field_accepted"],
                    "cell_accepted": evidence["cell_accepted"],
                    "density_path": str(density),
                    "evaluation_path": str(evaluation),
                }
            )
    fields = list(audit[0])
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(audit)


def boundary_geometry(
    rows: list[dict[str, object]], transition: str,
    odt_approximation: list[tuple[float, float]],
    *, critical_endpoint_exclusions: tuple[str, ...] = (),
) -> tuple[list[list[tuple[float, float]]], list[tuple[float, float]]]:
    computed = sorted(
        ((float(row["fA"]), float(row["chiN"])) for row in rows if row["transition"] == transition),
        key=lambda point: point[1],
    )
    if transition == "S/DIS" and odt_approximation:
        # The L/DIS trace is only a display approximation in the unresolved
        # near-critical composition range.  It must not be folded into the
        # solid, symbolized BCC/DIS boundary or its acceptance ledger.
        computed = sorted(
            (point for point in computed if abs(point[0] - 0.5) > 1.0e-9),
            key=lambda point: point[0],
        )
        left_curve = computed
        right_curve = [
            (1.0 - f_a, chi_n) for f_a, chi_n in reversed(computed)
        ]
        curves = [left_curve, right_curve]
    else:
        left_curve = list(reversed(computed))
        right_curve = [(1.0 - f_a, chi_n) for f_a, chi_n in computed]
        if (
            transition != "S/DIS"
            and computed
            and transition not in critical_endpoint_exclusions
        ):
            critical = next(
                (
                    (float(row["fA"]), float(row["chiN"]))
                    for row in rows
                    if row["transition"] == "S/DIS"
                    and abs(float(row["fA"]) - 0.5) <= 1.0e-9
                ),
                None,
            )
            if critical is not None:
                left_curve.append(critical)
                right_curve.insert(0, critical)
        curves = [left_curve, right_curve]
    return curves, computed


def odt_bridge(rows: list[dict[str, object]], odt_approximation: list[tuple[float, float]]) -> list[tuple[float, float]]:
    computed = [
        (float(row["fA"]), float(row["chiN"]))
        for row in rows
        if row["transition"] == "S/DIS"
        and row.get("status") == "accepted"
        and float(row["fA"]) < 0.5 - 1.0e-9
    ]
    if not computed or not odt_approximation:
        return []
    join = max(computed, key=lambda point: point[0])
    approximate_left = sorted(
        {
            (float(f_a), float(chi_n))
            for f_a, chi_n in odt_approximation
            if join[0] < float(f_a) <= 0.5 + 1.0e-9
        },
        key=lambda point: point[0],
    )
    if not approximate_left:
        return []
    left = [join, *approximate_left]
    return [*left, *((1.0 - f_a, chi_n) for f_a, chi_n in reversed(left[:-1]))]


def svg_marker(kind: str, x: float, y: float, color: str) -> str:
    common = f"class='computed-point' fill='white' stroke='{color}' stroke-width='1.6'"
    if kind == "circle":
        return f"<circle cx='{x:.2f}' cy='{y:.2f}' r='4.2' {common}/>"
    if kind == "square":
        return f"<rect x='{x-3.7:.2f}' y='{y-3.7:.2f}' width='7.4' height='7.4' {common}/>"
    if kind == "triangle":
        points = f"{x:.2f},{y-4.5:.2f} {x-4.1:.2f},{y+3.7:.2f} {x+4.1:.2f},{y+3.7:.2f}"
        return f"<polygon points='{points}' {common}/>"
    if kind == "diamond":
        points = f"{x:.2f},{y-4.4:.2f} {x-4.4:.2f},{y:.2f} {x:.2f},{y+4.4:.2f} {x+4.4:.2f},{y:.2f}"
        return f"<polygon points='{points}' {common}/>"
    return (
        f"<g class='computed-point' fill='none' stroke='{color}' stroke-width='1.8'>"
        f"<line x1='{x-4.2:.2f}' y1='{y:.2f}' x2='{x+4.2:.2f}' y2='{y:.2f}'/>"
        f"<line x1='{x:.2f}' y1='{y-4.2:.2f}' x2='{x:.2f}' y2='{y+4.2:.2f}'/></g>"
    )


def star_polygon(
    x: float, y: float, outer_radius: float = 5.0, inner_radius: float = 2.2
) -> tuple[tuple[float, float], ...]:
    return tuple(
        (
            x + (outer_radius if index % 2 == 0 else inner_radius)
            * cos(-pi / 2 + index * pi / 5),
            y + (outer_radius if index % 2 == 0 else inner_radius)
            * sin(-pi / 2 + index * pi / 5),
        )
        for index in range(10)
    )


def rows_with_estimated_triple_points(
    rows: list[dict[str, object]],
    triple_points: tuple[EstimatedTriplePoint, ...],
) -> list[dict[str, object]]:
    plotted = list(rows)
    for f_a, chi_n, connected_transitions in triple_points:
        plotted.extend(
            {
                "transition": transition,
                "fA": f_a,
                "chiN": chi_n,
                "status": "accepted",
                "estimated_display_only": True,
            }
            for transition in connected_transitions
        )
    return plotted


def rows_in_chi_window(
    rows: list[dict[str, object]], chi_n_min: float, chi_n_max: float
) -> list[dict[str, object]]:
    return [
        row for row in rows
        if chi_n_min <= float(row["chiN"]) <= chi_n_max
    ]


def accepted_legend_transitions(
    rows: list[dict[str, object]], transitions: tuple[str, ...]
) -> tuple[str, ...]:
    accepted = {
        str(row.get("transition"))
        for row in rows
        if row.get("status") == "accepted"
    }
    return tuple(transition for transition in transitions if transition in accepted)


def transition_display_tier(
    rows: list[dict[str, object]], transition: str
) -> str:
    """Return the least mature evidence tier displayed for a transition."""
    statuses = {
        str(row.get("status"))
        for row in rows
        if row.get("transition") == transition
    }
    if "provisional" in statuses:
        return "provisional"
    if "fixed_cell" in statuses:
        return "fixed_cell"
    if "confirmed" in statuses:
        return "confirmed"
    if "accepted" in statuses:
        return "accepted"
    return ""


def transition_tier_rows(
    rows: list[dict[str, object]], transition: str, tier: str
) -> list[dict[str, object]]:
    """Return one evidence tier without connecting it to another grid tier."""
    selected = [
        row for row in rows
        if row.get("transition") == transition and row.get("status") == tier
    ]
    if tier == "accepted" and transition == "S/DIS":
        selected.extend(
            row for row in rows
            if row.get("transition") == transition
            and row.get("status") == "critical_endpoint"
        )
    elif transition != "S/DIS":
        selected.extend(
            row for row in rows
            if row.get("transition") == "S/DIS"
            and row.get("status") == "critical_endpoint"
        )
    return selected


def transition_tier_segments(
    rows: list[dict[str, object]], transition: str, tier: str
) -> list[list[dict[str, object]]]:
    """Split one tier wherever another visible tier intervenes.

    Connecting all accepted rows after filtering out provisional rows makes a
    solid curve bridge an unresolved interval.  Preserve the full chiN order
    first, then form only contiguous runs of the requested tier.  The run at
    the low-chi end alone inherits the shared critical endpoint.
    """
    visible = {"accepted", "confirmed", "fixed_cell", "provisional"}
    ordered = sorted(
        (
            row for row in rows
            if row.get("transition") == transition
            and row.get("status") in visible
        ),
        key=lambda row: float(row["chiN"]),
    )
    slices: list[list[dict[str, object]]] = []
    for row in ordered:
        if (
            not slices
            or float(slices[-1][0]["chiN"]) != float(row["chiN"])
        ):
            slices.append([row])
        else:
            slices[-1].append(row)
    segments: list[list[dict[str, object]]] = []
    current: list[dict[str, object]] = []
    for slice_rows in slices:
        matching = [
            row for row in slice_rows if row.get("status") == tier
        ]
        if matching:
            current.extend(matching)
        elif current:
            segments.append(current)
            current = []
    if current:
        segments.append(current)

    if (
        transition != "S/DIS"
        and segments
        and slices
        and any(row.get("status") == tier for row in slices[0])
    ):
        critical = [
            row for row in rows
            if row.get("transition") == "S/DIS"
            and row.get("status") == "critical_endpoint"
        ]
        segments[0] = [*segments[0], *critical]
    return segments


def displayed_transition_tiers(
    rows: list[dict[str, object]], transition: str
) -> tuple[str, ...]:
    return tuple(
        tier for tier in ("accepted", "confirmed", "fixed_cell", "provisional")
        if any(
            row.get("transition") == transition and row.get("status") == tier
            for row in rows
        )
    )


def displayed_legend_transitions(
    rows: list[dict[str, object]], transitions: tuple[str, ...]
) -> tuple[str, ...]:
    return tuple(
        transition for transition in transitions
        if transition_display_tier(rows, transition)
    )


def displayed_legend_entries(
    rows: list[dict[str, object]], transitions: tuple[str, ...]
) -> tuple[tuple[str, str], ...]:
    """List every visible transition/tier so mixed grids stay explicit."""
    return tuple(
        (transition, tier)
        for transition in transitions
        for tier in displayed_transition_tiers(rows, transition)
    )


def render_svg(
    path: Path,
    rows: list[dict[str, object]],
    odt_approximation: list[tuple[float, float]],
    scft_reference: list[dict[str, object]],
    title: str = "BURP-TI continuous phase diagram",
    transitions: tuple[str, ...] = TRANSITIONS,
    *,
    critical_endpoint_exclusions: tuple[str, ...] = (),
    estimated_triple_points: tuple[EstimatedTriplePoint, ...] = (),
    chi_n_max: float = 64.0,
    rpa_stability_limit: tuple[tuple[float, float], ...] = (),
) -> None:
    width, height = 1000, 620
    left, right, top, bottom = 70.0, 680.0, 40.0, 560.0
    xmin, xmax, ymin, ymax = 0.0, 1.0, 8.0, chi_n_max

    def sx(value: float) -> float:
        return left + (value - xmin) / (xmax - xmin) * (right - left)

    def sy(value: float) -> float:
        return bottom - (value - ymin) / (ymax - ymin) * (bottom - top)

    plotted_rows = rows_in_chi_window(
        rows_with_estimated_triple_points(rows, estimated_triple_points),
        ymin,
        ymax,
    )

    parts = [
        f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' viewBox='0 0 {width} {height}'>",
        "<rect width='100%' height='100%' fill='white'/>",
        "<style>text{font-family:'DejaVu Sans',sans-serif}</style>",
        f"<defs><clipPath id='plot-clip'><rect x='{left}' y='{top}' width='{right-left}' height='{bottom-top}'/></clipPath></defs>",
    ]
    for f_a in (0.0, 0.2, 0.4, 0.6, 0.8, 1.0):
        x = sx(f_a)
        parts.append(f"<text x='{x:.2f}' y='{bottom+20}' text-anchor='middle' font-size='12'>{f_a:.1f}</text>")
    for chi_n in range(10, int(ymax) + 1, 10):
        y = sy(chi_n)
        parts.append(f"<text x='{left-10}' y='{y+4:.2f}' text-anchor='end' font-size='12'>{chi_n}</text>")
    center = sx(0.5)
    parts.append(f"<line x1='{center:.2f}' y1='{top}' x2='{center:.2f}' y2='{bottom}' stroke='#9ca3af' stroke-width='1' stroke-dasharray='4 4'/>")
    parts.append(f"<rect x='{left}' y='{top}' width='{right-left}' height='{bottom-top}' fill='none' stroke='#111827'/>")
    rpa_left = [
        point for point in sorted(rpa_stability_limit)
        if ymin <= point[1] <= ymax
    ]
    if rpa_left:
        rpa_curve = [
            *rpa_left,
            *((1.0 - f_a, chi_n) for f_a, chi_n in reversed(rpa_left[:-1])),
        ]
        points = " ".join(
            f"{sx(f_a):.2f},{sy(chi_n):.2f}" for f_a, chi_n in rpa_curve
        )
        parts.append(
            "<polyline class='rpa-stability-limit' "
            f"points='{points}' fill='none' stroke='{RPA_COLOR}' "
            "stroke-width='1.8' stroke-dasharray='1 3' stroke-linecap='round' "
            "clip-path='url(#plot-clip)'/>"
        )
    for transition, curve in scft_reference_curves(scft_reference):
        curve = [point for point in curve if ymin <= point[1] <= ymax]
        points = " ".join(f"{sx(f_a):.2f},{sy(chi_n):.2f}" for f_a, chi_n in curve)
        parts.append(
            f"<polyline class='scft-literature-overlay' data-transition='{transition}' "
            f"points='{points}' fill='none' stroke='{SCFT_COLOR}' stroke-width='1.6' "
            "stroke-dasharray='6 4' clip-path='url(#plot-clip)'/>"
        )
    for transition in transitions:
        color = COLORS[transition]
        for tier in displayed_transition_tiers(plotted_rows, transition):
            for tier_rows in transition_tier_segments(
                plotted_rows, transition, tier
            ):
                curves, computed = boundary_geometry(
                    tier_rows, transition, odt_approximation,
                    critical_endpoint_exclusions=critical_endpoint_exclusions,
                )
                for curve_index, curve in enumerate(curves):
                    if len(curve) >= 2:
                        side = "computed" if curve_index == 0 else "mirrored"
                        points = " ".join(f"{sx(f_a):.2f},{sy(chi_n):.2f}" for f_a, chi_n in curve)
                        css_class = (
                            "computed-boundary provisional-boundary"
                            if tier == "provisional"
                            else "computed-boundary fixed-cell-boundary"
                            if tier == "fixed_cell"
                            else "computed-boundary confirmed-boundary"
                            if tier == "confirmed"
                            else "computed-boundary"
                        )
                        dash = (
                            ""
                            if side == "mirrored"
                            else
                            " stroke-dasharray='2 4'"
                            if tier == "provisional"
                            else " stroke-dasharray='3 3'"
                            if tier == "fixed_cell"
                            else " stroke-dasharray='5 3'"
                            if tier == "confirmed"
                            else ""
                        )
                        parts.append(
                            f"<polyline class='{css_class}' data-transition='{transition}' "
                            f"data-side='{side}' "
                            f"points='{points}' fill='none' stroke='{color}' stroke-width='2.4' "
                            f"{dash.strip()}/>")
        for row in plotted_rows:
            if (
                row["transition"] == transition
                and row["status"] in {
                    "accepted", "confirmed", "fixed_cell", "provisional"
                }
                and not row.get("estimated_display_only")
                and not (transition == "S/DIS" and odt_approximation)
            ):
                parts.append(svg_marker(MARKERS[transition], sx(float(row["fA"])), sy(float(row["chiN"])), color))
    for f_a, chi_n, _connected_transitions in estimated_triple_points:
        for display_f_a in (f_a, 1.0 - f_a):
            points = " ".join(
                f"{x:.2f},{y:.2f}"
                for x, y in star_polygon(sx(display_f_a), sy(chi_n))
            )
            parts.append(
                f"<polygon class='estimated-triple-point' points='{points}' "
                f"fill='{ESTIMATED_POINT_COLOR}'/>"
            )
    bridge = odt_bridge(rows, odt_approximation)
    if bridge:
        points = " ".join(f"{sx(f_a):.2f},{sy(chi_n):.2f}" for f_a, chi_n in bridge)
        parts.append(f"<polyline class='odt-approximation' points='{points}' fill='none' stroke='{ODT_COLOR}' stroke-width='2.4' stroke-dasharray='7 4'/>")
    # Draw direct BCC/DIS evidence above the dashed approximation.  The
    # approximation, reflected branch, and critical endpoint stay symbol-free;
    # construction views may additionally expose provisional left-half roots.
    if odt_approximation:
        for row in rows:
            if row["transition"] == "S/DIS" and row["status"] in {
                "accepted", "confirmed", "fixed_cell", "provisional"
            }:
                parts.append(svg_marker("square", sx(float(row["fA"])), sy(float(row["chiN"])), COLORS["S/DIS"]))
    parts.extend(
        [
            f"<text x='{(left+right)/2}' y='{bottom+48}' text-anchor='middle' font-size='14'>f_A</text>",
            f"<text transform='translate(22 {(top+bottom)/2}) rotate(-90)' text-anchor='middle' font-size='14'>&#967;N</text>",
            f"<text x='{(left+right)/2}' y='{bottom+73}' text-anchor='middle' font-size='11' fill='#4b5563'>A/B symmetry shown by reflection about f_A = 0.5</text>",
        ]
    )
    legend_entries = displayed_legend_entries(rows, transitions)
    for index, (transition, tier) in enumerate(legend_entries):
        y = 92 + 34 * index
        dash = (
            " stroke-dasharray='2 4'"
            if tier == "provisional"
            else " stroke-dasharray='3 3'"
            if tier == "fixed_cell"
            else " stroke-dasharray='5 3'"
            if tier == "confirmed"
            else ""
        )
        label = (
            f"{transition} (provisional)"
            if tier == "provisional"
            else f"{transition} (fixed-cell)"
            if tier == "fixed_cell"
            else f"{transition} (confirmed)"
            if tier == "confirmed"
            else transition
        )
        parts.append(f"<line x1='710' y1='{y}' x2='750' y2='{y}' stroke='{COLORS[transition]}' stroke-width='2.4'{dash}/>")
        parts.append(svg_marker(MARKERS[transition], 730, y, COLORS[transition]))
        parts.append(f"<text x='762' y='{y+4}' font-size='13'>{label}</text>")
    if odt_approximation:
        y = 92 + 34 * len(legend_entries)
        parts.append(f"<line x1='710' y1='{y}' x2='750' y2='{y}' stroke='{ODT_COLOR}' stroke-width='2.4' stroke-dasharray='7 4'/>")
        parts.append(f"<text x='762' y='{y+4}' font-size='13'>L/DIS approx.</text>")
    legend_offset = len(legend_entries) + int(bool(odt_approximation))
    if estimated_triple_points:
        y = 92 + 34 * legend_offset
        points = " ".join(
            f"{x:.2f},{point_y:.2f}"
            for x, point_y in star_polygon(730, y)
        )
        parts.append(
            f"<polygon class='estimated-triple-point-legend' points='{points}' "
            f"fill='{ESTIMATED_POINT_COLOR}'/>"
        )
        parts.append(
            f"<text x='762' y='{y+4}' font-size='13'>estimated triple point</text>"
        )
        legend_offset += 1
    if rpa_stability_limit:
        y = 92 + 34 * legend_offset
        parts.append(
            f"<line x1='710' y1='{y}' x2='750' y2='{y}' "
            f"stroke='{RPA_COLOR}' stroke-width='1.8' "
            "stroke-dasharray='1 3' stroke-linecap='round'/>")
        parts.append(
            f"<text x='762' y='{y+4}' font-size='13'>RPA stability limit</text>"
        )
        legend_offset += 1
    y = 92 + 34 * legend_offset
    parts.append(f"<line x1='710' y1='{y}' x2='750' y2='{y}' stroke='{SCFT_COLOR}' stroke-width='1.6' stroke-dasharray='6 4'/>")
    parts.append(f"<text x='762' y='{y+4}' font-size='13'>SCFT literature</text>")
    parts.append("</svg>")
    path.write_text("\n".join(parts) + "\n")


def draw_marker(draw: ImageDraw.ImageDraw, kind: str, x: float, y: float, color: str) -> None:
    radius = 4.2
    if kind == "circle":
        draw.ellipse((x-radius, y-radius, x+radius, y+radius), fill="white", outline=color, width=2)
    elif kind == "square":
        draw.rectangle((x-radius, y-radius, x+radius, y+radius), fill="white", outline=color, width=2)
    elif kind == "triangle":
        draw.polygon(((x, y-4.5), (x-4.1, y+3.7), (x+4.1, y+3.7)), fill="white", outline=color)
        draw.line(((x, y-4.5), (x-4.1, y+3.7), (x+4.1, y+3.7), (x, y-4.5)), fill=color, width=2)
    elif kind == "diamond":
        draw.polygon(((x, y-4.4), (x-4.4, y), (x, y+4.4), (x+4.4, y)), fill="white", outline=color)
        draw.line(((x, y-4.4), (x-4.4, y), (x, y+4.4), (x+4.4, y), (x, y-4.4)), fill=color, width=2)
    else:
        draw.line((x-4.2, y, x+4.2, y), fill=color, width=2)
        draw.line((x, y-4.2, x, y+4.2), fill=color, width=2)


def draw_star(draw: ImageDraw.ImageDraw, x: float, y: float) -> None:
    draw.polygon(star_polygon(x, y), fill=ESTIMATED_POINT_COLOR)


def draw_dashed_polyline(
    draw: ImageDraw.ImageDraw,
    points: list[tuple[float, float]],
    color: str,
    width: int = 3,
    dash_length: float = 7.0,
    gap_length: float = 4.0,
) -> None:
    draw_dash = True
    interval_remaining = dash_length
    for (x0, y0), (x1, y1) in zip(points, points[1:]):
        length = hypot(x1 - x0, y1 - y0)
        if length == 0.0:
            continue
        position = 0.0
        while position < length:
            step = min(interval_remaining, length - position)
            end = position + step
            if draw_dash:
                a, b = position / length, end / length
                draw.line(
                    (x0 + a * (x1 - x0), y0 + a * (y1 - y0), x0 + b * (x1 - x0), y0 + b * (y1 - y0)),
                    fill=color,
                    width=width,
                )
            position = end
            interval_remaining -= step
            if interval_remaining <= 1.0e-9:
                draw_dash = not draw_dash
                interval_remaining = dash_length if draw_dash else gap_length


def render_png(
    path: Path,
    rows: list[dict[str, object]],
    odt_approximation: list[tuple[float, float]],
    scft_reference: list[dict[str, object]],
    title: str = "BURP-TI continuous phase diagram",
    transitions: tuple[str, ...] = TRANSITIONS,
    *,
    critical_endpoint_exclusions: tuple[str, ...] = (),
    estimated_triple_points: tuple[EstimatedTriplePoint, ...] = (),
    chi_n_max: float = 64.0,
    rpa_stability_limit: tuple[tuple[float, float], ...] = (),
) -> None:
    width, height = 1000, 620
    left, right, top, bottom = 70.0, 680.0, 40.0, 560.0
    xmin, xmax, ymin, ymax = 0.0, 1.0, 8.0, chi_n_max
    sx = lambda value: left + (value - xmin) / (xmax - xmin) * (right - left)
    sy = lambda value: bottom - (value - ymin) / (ymax - ymin) * (bottom - top)
    plotted_rows = rows_in_chi_window(
        rows_with_estimated_triple_points(rows, estimated_triple_points),
        ymin,
        ymax,
    )
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)
    font = ImageFont.truetype("DejaVuSans.ttf", 12)
    small_font = ImageFont.truetype("DejaVuSans.ttf", 10)
    for f_a in (0.0, 0.2, 0.4, 0.6, 0.8, 1.0):
        x = sx(f_a)
        draw.text((x, bottom + 15), f"{f_a:.1f}", fill="#111827", anchor="mm", font=font)
    for chi_n in range(10, int(ymax) + 1, 10):
        y = sy(chi_n)
        draw.text((left - 8, y), str(chi_n), fill="#111827", anchor="rm", font=font)
    center = sx(0.5)
    for y in range(int(top), int(bottom), 8):
        draw.line((center, y, center, min(y + 4, bottom)), fill="#9ca3af")
    draw.rectangle((left, top, right, bottom), outline="#111827")
    rpa_left = [
        point for point in sorted(rpa_stability_limit)
        if ymin <= point[1] <= ymax
    ]
    if rpa_left:
        rpa_curve = [
            *rpa_left,
            *((1.0 - f_a, chi_n) for f_a, chi_n in reversed(rpa_left[:-1])),
        ]
        draw_dashed_polyline(
            draw,
            [(sx(f_a), sy(chi_n)) for f_a, chi_n in rpa_curve],
            RPA_COLOR,
            width=2,
            dash_length=1.0,
            gap_length=3.0,
        )
    for _transition, curve in scft_reference_curves(scft_reference):
        curve = [point for point in curve if ymin <= point[1] <= ymax]
        points = [(sx(f_a), sy(chi_n)) for f_a, chi_n in curve]
        draw_dashed_polyline(draw, points, SCFT_COLOR, width=2, dash_length=6.0, gap_length=4.0)
    for transition in transitions:
        for tier in displayed_transition_tiers(plotted_rows, transition):
            for tier_rows in transition_tier_segments(
                plotted_rows, transition, tier
            ):
                curves, computed = boundary_geometry(
                    tier_rows, transition, odt_approximation,
                    critical_endpoint_exclusions=critical_endpoint_exclusions,
                )
                for curve_index, curve in enumerate(curves):
                    mirrored = curve_index == 1
                    points = [(sx(f_a), sy(chi_n)) for f_a, chi_n in curve]
                    if len(points) >= 2:
                        if (
                            not mirrored
                            and tier in {
                                "confirmed", "fixed_cell", "provisional"
                            }
                        ):
                            draw_dashed_polyline(
                                draw, points, COLORS[transition], width=3,
                                dash_length=(
                                    2.0 if tier == "provisional"
                                    else 3.0 if tier == "fixed_cell" else 5.0
                                ),
                                gap_length=(
                                    4.0 if tier == "provisional" else 3.0
                                ),
                            )
                        else:
                            draw.line(
                                points, fill=COLORS[transition], width=3,
                                joint="curve",
                            )
        for row in plotted_rows:
            if (
                row["transition"] == transition
                and row["status"] in {
                    "accepted", "confirmed", "fixed_cell", "provisional"
                }
                and not row.get("estimated_display_only")
                and not (transition == "S/DIS" and odt_approximation)
            ):
                draw_marker(draw, MARKERS[transition], sx(float(row["fA"])), sy(float(row["chiN"])), COLORS[transition])
    for f_a, chi_n, _connected_transitions in estimated_triple_points:
        draw_star(draw, sx(f_a), sy(chi_n))
        draw_star(draw, sx(1.0 - f_a), sy(chi_n))
    bridge = [(sx(f_a), sy(chi_n)) for f_a, chi_n in odt_bridge(rows, odt_approximation)]
    draw_dashed_polyline(draw, bridge, ODT_COLOR)
    if odt_approximation:
        for row in rows:
            if row["transition"] == "S/DIS" and row["status"] in {
                "accepted", "confirmed", "fixed_cell", "provisional"
            }:
                draw_marker(draw, "square", sx(float(row["fA"])), sy(float(row["chiN"])), COLORS["S/DIS"])
    draw.text(((left + right) / 2, bottom + 42), "f_A", fill="#111827", anchor="mm", font=font)
    y_label = Image.new("RGBA", (90, 30), (255, 255, 255, 0))
    ImageDraw.Draw(y_label).text((45, 15), "χN", fill="#111827", anchor="mm", font=font)
    y_label = y_label.rotate(90, expand=True)
    image.paste(y_label, (4, int((top + bottom - y_label.height) / 2)), y_label)
    draw.text(((left + right) / 2, bottom + 64), "A/B symmetry shown by reflection about f_A = 0.5", fill="#4b5563", anchor="mm", font=small_font)
    legend_entries = displayed_legend_entries(rows, transitions)
    for index, (transition, tier) in enumerate(legend_entries):
        y = 92 + 34 * index
        if tier in {"confirmed", "fixed_cell", "provisional"}:
            draw_dashed_polyline(
                draw, [(710, y), (750, y)], COLORS[transition], width=3,
                dash_length=(
                    2.0 if tier == "provisional"
                    else 3.0 if tier == "fixed_cell" else 5.0
                ),
                gap_length=4.0 if tier == "provisional" else 3.0,
            )
        else:
            draw.line((710, y, 750, y), fill=COLORS[transition], width=3)
        draw_marker(draw, MARKERS[transition], 730, y, COLORS[transition])
        label = (
            f"{transition} (provisional)"
            if tier == "provisional"
            else f"{transition} (fixed-cell)"
            if tier == "fixed_cell"
            else f"{transition} (confirmed)"
            if tier == "confirmed"
            else transition
        )
        draw.text((762, y), label, fill="#111827", anchor="lm", font=font)
    if odt_approximation:
        y = 92 + 34 * len(legend_entries)
        for x in range(710, 750, 10):
            draw.line((x, y, min(x + 6, 750), y), fill=ODT_COLOR, width=3)
        draw.text((762, y), "L/DIS approx.", fill="#111827", anchor="lm", font=font)
    legend_offset = len(legend_entries) + int(bool(odt_approximation))
    if estimated_triple_points:
        y = 92 + 34 * legend_offset
        draw_star(draw, 730, y)
        draw.text(
            (762, y), "estimated triple point",
            fill="#111827", anchor="lm", font=font,
        )
        legend_offset += 1
    if rpa_stability_limit:
        y = 92 + 34 * legend_offset
        draw_dashed_polyline(
            draw, [(710, y), (750, y)], RPA_COLOR, width=2,
            dash_length=1.0, gap_length=3.0,
        )
        draw.text(
            (762, y), "RPA stability limit",
            fill="#111827", anchor="lm", font=font,
        )
        legend_offset += 1
    y = 92 + 34 * legend_offset
    for x in range(710, 750, 10):
        draw.line((x, y, min(x + 6, 750), y), fill=SCFT_COLOR, width=2)
    draw.text((762, y), "SCFT literature", fill="#111827", anchor="lm", font=font)
    image.save(path)


def write_summary(
    path: Path,
    rows: list[dict[str, object]],
    odt_approximation: list[tuple[float, float]],
    scft_reference: list[dict[str, object]],
    scft_source: Path,
    scft_source_svg: Path,
) -> None:
    counts = {transition: sum(row["transition"] == transition for row in rows) for transition in COLORS}
    sc_chi = {float(row["chiN"]) for row in rows if row["transition"] == "S/C"}
    sdis_chi = {float(row["chiN"]) for row in rows if row["transition"] == "S/DIS"}
    closure_values = [float(row["chiN"]) for row in rows if "sphere-pocket closure" in str(row["acceptance_basis"])]
    closure_chi = min(closure_values) if closure_values else -float("inf")
    pending_sc = [chi for chi in EXPECTED_SC_CHI if chi >= closure_chi and chi not in sc_chi]
    pending_sdis = [chi for chi in EXPECTED_SDIS_CHI if chi >= closure_chi and chi not in sdis_chi]
    scft_curves = scft_reference_curves(scft_reference)
    lines = [
        "# BURP-TI publication phase diagram",
        "",
        "This assembly contains accepted roots plus the shared A/B-symmetry critical endpoint; incomplete traces remain gaps. The ODT is defined by the validated S/DIS boundary because BCC is the first stable ordered phase above disorder. All reliable BCC/DIS roots are connected as one solid S/DIS curve. The computed f_A <= 0.5 half is reflected about f_A = 0.5 using A/B diblock symmetry, with symbols only on computed roots. " + ("The previous L/DIS approximation remains only beyond the computed S/DIS range." if odt_approximation else "The fixed-composition S/DIS extension reaches f_A = 0.49 and connects to the critical endpoint, so no L/DIS approximation is used."),
        "",
        "| boundary | boundary entries |",
        "| --- | ---: |",
        *[f"| {transition} | {counts[transition]} |" for transition in ("S/DIS", "S/C", "G/C", "G/L", "C/L")],
        "",
        f"Pending S/C chiN targets: {', '.join(map(str, pending_sc)) or 'none'}.",
        "",
        f"Pending S/DIS chiN targets: {', '.join(map(str, pending_sdis)) or 'none'}.",
        "",
        f"The plot overlays {len(scft_curves)} dashed neutral-gray SCFT literature curves (S/DIS, C/S, L/C, and the low-chiN gyroid-region boundary) from `{scft_source_svg}`. The plotted coordinates are reproduced from `{scft_source}` and exported in `scft_reference_boundaries.csv`. These curves are comparison references only: they carry no symbols and are not counted as accepted BURP-TI roots.",
        "",
        "S/C points require a gate-passing BCC/CYL crossing below the accepted cylinder-side G/C or C/L envelope. S/DIS points require a validated ordered/disordered bracket with the ordered BCC endpoint passing all gates.",
    ]
    path.write_text("\n".join(lines) + "\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, default=DEFAULT_BASE)
    parser.add_argument("--sc", type=Path, default=DEFAULT_SC)
    parser.add_argument("--sdis", type=Path, default=DEFAULT_SDIS)
    parser.add_argument("--sdis-refinement", type=Path, default=DEFAULT_SDIS_REFINEMENT)
    parser.add_argument("--sdis-fixed-fa", type=Path, default=DEFAULT_SDIS_FIXED_FA)
    parser.add_argument("--odt", type=Path, default=DEFAULT_ODT)
    parser.add_argument("--scft-reference", type=Path, default=DEFAULT_SCFT_REFERENCE)
    parser.add_argument("--scft-source-svg", type=Path, default=DEFAULT_SCFT_SOURCE_SVG)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    base = accepted_base(args.base / "stable_phase_boundaries.csv")
    sc = accepted_sc(args.sc, base)
    sdis = accepted_sdis(args.sdis, args.sdis_refinement)
    fixed_fa = accepted_sdis_fixed_fa(args.sdis_fixed_fa)
    if complete_fixed_fa_extension(fixed_fa):
        sdis = [*sdis, *fixed_fa, critical_sdis_endpoint(args.odt)]
    else:
        sdis = [*sdis, *fixed_fa]
    rows = sorted(base + sc + sdis, key=lambda row: (str(row["transition"]), float(row["chiN"])))
    scft_reference = load_scft_reference(args.scft_reference)
    odt_approximation = [] if complete_fixed_fa_extension(fixed_fa) else remove_odt_overlap(read_odt_approximation(args.odt), sdis)
    write_csv(args.outdir / "accepted_phase_boundaries.csv", rows)
    write_odt_approximation(args.outdir / "odt_approximation.csv", odt_approximation, args.odt)
    write_scft_reference_csv(
        args.outdir / "scft_reference_boundaries.csv",
        scft_reference,
        args.scft_reference,
        args.scft_source_svg,
    )
    write_sc_identity_audit(args.outdir / "sc_identity_audit.csv", args.sc)
    render_svg(args.outdir / "phase_diagram.svg", rows, odt_approximation, scft_reference)
    render_png(args.outdir / "phase_diagram.png", rows, odt_approximation, scft_reference)
    write_summary(
        args.outdir / "phase_diagram.md",
        rows,
        odt_approximation,
        scft_reference,
        args.scft_reference,
        args.scft_source_svg,
    )
    print(
        f"accepted boundaries: base={len(base)} S/C={len(sc)} S/DIS={len(sdis)} "
        f"fixed-fA={len(fixed_fa)} SCFT-curves={len(scft_reference_curves(scft_reference))}"
    )


if __name__ == "__main__":
    main()
