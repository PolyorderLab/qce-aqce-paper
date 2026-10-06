#!/usr/bin/env python3
"""Assemble gate-passing BVK2 boundary roots into a publication diagram.

Mirrors the BURP-TI publication assembly and reuses its geometry/render
functions verbatim, so the two diagrams are directly comparable: the
computed ``f_A <= 0.5`` half carries symbols, the right half is the A/B
symmetric reflection drawn without symbols, and the dashed neutral-gray
SCFT literature curves are an orientation overlay that is never counted as
a BVK2 root.

Only rows whose unit reported ``status=accepted`` are admitted, and each is
re-screened here against the campaign gates (bracket width, root residual)
so a loosened unit run cannot quietly promote a point.

Safe to call at any time, including mid-campaign: it renders whatever has
landed so far, which is what makes the continuous re-render loop possible.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import subprocess
from math import isclose, pi
from pathlib import Path

from assemble_burp_ti_publication_phase_diagram import (
    boundary_geometry,
    load_scft_reference,
    render_png,
    scft_reference_curves,
    write_csv,
    write_scft_reference_csv,
)
from phase_boundary_acceptance import (
    ACCEPTED as CLAIM_ACCEPTED,
    BOUNDARY_BRACKET,
    PROVISIONAL as CLAIM_PROVISIONAL,
    REJECTED as CLAIM_REJECTED,
    require_boundary_claim,
)


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent
RESULTS = PROJECT / "results"
DEFAULT_CAMPAIGN = RESULTS / "bvk2_gc_gl_boundary_campaign"
DEFAULT_SDIS = RESULTS / "bvk2_sdis_boundary_campaign"
DEFAULT_SDIS_FIXED_FA_REFINEMENT = (
    RESULTS / "bvk2_sdis_fixed_fa_chi_refinement"
)
DEFAULT_SC_FIXED_FA_REFINEMENT = (
    RESULTS / "bvk2_sc_fixed_fa_chi_continuation"
)
# SCFT overlay: vector paths from Li and Liu, JCP 154, 014903 (2021),
# Fig. 1(a).  This source is preferred for the manuscript because its native
# PostScript paths give smooth, directly calibrated curves and retain the
# close-packed-sphere and O70/Fddd pockets without raster tracing.
DEFAULT_SCFT_REFERENCE = (
    RESULTS / "jcp2021_ab_phase_diagram_digitization"
    / "jcp2021_ab_phase_boundaries.csv"
)
# Original vector source used by the extraction audit.
DEFAULT_SCFT_SOURCE_SVG = (
    PROJECT.parent / "references" / "AB_phase_diagram_fig1_jcp2021.eps"
)
DEFAULT_OUTDIR = RESULTS / "bvk2_publication_phase_diagram"
DEFAULT_POLYORDER_PROJECT = PROJECT / "reproducibility" / "polyorder"
DEFAULT_RPA_SCRIPT = HERE / "write_bvk2_rpa_stability_limit.jl"
DEFAULT_CORRECTIONS = RESULTS / "bvk2_richardson_correction"
DEFAULT_UD_THETA_SC = (
    RESULTS / "bvk2_ud_theta_sc_highchi" / "accepted_boundaries.csv"
)
DEFAULT_UD_THETA_GC = (
    RESULTS / "bvk2_ud_theta_g_pocket_highchi_ledger"
    / "accepted_gc_boundaries.csv"
)
DEFAULT_UD_THETA_GC_DIRECT = (
    RESULTS / "bvk2_ud_theta_gc_chi54_direct" / "boundary_point.csv"
)
DEFAULT_UD_THETA_GC_FILTERED_R6 = (
    RESULTS / "bvk2_ud_theta_gc_filtered_r6" / "boundary_points.csv"
)
DEFAULT_UD_THETA_GL = (
    RESULTS / "bvk2_ud_theta_g_pocket_highchi_ledger"
    / "accepted_gl_boundaries.csv"
)
DEFAULT_UD_THETA_GL_DIRECT = (
    RESULTS / "bvk2_ud_theta_gl_direct_fine"
    / "boundary_candidates.csv"
)
DEFAULT_UD_THETA_GL_MATCHED = (
    RESULTS / "bvk2_ud_theta_gl_matched_spacing_audit"
    / "matched_boundary_promotions.csv"
)
DEFAULT_O70_GL_CHI15_AUDIT = (
    RESULTS / "bvk2_o70_gl_triple_chi15_matched"
    / "thermodynamic_audit.csv"
)
DEFAULT_UD_THETA_GL_LOWCHI = RESULTS / "bvk2_ud_theta_gl_lowchi_direct"
DEFAULT_FCC_O70 = (
    RESULTS / "bvk2_fcc_o70_boundary_campaign" / "boundary_brackets.csv"
)
DEFAULT_STRICT_FCC = (
    RESULTS / "bvk2_fcc_o70_boundary_campaign"
    / "matched_dx_fcc_boundary_brackets.csv"
)
DEFAULT_STRICT_FCC_POINTS = (
    RESULTS / "bvk2_fcc_o70_boundary_campaign"
    / "matched_dx_fcc_boundary_points.csv"
)

PAIRWISE_TRANSITIONS = ("S/C", "S/DIS", "C/L", "G/C", "G/L")
EXTENDED_TRANSITIONS = ("FCC/BCC", "FCC/DIS", "C/O70", "G/O70", "O70/L")
FCC_O70_CAMPAIGN_TRANSITIONS = (*EXTENDED_TRANSITIONS, "G/L")
TRANSITIONS = ("S/C", "S/DIS", "G/C", "G/L", *EXTENDED_TRANSITIONS)
TITLE = "BVK2 continuous phase diagram"
LOWER_O70_TRIPLE_FILENAME = "estimated_g_c_o70_triple_point.csv"
UPPER_O70_TRIPLE_FILENAME = "estimated_g_o70_l_triple_point.csv"
RPA_FILENAME = "rpa_stability_limit.csv"
FCC_TRIPLE_FILENAME = "fcc_region_zoom_triple_point.csv"
RPA_SCHEMA = "polyorder-rpa-ab-diblock-stability-v1"

# Campaign acceptance gates.  These MUST match the unit scripts' tolerances
# (run_bvk2_gc_gl_boundary_unit.jl B2B_DF_TOL / B2B_DF_TOL_LOW_CHI): when the
# assembler demanded 0.001 below chiN 15 while the unit converged to 0.002,
# the low-chiN roots were computed and then silently discarded here.
# Composition resolution is 0.002 everywhere (user 2026-07-17: 0.002 is
# acceptable, not a hard limit).  Each row still records its achieved
# bracket_width, so tighter roots remain visible as such.
DF_TOL = 0.002
DF_TOL_LOW_CHI = 0.002
LOW_CHI = 15.0
CHIN_TOL = 0.1
ROOT_RESIDUAL_TOL = 1.0e-3
LAM_CELL_STRESS_TARGET = 1.0e-3
LAM_CELL_STRESS_MAX = 2.0e-3
GYR_CELL_STRESS_TARGET = 1.0e-3
GYR_CELL_STRESS_MAX = 2.0e-3
FILTERED_R6_GC_CHI_SLICES = (
    30.0, 32.0, 34.0, 37.0, 40.0, 42.0, 45.0, 48.0, 50.0,
)
SDIS_FIXED_FA_REFINEMENT_SLICES = (
    0.40, 0.41, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49,
)
SDIS_FIXED_FA_CHIN_TARGET = 0.01
SC_FIXED_FA_CHIN_TARGET = 0.05


def finite_float(value: object) -> float | None:
    try:
        number = float(str(value))
    except (TypeError, ValueError):
        return None
    return number if number == number and abs(number) != float("inf") else None


def truth(value: object) -> bool:
    return str(value).strip().lower() in ("true", "1", "yes", "y")


def read_rows(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def resolve_provenance_path(value: object, ledger: Path) -> Path:
    path = Path(str(value))
    if path.is_absolute():
        return path.resolve()
    if path.parts and path.parts[0] == "results":
        return (PROJECT / path).resolve()
    return (ledger.parent / path).resolve()


def artifact_fingerprint_path(artifact: Path) -> Path:
    if artifact.is_file():
        return artifact
    for name in ("theta.csv", "density.csv"):
        candidate = artifact / name
        if candidate.is_file():
            return candidate
    raise ValueError(f"{artifact}: missing theta/density provenance artifact")


def read_points(directory: Path) -> list[dict[str, str]]:
    points = directory / "points"
    if not points.is_dir():
        return []
    rows: list[dict[str, str]] = []
    for path in sorted(points.glob("*.csv")):
        rows.extend(read_rows(path))
    return rows


def compute_rpa_stability_limit(
    output: Path,
    *,
    project: Path = DEFAULT_POLYORDER_PROJECT,
    script: Path = DEFAULT_RPA_SCRIPT,
) -> None:
    """Generate the analytical AB-diblock spinodal with Polyorder.RPA."""
    subprocess.run(
        [
            "julia",
            "--startup-file=no",
            f"--project={project.resolve()}",
            str(script.resolve()),
            "--output",
            str(output.resolve()),
        ],
        cwd=PROJECT,
        check=True,
    )


def load_rpa_stability_limit(path: Path) -> tuple[tuple[float, float], ...]:
    """Load and validate the Polyorder.RPA stability-limit artifact."""
    points: list[tuple[float, float]] = []
    seen: set[float] = set()
    for row in read_rows(path):
        if row.get("schema") != RPA_SCHEMA:
            raise ValueError(f"{path}: unexpected RPA schema")
        f_a = finite_float(row.get("fA"))
        chi_n = finite_float(row.get("chiN_spinodal"))
        d_star = finite_float(row.get("D_star_Rg"))
        k_star = finite_float(row.get("k_star_Rg_inv"))
        if None in (f_a, chi_n, d_star, k_star):
            raise ValueError(f"{path}: non-finite RPA field")
        if not (0.0 < f_a <= 0.5 and chi_n > 0.0 and d_star > 0.0):
            raise ValueError(f"{path}: non-physical RPA row")
        if not isclose(k_star, 2.0 * pi / d_star, rel_tol=1.0e-10):
            raise ValueError(f"{path}: inconsistent RPA k*/D* pair")
        if f_a in seen:
            raise ValueError(f"{path}: duplicate RPA composition {f_a}")
        seen.add(f_a)
        points.append((f_a, chi_n))
    points.sort()
    if len(points) < 2 or not isclose(points[-1][0], 0.5, abs_tol=1.0e-12):
        raise ValueError(f"{path}: RPA curve must reach f_A=0.5")
    if not isclose(points[-1][1], 10.4949, abs_tol=5.0e-4):
        raise ValueError(f"{path}: invalid symmetric-diblock RPA anchor")
    return tuple(points)


def load_estimated_triple_point(
    path: Path, connected_transitions: tuple[str, ...]
) -> tuple[float, float, tuple[str, ...]]:
    rows = read_rows(path)
    if len(rows) != 1:
        raise RuntimeError(f"expected one estimated triple point in {path}")
    row = rows[0]
    if (
        row.get("claim_kind") != "estimated_triple_point"
        or row.get("status") != "estimated_display_only"
    ):
        raise RuntimeError(f"invalid estimated triple-point contract in {path}")
    f_a = finite_float(row.get("fA"))
    chi_n = finite_float(row.get("chiN"))
    if f_a is None or chi_n is None:
        raise RuntimeError(f"non-finite estimated triple point in {path}")
    return f_a, chi_n, connected_transitions


def load_fcc_o70_boundaries(path: Path) -> list[dict[str, object]]:
    """Load accepted O70 brackets; FCC comes from the strict matched-dx ledger."""
    accepted: list[dict[str, object]] = []
    for row in read_rows(path):
        transition = str(row.get("boundary", ""))
        f_a = finite_float(row.get("fA_estimate"))
        chi_n = finite_float(row.get("chiN"))
        width = finite_float(row.get("width"))
        if (
            row.get("status") != "accepted"
            or transition not in FCC_O70_CAMPAIGN_TRANSITIONS
            or transition in {"BCC/FCC", "FCC/BCC", "FCC/DIS"}
            or f_a is None
            or chi_n is None
            or width is None
        ):
            continue
        accepted.append({
            "transition": transition,
            "chiN": chi_n,
            "fA": f_a,
            "bracket_width": width,
            "claim_kind": BOUNDARY_BRACKET,
            "status": "accepted",
            "acceptance_basis": "accepted_signed_fcc_o70_bracket",
            "source": str(path),
        })
    return accepted


def load_strict_fcc_boundaries(path: Path) -> list[dict[str, object]]:
    """Load accepted equilibrium FCC edges from the matched-dx ledger."""
    transition_names = {"FCC/DIS": "FCC/DIS", "FCC/BCC": "FCC/BCC"}
    accepted: list[dict[str, object]] = []
    for row in read_rows(path):
        source_transition = str(row.get("boundary", ""))
        f_a = finite_float(row.get("fA_estimate"))
        chi_n = finite_float(row.get("chiN"))
        width = finite_float(row.get("width"))
        if (
            source_transition not in transition_names
            or row.get("status") != "accepted"
            or not truth(row.get("equilibrium_boundary"))
            or not truth(row.get("resolution_uncertainty_pass"))
            or f_a is None
            or chi_n is None
            or width is None
        ):
            continue
        accepted.append({
            "transition": transition_names[source_transition],
            "chiN": chi_n,
            "fA": f_a,
            "bracket_width": width,
            "claim_kind": BOUNDARY_BRACKET,
            "status": "accepted",
            "acceptance_basis": "accepted_strict_matched_dx_equilibrium_boundary",
            "source": str(path),
        })
    return accepted


def load_fcc_triple_point(
    path: Path,
) -> tuple[float, float, tuple[str, ...]]:
    """Load the display-only FCC/BCC/DIS pocket-closure estimate."""
    rows = read_rows(path)
    if len(rows) != 1:
        raise RuntimeError(f"expected one FCC triple-point estimate in {path}")
    row = rows[0]
    if row.get("status") != "estimated_linear_pocket_closure":
        raise RuntimeError(f"invalid FCC triple-point contract in {path}")
    f_a = finite_float(row.get("fA_estimate"))
    chi_n = finite_float(row.get("chiN_estimate"))
    if f_a is None or chi_n is None:
        raise RuntimeError(f"non-finite FCC triple point in {path}")
    return f_a, chi_n, ("S/DIS", "FCC/DIS", "FCC/BCC")


def stable_fcc_plot_rows(
    rows: list[dict[str, object]], triple_chi_n: float
) -> list[dict[str, object]]:
    """Exclude pairwise FCC continuations outside their stable regions."""
    return [
        row for row in rows
        if not (
            row.get("transition") == "S/DIS"
            and float(row["chiN"]) > triple_chi_n
        )
        and not (
            row.get("transition") in {"FCC/DIS", "FCC/BCC"}
            and float(row["chiN"]) < triple_chi_n
        )
    ]


def open_o70_slices(path: Path) -> set[float]:
    """Return slices where signed O70 edges enclose a nonzero interval.

    Provisional edges are sufficient to establish that a pairwise G/L root
    lies inside the O70 pocket, even though they remain ineligible for the
    accepted publication boundary itself.
    """
    edges: dict[str, dict[float, float]] = {"G/O70": {}, "O70/L": {}}
    for row in read_rows(path):
        transition = str(row.get("boundary", ""))
        chi_n = finite_float(row.get("chiN"))
        f_a = finite_float(row.get("fA_estimate"))
        if (
            transition not in edges
            or row.get("status") not in {"accepted", "provisional"}
            or chi_n is None
            or f_a is None
        ):
            continue
        edges[transition][chi_n] = f_a
    return {
        chi_n
        for chi_n in edges["G/O70"].keys() & edges["O70/L"].keys()
        if edges["G/O70"][chi_n] < edges["O70/L"][chi_n]
    }


def filter_metastable_o70_pairwise(
    rows: list[dict[str, object]],
    additional_open_slices: set[float] | None = None,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    """Remove pairwise O70 roots that do not bound a stable O70 interval."""
    co70 = {
        float(row["chiN"]): float(row["fA"])
        for row in rows
        if row.get("transition") == "C/O70"
        and row.get("status") == "accepted"
    }
    go70 = {
        float(row["chiN"]): float(row["fA"])
        for row in rows
        if row.get("transition") == "G/O70"
        and row.get("status") == "accepted"
    }
    o70l = {
        float(row["chiN"]): float(row["fA"])
        for row in rows
        if row.get("transition") == "O70/L"
        and row.get("status") == "accepted"
    }
    closed_o70_slices = {
        chi_n for chi_n in go70.keys() & o70l.keys()
        if go70[chi_n] >= o70l[chi_n]
    }
    open_o70_slices = {
        chi_n for chi_n in go70.keys() & o70l.keys()
        if go70[chi_n] < o70l[chi_n]
    }
    reversed_co70_slices = {
        chi_n for chi_n in co70.keys() & o70l.keys()
        if co70[chi_n] >= o70l[chi_n]
    }
    open_o70_slices.update(additional_open_slices or set())
    metastable = [
        row for row in rows
        if (
            row.get("transition") in {"G/O70", "G/C"}
            and float(row["chiN"]) in co70
        ) or (
            row.get("transition") in {"G/O70", "O70/L"}
            and float(row["chiN"]) in closed_o70_slices
        ) or (
            row.get("transition") in {"C/O70", "O70/L"}
            and float(row["chiN"]) in reversed_co70_slices
        ) or (
            row.get("transition") == "G/L"
            and float(row["chiN"]) in open_o70_slices
        )
    ]
    return [row for row in rows if row not in metastable], metastable


def apply_ud_theta_sc_overrides(
    rows: list[dict[str, object]],
    ledger: Path,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    """Replace same-chiN S/C rows with validated UD-theta roots.

    The dedicated ledger assembler has already checked the signed 4x bracket.
    Recheck its immutable contract here so a hand-edited CSV cannot quietly
    enter the publication curve.
    """
    if not ledger.is_file():
        return rows, []
    overrides: list[dict[str, object]] = []
    for source in read_rows(ledger):
        if source.get("boundary") != "S/C":
            raise ValueError(f"{ledger}: only S/C overrides are permitted")
        if source.get("schema") != (
            "ud-theta-spectral-adaptive-k-oversampled-v1"
        ):
            raise ValueError(f"{ledger}: unexpected UD-theta schema")
        if source.get("status") != CLAIM_ACCEPTED:
            raise ValueError(f"{ledger}: non-accepted override row")
        require_boundary_claim(source)
        chi_n = finite_float(source.get("chiN"))
        f_a = finite_float(source.get("fA_secant"))
        width = finite_float(source.get("bracket_width"))
        delta_lower = finite_float(source.get("deltaF_lower_4x"))
        delta_upper = finite_float(source.get("deltaF_upper_4x"))
        audit_factor = finite_float(source.get("audit_factor"))
        if None in (
            chi_n, f_a, width, delta_lower, delta_upper, audit_factor
        ):
            raise ValueError(f"{ledger}: nonfinite override field")
        if width > DF_TOL + 1.0e-12:
            raise ValueError(f"{ledger}: S/C bracket exceeds {DF_TOL}")
        if delta_lower * delta_upper > 0.0:
            raise ValueError(f"{ledger}: S/C 4x endpoint signs do not bracket")
        if audit_factor < 4:
            raise ValueError(f"{ledger}: S/C audit factor is below 4")
        overrides.append({
            "transition": "S/C",
            "chiN": chi_n,
            "fA": f_a,
            "bracket_width": width,
            "claim_kind": BOUNDARY_BRACKET,
            "status": "accepted",
            "acceptance_basis": (
                "BVK2 Uneyama-Doi theta-before-discretization root; "
                "3x optimization; signed 4x endpoint audit; fixed-chiN "
                "fA bracket"
            ),
            "source": source.get("source", str(ledger.resolve())),
        })
    replaced = {("S/C", float(row["chiN"])) for row in overrides}
    kept = [
        row for row in rows
        if (str(row["transition"]), float(row["chiN"])) not in replaced
    ]
    kept.extend(overrides)
    kept.sort(key=lambda item: (
        str(item["transition"]), float(item["chiN"] or 0.0)
    ))
    return kept, overrides


def apply_ud_theta_order_order_override(
    rows: list[dict[str, object]],
    ledger: Path,
    *,
    source_boundary: str,
    transition: str,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    """Replace one order-order row with a validated UD-theta bracket."""
    if not ledger.is_file():
        return rows, []
    overrides: list[dict[str, object]] = []
    for source in read_rows(ledger):
        normalized = source.get("boundary", "").replace("_", "/")
        expected = source_boundary.replace("_", "/")
        if normalized != expected:
            raise ValueError(
                f"{ledger}: expected boundary {source_boundary}"
            )
        if source.get("schema") != (
            "ud-theta-spectral-adaptive-k-oversampled-v1"
        ):
            raise ValueError(f"{ledger}: unexpected UD-theta schema")
        if source.get("status") != CLAIM_ACCEPTED:
            raise ValueError(f"{ledger}: non-accepted override row")
        require_boundary_claim(source)
        chi_n = finite_float(source.get("chiN"))
        f_a = finite_float(source.get("fA_secant"))
        width = finite_float(source.get("bracket_width"))
        audit_factor = finite_float(source.get("audit_factor"))
        if None in (chi_n, f_a, width, audit_factor):
            raise ValueError(f"{ledger}: nonfinite override field")
        audit = int(audit_factor)
        if audit < 3:
            raise ValueError(f"{ledger}: audit factor is below 3")
        optimizer_factor = finite_float(source.get("optimizer_factor"))
        if optimizer_factor != 2:
            raise ValueError(f"{ledger}: optimizer factor is not 2")
        if source.get("gyr_resolution") != "96x96x96":
            raise ValueError(f"{ledger}: GYR resolution is not 96^3")
        if source.get("cyl_resolution") != "48x84":
            raise ValueError(f"{ledger}: CYL resolution is not 48x84")
        if source.get("energy_contract") != (
            "matched_factor_3_gyr_minus_cyl"
        ):
            raise ValueError(f"{ledger}: unexpected G/C energy contract")
        delta_lower = finite_float(source.get(f"deltaF_lower_{audit}x"))
        delta_upper = finite_float(source.get(f"deltaF_upper_{audit}x"))
        if None in (delta_lower, delta_upper):
            raise ValueError(f"{ledger}: missing signed audit deltas")
        if width > DF_TOL + 1.0e-12:
            raise ValueError(f"{ledger}: bracket exceeds {DF_TOL}")
        if delta_lower * delta_upper > 0.0:
            raise ValueError(f"{ledger}: endpoint signs do not bracket")
        overrides.append({
            "transition": transition,
            "chiN": chi_n,
            "fA": f_a,
            "bracket_width": width,
            "claim_kind": BOUNDARY_BRACKET,
            "status": "accepted",
            "acceptance_basis": (
                "BVK2 Uneyama-Doi theta-before-discretization root; "
                f"signed {audit}x endpoint audit; fixed-chiN fA bracket"
            ),
            "source": str(ledger.resolve()),
        })
    replaced = {
        (transition, float(row["chiN"])) for row in overrides
    }
    kept = [
        row for row in rows
        if (str(row["transition"]), float(row["chiN"])) not in replaced
    ]
    kept.extend(overrides)
    kept.sort(key=lambda item: (
        str(item["transition"]), float(item["chiN"] or 0.0)
    ))
    return kept, overrides


def apply_ud_theta_direct_gc_override(
    rows: list[dict[str, object]],
    ledger: Path,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    """Replace G/C rows with accepted direct GYR112/CYL96x168 brackets."""
    if not ledger.is_file():
        return rows, []
    overrides: list[dict[str, object]] = []
    for source in read_rows(ledger):
        if source.get("boundary", "").replace("_", "/") != "G/C":
            raise ValueError(f"{ledger}: expected G/C boundary")
        if source.get("schema") != (
            "ud-theta-spectral-adaptive-k-oversampled-v1"
        ):
            raise ValueError(f"{ledger}: unexpected direct G/C schema")
        require_boundary_claim(source)
        if source.get("status") != CLAIM_ACCEPTED:
            raise ValueError(f"{ledger}: non-accepted direct G/C row")
        chi_n = finite_float(source.get("chiN"))
        f_a = finite_float(source.get("fA_secant", source.get("fA")))
        width = finite_float(source.get("bracket_width"))
        delta_lower = finite_float(source.get("deltaF_lower_3x"))
        delta_upper = finite_float(source.get("deltaF_upper_3x"))
        optimizer_factor = finite_float(source.get("optimizer_factor"))
        audit_factor = finite_float(source.get("audit_factor"))
        if None in (
            chi_n,
            f_a,
            width,
            delta_lower,
            delta_upper,
            optimizer_factor,
            audit_factor,
        ):
            raise ValueError(f"{ledger}: incomplete direct G/C row")
        if optimizer_factor != 2 or audit_factor != 3:
            raise ValueError(f"{ledger}: expected factor-2/factor-3 contract")
        if source.get("gyr_resolution") != "112x112x112":
            raise ValueError(f"{ledger}: direct G/C GYR grid is not 112^3")
        if source.get("cyl_resolution") != "96x168":
            raise ValueError(f"{ledger}: direct G/C CYL grid is not 96x168")
        if source.get("energy_contract") != (
            "direct_factor_3_gyr112_minus_cyl96x168"
        ):
            raise ValueError(f"{ledger}: unexpected direct G/C energy contract")
        if width > DF_TOL + 1.0e-12:
            raise ValueError(f"{ledger}: direct G/C bracket exceeds {DF_TOL}")
        if delta_lower * delta_upper > 0.0:
            raise ValueError(f"{ledger}: direct G/C endpoints do not bracket")
        overrides.append({
            "transition": "G/C",
            "chiN": chi_n,
            "fA": f_a,
            "bracket_width": width,
            "claim_kind": BOUNDARY_BRACKET,
            "status": CLAIM_ACCEPTED,
            "acceptance_basis": (
                "BVK2 Uneyama-Doi theta-before-discretization; direct "
                "GYR 112^3 and CYL 96x168 factor-3 energies; signed "
                "fixed-chiN fA bracket"
            ),
            "source": str(ledger.resolve()),
        })
    replaced = {("G/C", float(row["chiN"])) for row in overrides}
    kept = [
        row for row in rows
        if (str(row["transition"]), float(row["chiN"])) not in replaced
    ]
    kept.extend(overrides)
    kept.sort(key=lambda item: (
        str(item["transition"]), float(item["chiN"] or 0.0)
    ))
    return kept, overrides


def apply_ud_theta_filtered_r6_gc_override(
    rows: list[dict[str, object]],
    ledger: Path,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    """Replace production G/C slices with the global-filter r=6 ledger.

    These rows are the final G/C authority at chiN=30--50.  The campaign
    compares direct factor-3 energies from GYR96 and CYL128x224 using one
    global adaptive-curvature sensor-filter ratio, and publishes only complete
    signed brackets.  Rechecking the immutable contract here prevents a
    partial campaign or a hand-edited coordinate from entering the figure.
    """
    if not ledger.is_file():
        return rows, []
    overrides: list[dict[str, object]] = []
    for source in read_rows(ledger):
        if source.get("transition") != "G/C":
            raise ValueError(f"{ledger}: only G/C overrides are permitted")
        require_boundary_claim(source)
        if source.get("status") != CLAIM_ACCEPTED:
            raise ValueError(f"{ledger}: non-accepted filtered-r6 G/C row")
        if source.get("status_reason") != "signed_filtered_bracket":
            raise ValueError(f"{ledger}: unexpected filtered-r6 status reason")
        if source.get("energy_contract") != (
            "same_schema_filtered_gyr96_minus_cyl128x224_r6"
        ):
            raise ValueError(f"{ledger}: unexpected filtered-r6 energy contract")

        chi_n = finite_float(source.get("chiN"))
        f_a = finite_float(source.get("coordinate"))
        half_width = finite_float(source.get("coordinate_half_width"))
        lower = finite_float(source.get("coordinate_lower"))
        upper = finite_float(source.get("coordinate_upper"))
        delta_lower = finite_float(source.get("delta_lower"))
        delta_upper = finite_float(source.get("delta_upper"))
        filter_ratio = finite_float(source.get("sensor_filter_ratio"))
        if None in (
            chi_n,
            f_a,
            half_width,
            lower,
            upper,
            delta_lower,
            delta_upper,
            filter_ratio,
        ):
            raise ValueError(f"{ledger}: incomplete filtered-r6 G/C row")
        width = 2.0 * half_width
        interval_width = upper - lower
        if width <= 0.0 or width > DF_TOL + 1.0e-12:
            raise ValueError(f"{ledger}: filtered-r6 bracket exceeds {DF_TOL}")
        if not isclose(interval_width, width, rel_tol=0.0, abs_tol=1.0e-12):
            raise ValueError(f"{ledger}: inconsistent filtered-r6 half width")
        if not lower <= f_a <= upper:
            raise ValueError(f"{ledger}: filtered-r6 root lies outside bracket")
        if delta_lower * delta_upper > 0.0:
            raise ValueError(f"{ledger}: filtered-r6 endpoint signs do not bracket")
        if not isclose(filter_ratio, 6.0, rel_tol=0.0, abs_tol=1.0e-12):
            raise ValueError(f"{ledger}: sensor-filter ratio is not 6")

        overrides.append({
            "transition": "G/C",
            "chiN": chi_n,
            "fA": f_a,
            "bracket_width": width,
            "claim_kind": BOUNDARY_BRACKET,
            "status": CLAIM_ACCEPTED,
            "acceptance_basis": (
                "BVK2 Uneyama-Doi theta-before-discretization; global "
                "adaptive-curvature sensor filter r=6; direct factor-3 "
                "GYR 96^3 and CYL 128x224 energies; signed fixed-chiN "
                "fA bracket"
            ),
            "source": str(ledger.resolve()),
        })

    actual_slices = tuple(sorted(float(row["chiN"]) for row in overrides))
    if actual_slices != FILTERED_R6_GC_CHI_SLICES:
        raise ValueError(
            f"{ledger}: expected complete filtered-r6 G/C slices "
            f"{FILTERED_R6_GC_CHI_SLICES}, got {actual_slices}"
        )
    replaced = {("G/C", float(row["chiN"])) for row in overrides}
    kept = [
        row for row in rows
        if (str(row["transition"]), float(row["chiN"])) not in replaced
    ]
    kept.extend(overrides)
    kept.sort(key=lambda item: (
        str(item["transition"]), float(item["chiN"] or 0.0)
    ))
    return kept, overrides


def apply_ud_theta_corrected_gl_override(
    rows: list[dict[str, object]],
    ledger: Path,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    """Replace G/L rows with validated UD-theta continuum-energy roots."""
    if not ledger.is_file():
        return rows, []
    overrides: list[dict[str, object]] = []
    for source in read_rows(ledger):
        if source.get("boundary", "").replace("_", "/") != "G/L":
            raise ValueError(f"{ledger}: only G/L overrides are permitted")
        if source.get("schema") != (
            "ud-theta-spectral-adaptive-k-oversampled-v1"
        ):
            raise ValueError(f"{ledger}: unexpected UD-theta schema")
        if source.get("status") != CLAIM_ACCEPTED:
            raise ValueError(f"{ledger}: non-accepted override row")
        require_boundary_claim(source)
        chi_n = finite_float(source.get("chiN"))
        f_a = finite_float(source.get("fA_secant"))
        width = finite_float(source.get("bracket_width"))
        delta_lower = finite_float(source.get("deltaF_lower_corrected"))
        delta_upper = finite_float(source.get("deltaF_upper_corrected"))
        order = finite_float(source.get("richardson_order"))
        lam_resolution = finite_float(source.get("lam_resolution"))
        if None in (
            chi_n, f_a, width, delta_lower, delta_upper, order,
            lam_resolution,
        ):
            raise ValueError(f"{ledger}: nonfinite override field")
        if width > DF_TOL + 1.0e-12:
            raise ValueError(f"{ledger}: G/L bracket exceeds {DF_TOL}")
        if delta_lower * delta_upper > 0.0:
            raise ValueError(
                f"{ledger}: corrected endpoint signs do not bracket"
            )
        if source.get("gyr_resolution_pair") != (
            "96x96x96;112x112x112"
        ):
            raise ValueError(f"{ledger}: unexpected GYR resolution pair")
        if int(lam_resolution) != 1024:
            raise ValueError(f"{ledger}: LAM resolution is not 1024")
        if abs(order - 2.0) > 1.0e-12:
            raise ValueError(f"{ledger}: Richardson order is not 2")
        if source.get("energy_contract") != (
            "gyr_richardson_order_2_minus_lam1024"
        ):
            raise ValueError(f"{ledger}: unexpected G/L energy contract")
        overrides.append({
            "transition": "G/L",
            "chiN": chi_n,
            "fA": f_a,
            "bracket_width": width,
            "claim_kind": BOUNDARY_BRACKET,
            "status": "accepted",
            "acceptance_basis": (
                "BVK2 Uneyama-Doi theta-before-discretization root; "
                "GYR 96^3/112^3 second-order Richardson correction; "
                "direct LAM N=1024; signed corrected fixed-chiN fA bracket"
            ),
            "source": str(ledger.resolve()),
        })
    replaced = {("G/L", float(row["chiN"])) for row in overrides}
    kept = [
        row for row in rows
        if (str(row["transition"]), float(row["chiN"])) not in replaced
    ]
    kept.extend(overrides)
    kept.sort(key=lambda item: (
        str(item["transition"]), float(item["chiN"] or 0.0)
    ))
    return kept, overrides


def apply_ud_theta_direct_gl_override(
    rows: list[dict[str, object]],
    ledger: Path,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    """Replace G/L rows with direct GYR112/LAM1024 coordinates."""
    if not ledger.is_file():
        return rows, []
    overrides: list[dict[str, object]] = []
    for source in read_rows(ledger):
        if source.get("boundary", "").replace("_", "/") != "G/L":
            raise ValueError(f"{ledger}: only G/L overrides are permitted")
        if source.get("schema") != (
            "ud-theta-spectral-adaptive-k-oversampled-v1"
        ):
            raise ValueError(f"{ledger}: unexpected UD-theta schema")
        if source.get("gyr_resolution") != "112x112x112":
            raise ValueError(f"{ledger}: GYR resolution is not 112^3")
        if source.get("lam_resolution") != "1024":
            raise ValueError(f"{ledger}: LAM resolution is not 1024")
        if source.get("energy_contract") != "gyr112_minus_lam1024":
            raise ValueError(f"{ledger}: unexpected direct G/L contract")
        source_path = Path(source.get("source", ""))
        if not source_path.is_file():
            raise ValueError(f"{ledger}: missing direct G/L source")
        if source.get("source_sha256") != sha256(source_path):
            raise ValueError(f"{ledger}: direct G/L source hash mismatch")
        optimizer_factor = finite_float(source.get("optimizer_factor"))
        if optimizer_factor not in {2, 3}:
            raise ValueError(f"{ledger}: optimizer factor is not 2 or 3")
        audit = finite_float(source.get("audit_factor"))
        if audit is None or audit < 3:
            raise ValueError(f"{ledger}: audit factor is below 3")

        chi_n = finite_float(source.get("chiN"))
        f_a = finite_float(source.get("fA_secant"))
        f_lower = finite_float(source.get("fA_lower"))
        f_upper = finite_float(source.get("fA_upper"))
        width = finite_float(source.get("bracket_width"))
        delta_lower = finite_float(source.get("deltaF_lower_direct"))
        delta_upper = finite_float(source.get("deltaF_upper_direct"))
        lam_stress_max = finite_float(source.get("lam_cell_stress_max"))
        gyr_stress_max = finite_float(source.get("gyr_cell_stress_max"))
        if None in (
            chi_n,
            f_a,
            f_lower,
            f_upper,
            width,
            delta_lower,
            delta_upper,
            lam_stress_max,
            gyr_stress_max,
        ):
            raise ValueError(f"{ledger}: nonfinite direct G/L field")
        if abs((f_upper - f_lower) - width) > 1.0e-10:
            raise ValueError(f"{ledger}: inconsistent endpoint width")

        status = source.get("status")
        claim_kind = source.get("claim_kind")
        if claim_kind == BOUNDARY_BRACKET and status == CLAIM_ACCEPTED:
            if gyr_stress_max > GYR_CELL_STRESS_MAX:
                raise ValueError(
                    f"{ledger}: accepted direct G/L GYR cell stress "
                    f"{gyr_stress_max:g} exceeds {GYR_CELL_STRESS_MAX:g}"
                )
            if lam_stress_max > LAM_CELL_STRESS_MAX:
                raise ValueError(
                    f"{ledger}: accepted direct G/L LAM cell stress "
                    f"{lam_stress_max:g} exceeds {LAM_CELL_STRESS_MAX:g}"
                )
            if width > DF_TOL + 1.0e-12:
                raise ValueError(
                    f"{ledger}: G/L bracket exceeds {DF_TOL}"
                )
            if delta_lower * delta_upper > 0.0:
                raise ValueError(
                    f"{ledger}: direct endpoint signs do not bracket"
                )
            if not f_lower <= f_a <= f_upper:
                raise ValueError(
                    f"{ledger}: direct secant lies outside bracket"
                )
            stress_tier = (
                "target"
                if (
                    lam_stress_max <= LAM_CELL_STRESS_TARGET
                    and gyr_stress_max <= GYR_CELL_STRESS_TARGET
                )
                else "fallback"
            )
            gate = (
                "signed direct fixed-chiN fA bracket; "
                f"GYR/LAM cell stresses "
                f"{gyr_stress_max:.3g}/{lam_stress_max:.3g} "
                f"({stress_tier}, target {LAM_CELL_STRESS_TARGET:g}, "
                f"maximum {LAM_CELL_STRESS_MAX:g})"
            )
            output_status = "accepted"
        elif claim_kind == BOUNDARY_BRACKET and status == CLAIM_PROVISIONAL:
            gate = (
                f"provisional: {source.get('status_reason', 'incomplete_evidence')}; "
                f"current maximum endpoint stresses "
                f"{gyr_stress_max:.3g}/{lam_stress_max:.3g}; "
                f"target {GYR_CELL_STRESS_TARGET:g}, "
                f"maximum {GYR_CELL_STRESS_MAX:g}"
            )
            output_status = "provisional"
        elif claim_kind == BOUNDARY_BRACKET and status == CLAIM_REJECTED:
            gate = (
                f"rejected: {source.get('status_reason', 'invalid_boundary')}"
            )
            output_status = "rejected"
        else:
            raise ValueError(
                f"{ledger}: noncanonical direct G/L claim "
                f"({claim_kind=}, {status=})"
            )

        overrides.append({
            "transition": "G/L",
            "chiN": chi_n,
            "fA": f_a,
            "bracket_width": width,
            "claim_kind": BOUNDARY_BRACKET,
            "status": output_status,
            "acceptance_basis": (
                "BVK2 Uneyama-Doi theta-before-discretization; "
                "direct GYR 112^3 and LAM N=1024 "
                f"optimizer-factor-{optimizer_factor:g}, factor-3 energies; "
                f"{gate}"
            ),
            "source": source.get("source", str(ledger.resolve())),
        })
    replaced = {("G/L", float(row["chiN"])) for row in overrides}
    kept = [
        row for row in rows
        if (str(row["transition"]), float(row["chiN"])) not in replaced
    ]
    kept.extend(overrides)
    kept.sort(key=lambda item: (
        str(item["transition"]), float(item["chiN"] or 0.0)
    ))
    return kept, overrides


def resolution_axis(value: object, *, phase: str) -> int:
    text = str(value).strip().lower().replace("^3", "")
    parts = text.split("x")
    try:
        axes = [int(part) for part in parts]
    except ValueError as error:
        raise ValueError(f"invalid {phase} resolution {value!r}") from error
    if not axes or any(axis <= 0 for axis in axes):
        raise ValueError(f"invalid {phase} resolution {value!r}")
    if phase == "GYR" and any(axis != axes[0] for axis in axes):
        raise ValueError(f"non-cubic GYR resolution {value!r}")
    if phase == "LAM" and len(axes) != 1:
        raise ValueError(f"non-1D LAM resolution {value!r}")
    return axes[0]


def validate_matched_gl_endpoint(
    endpoint: Path,
    expected_sha256: str,
    *,
    ledger: Path,
    chi_n: float,
    f_a: float,
    delta_f: float,
) -> dict[str, object]:
    if not endpoint.is_file():
        raise ValueError(f"{ledger}: missing matched G/L endpoint {endpoint}")
    if sha256(endpoint) != expected_sha256:
        raise ValueError(f"{ledger}: matched G/L endpoint hash mismatch")
    records = read_rows(endpoint)
    if len(records) != 1:
        raise ValueError(f"{endpoint}: expected exactly one endpoint record")
    source = records[0]
    if source.get("claim_kind") != "phase_endpoint":
        raise ValueError(f"{endpoint}: expected phase_endpoint claim")
    if source.get("status") != CLAIM_ACCEPTED:
        raise ValueError(f"{endpoint}: endpoint is not accepted")
    if source.get("schema") != (
        "ud-theta-spectral-adaptive-k-oversampled-v1"
    ):
        raise ValueError(f"{endpoint}: unexpected endpoint schema")
    if source.get("status_reason") != "endpoint_gates_pass":
        raise ValueError(f"{endpoint}: endpoint gates did not pass")
    if not truth(source.get("nearest_even_optimal")):
        raise ValueError(f"{endpoint}: matched LAM grid is not nearest-even")

    fields = {
        name: finite_float(source.get(name))
        for name in (
            "chiN", "fA", "gyr_energy3", "gyr_length", "gyr_n",
            "gyr_dx", "gyr_stress", "matched_lam_energy3",
            "matched_lam_length", "matched_lam_n", "matched_lam_dx",
            "spacing_mismatch", "matched_lam_stress", "matched_lam_r2",
            "matched_lam_rinf", "deltaF",
        )
    }
    if any(value is None for value in fields.values()):
        raise ValueError(f"{endpoint}: nonfinite phase endpoint field")
    values = {name: float(value) for name, value in fields.items()}
    if not isclose(values["chiN"], chi_n, abs_tol=1.0e-10):
        raise ValueError(f"{endpoint}: endpoint chiN does not match bracket")
    if not isclose(values["fA"], f_a, abs_tol=1.0e-10):
        raise ValueError(f"{endpoint}: endpoint fA does not match bracket")
    if not isclose(values["deltaF"], delta_f, abs_tol=1.0e-10):
        raise ValueError(f"{endpoint}: endpoint deltaF does not match bracket")
    if not isclose(
        values["gyr_energy3"] - values["matched_lam_energy3"],
        values["deltaF"],
        abs_tol=1.0e-10,
    ):
        raise ValueError(f"{endpoint}: endpoint energies do not reproduce deltaF")

    gyr_n = int(values["gyr_n"])
    lam_n = int(values["matched_lam_n"])
    if gyr_n != values["gyr_n"] or lam_n != values["matched_lam_n"]:
        raise ValueError(f"{endpoint}: nonintegral phase grid")
    if gyr_n <= 0 or lam_n < 8 or lam_n % 2:
        raise ValueError(f"{endpoint}: invalid matched phase grid")
    computed_gyr_dx = values["gyr_length"] / gyr_n
    computed_lam_dx = values["matched_lam_length"] / lam_n
    computed_mismatch = abs(computed_gyr_dx - computed_lam_dx) / (
        (computed_gyr_dx + computed_lam_dx) / 2.0
    )
    if not isclose(values["gyr_dx"], computed_gyr_dx, abs_tol=1.0e-10):
        raise ValueError(f"{endpoint}: inconsistent GYR spacing")
    if not isclose(values["matched_lam_dx"], computed_lam_dx, abs_tol=1.0e-10):
        raise ValueError(f"{endpoint}: inconsistent matched LAM spacing")
    if not isclose(
        values["spacing_mismatch"], computed_mismatch, abs_tol=1.0e-10
    ):
        raise ValueError(f"{endpoint}: inconsistent spacing mismatch")
    if computed_mismatch > 0.03 + 1.0e-12:
        raise ValueError(f"{endpoint}: physical spacing mismatch exceeds 3%")
    if abs(values["gyr_stress"]) > GYR_CELL_STRESS_MAX:
        raise ValueError(f"{endpoint}: GYR cell stress exceeds endpoint gate")
    if abs(values["matched_lam_stress"]) > LAM_CELL_STRESS_MAX:
        raise ValueError(f"{endpoint}: LAM cell stress exceeds endpoint gate")
    if values["matched_lam_r2"] > LAM_CELL_STRESS_MAX:
        raise ValueError(f"{endpoint}: LAM R2 exceeds endpoint gate")
    if values["matched_lam_rinf"] > LAM_CELL_STRESS_MAX:
        raise ValueError(f"{endpoint}: LAM Rinf exceeds endpoint gate")

    gyr_evidence = resolve_provenance_path(source.get("gyr_evidence", ""), endpoint)
    gyr_artifact = resolve_provenance_path(source.get("gyr_artifact", ""), endpoint)
    lam_result = resolve_provenance_path(source.get("matched_lam_result", ""), endpoint)
    lam_artifact = resolve_provenance_path(source.get("matched_lam_artifact", ""), endpoint)
    if not gyr_evidence.is_file() or not gyr_artifact.exists():
        raise ValueError(f"{endpoint}: missing GYR provenance")
    if not lam_result.is_file() or not lam_artifact.exists():
        raise ValueError(f"{endpoint}: missing matched LAM provenance")
    gyr_fingerprint = source.get("gyr_fingerprint", "")
    if not gyr_fingerprint or sha256(
        artifact_fingerprint_path(gyr_artifact)
    ) != gyr_fingerprint:
        raise ValueError(f"{endpoint}: GYR artifact fingerprint mismatch")
    gyr_rows = read_rows(gyr_evidence)
    if len(gyr_rows) != 1:
        raise ValueError(f"{endpoint}: invalid GYR evidence record")
    gyr_source = gyr_rows[0]
    if gyr_source.get("status") and not gyr_source["status"].startswith(
        "accepted"
    ):
        raise ValueError(f"{endpoint}: GYR evidence is not accepted")
    if gyr_source.get("phase") not in (None, "", "GYR"):
        raise ValueError(f"{endpoint}: GYR evidence phase mismatch")
    for names, expected in (
        (("chiN",), values["chiN"]),
        (("fA",), values["fA"]),
        (("length", "lengths"), values["gyr_length"]),
        (("cell_stress",), values["gyr_stress"]),
        (("energy_density_factor3",), values["gyr_energy3"]),
    ):
        present = next((name for name in names if gyr_source.get(name)), None)
        if present is None:
            continue
        text = gyr_source[present].split(";")[0]
        actual = finite_float(text)
        if present == "cell_stress" and actual is not None:
            actual = abs(actual)
            expected = abs(expected)
        tolerance = 1.0e-9 if present in ("cell_stress", "energy_density_factor3") else 1.0e-10
        if actual is None or not isclose(actual, expected, abs_tol=tolerance):
            raise ValueError(f"{endpoint}: GYR evidence {present} mismatch")

    lam_rows = read_rows(lam_result)
    if len(lam_rows) != 1 or not lam_rows[0].get("status", "").startswith(
        "accepted"
    ):
        raise ValueError(f"{endpoint}: matched LAM result is not accepted")
    lam_source = lam_rows[0]
    result_artifact = resolve_provenance_path(
        lam_source.get("artifact", ""), lam_result
    )
    if result_artifact != lam_artifact:
        raise ValueError(f"{endpoint}: matched LAM artifact provenance mismatch")
    result_checks = (
        ("chiN", values["chiN"]),
        ("fA", values["fA"]),
        ("nx", values["matched_lam_n"]),
        ("length", values["matched_lam_length"]),
        ("cell_stress", values["matched_lam_stress"]),
        ("energy_density_factor3", values["matched_lam_energy3"]),
    )
    for name, expected in result_checks:
        actual = finite_float(lam_source.get(name))
        if actual is None or not isclose(actual, expected, abs_tol=1.0e-10):
            raise ValueError(f"{endpoint}: matched LAM result {name} mismatch")
    lam_fingerprint = source.get("matched_lam_fingerprint", "")
    if lam_fingerprint and sha256(
        artifact_fingerprint_path(lam_artifact)
    ) != lam_fingerprint:
        raise ValueError(f"{endpoint}: matched LAM artifact fingerprint mismatch")

    return {
        "gyr_n": gyr_n,
        "lam_n": lam_n,
        "spacing_mismatch": computed_mismatch,
        "gyr_stress": abs(values["gyr_stress"]),
        "lam_stress": abs(values["matched_lam_stress"]),
    }


def apply_ud_theta_matched_gl_override(
    rows: list[dict[str, object]],
    ledger: Path,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    """Promote rigorously validated matched-spacing GYR/LAM brackets."""
    if not ledger.is_file():
        return rows, []
    overrides: list[dict[str, object]] = []
    for source in read_rows(ledger):
        if source.get("boundary", "").replace("_", "/") != "G/L":
            raise ValueError(f"{ledger}: only G/L promotions are permitted")
        if source.get("schema") != (
            "ud-theta-spectral-adaptive-k-oversampled-v1"
        ):
            raise ValueError(f"{ledger}: unexpected matched G/L schema")
        if source.get("claim_kind") != BOUNDARY_BRACKET:
            raise ValueError(f"{ledger}: expected boundary_bracket claim")
        if source.get("status") != CLAIM_ACCEPTED:
            raise ValueError(f"{ledger}: matched G/L row is not accepted")
        if source.get("status_reason") != "signed_matched_spacing_bracket":
            raise ValueError(f"{ledger}: matched G/L bracket gates did not pass")
        require_boundary_claim(source)
        if source.get("energy_contract") != (
            "gyr_factor3_minus_matched_lam_factor3"
        ):
            raise ValueError(f"{ledger}: unexpected matched G/L energy contract")

        numeric = {
            name: finite_float(source.get(name))
            for name in (
                "chiN", "fA_secant", "fA_lower", "fA_upper",
                "bracket_width", "deltaF_lower", "deltaF_upper",
                "max_spacing_mismatch", "gyr_cell_stress_max",
                "lam_cell_stress_max",
            )
        }
        if any(value is None for value in numeric.values()):
            raise ValueError(f"{ledger}: nonfinite matched G/L field")
        values = {name: float(value) for name, value in numeric.items()}
        if any(values[name] < 0.0 for name in (
            "bracket_width", "max_spacing_mismatch",
            "gyr_cell_stress_max", "lam_cell_stress_max",
        )):
            raise ValueError(f"{ledger}: negative matched G/L gate metric")
        if values["fA_lower"] >= values["fA_upper"]:
            raise ValueError(f"{ledger}: invalid matched G/L bracket order")
        if not isclose(
            values["fA_upper"] - values["fA_lower"],
            values["bracket_width"],
            abs_tol=1.0e-10,
        ):
            raise ValueError(f"{ledger}: inconsistent matched G/L width")
        if values["bracket_width"] > DF_TOL + 1.0e-12:
            raise ValueError(f"{ledger}: matched G/L bracket exceeds {DF_TOL}")
        if values["deltaF_lower"] * values["deltaF_upper"] > 0.0:
            raise ValueError(f"{ledger}: matched G/L endpoints do not bracket")
        if not values["fA_lower"] <= values["fA_secant"] <= values["fA_upper"]:
            raise ValueError(f"{ledger}: matched G/L secant lies outside bracket")
        if values["max_spacing_mismatch"] > 0.03 + 1.0e-12:
            raise ValueError(f"{ledger}: matched G/L spacing mismatch exceeds 3%")
        if values["gyr_cell_stress_max"] > GYR_CELL_STRESS_MAX:
            raise ValueError(f"{ledger}: matched G/L GYR stress exceeds gate")
        if values["lam_cell_stress_max"] > LAM_CELL_STRESS_MAX:
            raise ValueError(f"{ledger}: matched G/L LAM stress exceeds gate")

        lower_path = resolve_provenance_path(source.get("lower_endpoint", ""), ledger)
        upper_path = resolve_provenance_path(source.get("upper_endpoint", ""), ledger)
        lower = validate_matched_gl_endpoint(
            lower_path,
            source.get("lower_endpoint_sha256", ""),
            ledger=ledger,
            chi_n=values["chiN"],
            f_a=values["fA_lower"],
            delta_f=values["deltaF_lower"],
        )
        upper = validate_matched_gl_endpoint(
            upper_path,
            source.get("upper_endpoint_sha256", ""),
            ledger=ledger,
            chi_n=values["chiN"],
            f_a=values["fA_upper"],
            delta_f=values["deltaF_upper"],
        )
        gyr_resolution = resolution_axis(source.get("gyr_resolution"), phase="GYR")
        lam_resolution = resolution_axis(source.get("lam_resolution"), phase="LAM")
        if lower["gyr_n"] != upper["gyr_n"] or lower["gyr_n"] != gyr_resolution:
            raise ValueError(f"{ledger}: GYR grid changes across matched bracket")
        if lower["lam_n"] != upper["lam_n"] or lower["lam_n"] != lam_resolution:
            raise ValueError(f"{ledger}: matched LAM grid changes across bracket")
        if not isclose(
            values["max_spacing_mismatch"],
            max(lower["spacing_mismatch"], upper["spacing_mismatch"]),
            abs_tol=1.0e-10,
        ):
            raise ValueError(f"{ledger}: endpoint spacing maximum mismatch")
        if not isclose(
            values["gyr_cell_stress_max"],
            max(lower["gyr_stress"], upper["gyr_stress"]),
            abs_tol=1.0e-10,
        ):
            raise ValueError(f"{ledger}: endpoint GYR stress maximum mismatch")
        if not isclose(
            values["lam_cell_stress_max"],
            max(lower["lam_stress"], upper["lam_stress"]),
            abs_tol=1.0e-10,
        ):
            raise ValueError(f"{ledger}: endpoint LAM stress maximum mismatch")

        provenance = resolve_provenance_path(source.get("source", ""), ledger)
        if not provenance.is_file():
            raise ValueError(f"{ledger}: missing matched G/L source provenance")
        source_hash = source.get("source_sha256", "")
        if source_hash and sha256(provenance) != source_hash:
            raise ValueError(f"{ledger}: matched G/L source hash mismatch")
        overrides.append({
            "transition": "G/L",
            "chiN": values["chiN"],
            "fA": values["fA_secant"],
            "bracket_width": values["bracket_width"],
            "claim_kind": BOUNDARY_BRACKET,
            "status": CLAIM_ACCEPTED,
            "acceptance_basis": (
                "BVK2 Uneyama-Doi theta-before-discretization; signed "
                f"factor-3 matched-spacing GYR {gyr_resolution}^3 and "
                f"LAM N={lam_resolution} bracket; maximum spacing mismatch "
                f"{values['max_spacing_mismatch']:.3%}"
            ),
            "source": str(ledger.resolve()),
        })

    replaced = {("G/L", float(row["chiN"])) for row in overrides}
    kept = [
        row for row in rows
        if (str(row["transition"]), float(row["chiN"])) not in replaced
    ]
    kept.extend(overrides)
    kept.sort(key=lambda item: (
        str(item["transition"]), float(item["chiN"] or 0.0)
    ))
    return kept, overrides


def apply_o70_chi15_lam1024_gl_override(
    rows: list[dict[str, object]],
    audit: Path,
    boundary_ledger: Path,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    """Restore the stable chiN=15 G/L root under the LAM1024 contract.

    Matched-spacing GYR/LAM promotions are normally authoritative, but the
    O70 campaign uses LAM1024 for every O70/L slice.  Once its two chiN=15
    O70 edges cross, the direct GYR64/LAM1024 bracket is the thermodynamic
    G/L boundary and must supersede the matched-LAM26 diagnostic.
    """
    if not audit.is_file() and not boundary_ledger.is_file():
        return rows, []
    if not audit.is_file():
        raise ValueError(f"{audit}: missing chiN=15 O70 closure audit")
    if not boundary_ledger.is_file():
        raise ValueError(f"{boundary_ledger}: missing FCC/O70 boundary ledger")

    audit_rows = read_rows(audit)
    if len(audit_rows) != 1:
        raise ValueError(f"{audit}: expected exactly one thermodynamic audit row")
    closure = audit_rows[0]
    chi_n = finite_float(closure.get("chiN"))
    direct_root = finite_float(closure.get("direct_g_l_root"))
    g_o70_root = finite_float(closure.get("g_o70_root"))
    o70_l_root = finite_float(closure.get("o70_l_root"))
    pocket_width = finite_float(closure.get("signed_o70_pocket_width"))
    if (
        closure.get("claim_kind") != "thermodynamic_ordering"
        or closure.get("status") != CLAIM_ACCEPTED
        or chi_n is None
        or not isclose(chi_n, 15.0, abs_tol=1.0e-12)
        or None in (direct_root, g_o70_root, o70_l_root, pocket_width)
    ):
        raise ValueError(f"{audit}: invalid chiN=15 thermodynamic audit")
    if truth(closure.get("o70_stable_at_chi15")):
        raise ValueError(f"{audit}: O70 is still stable at chiN=15")
    if not truth(closure.get("g_l_between_metastable_o70_edges")):
        raise ValueError(f"{audit}: direct G/L root is not between O70 edges")
    if not (
        pocket_width < 0.0
        and o70_l_root < direct_root < g_o70_root
        and isclose(
            g_o70_root - o70_l_root,
            -pocket_width,
            abs_tol=1.0e-12,
        )
    ):
        raise ValueError(f"{audit}: inconsistent crossed O70-edge ordering")

    candidates = []
    for source in read_rows(boundary_ledger):
        source_chi_n = finite_float(source.get("chiN"))
        if (
            source.get("boundary", "").replace("_", "/") == "G/L"
            and source_chi_n is not None
            and isclose(source_chi_n, 15.0, abs_tol=1.0e-12)
            and source.get("claim_kind") == BOUNDARY_BRACKET
            and source.get("status") == CLAIM_ACCEPTED
        ):
            candidates.append(source)
    if len(candidates) != 1:
        raise ValueError(
            f"{boundary_ledger}: expected one accepted chiN=15 G/L bracket"
        )
    source = candidates[0]
    f_a = finite_float(source.get("fA_estimate"))
    width = finite_float(source.get("width"))
    lower = finite_float(source.get("fA_lower"))
    upper = finite_float(source.get("fA_upper"))
    delta_lower = finite_float(source.get("delta_lower"))
    delta_upper = finite_float(source.get("delta_upper"))
    if None in (f_a, width, lower, upper, delta_lower, delta_upper):
        raise ValueError(f"{boundary_ledger}: nonfinite chiN=15 G/L bracket")
    if (
        not isclose(f_a, direct_root, abs_tol=1.0e-12)
        or lower >= upper
        or not isclose(upper - lower, width, abs_tol=1.0e-12)
        or width > DF_TOL + 1.0e-12
        or delta_lower * delta_upper > 0.0
        or not lower <= f_a <= upper
    ):
        raise ValueError(f"{boundary_ledger}: invalid chiN=15 G/L bracket")

    override = {
        "transition": "G/L",
        "chiN": 15.0,
        "fA": f_a,
        "bracket_width": width,
        "claim_kind": BOUNDARY_BRACKET,
        "status": CLAIM_ACCEPTED,
        "acceptance_basis": (
            "accepted direct GYR64/LAM1024 signed bracket; final chiN=15 "
            "authority after LAM1024-consistent O70 closure audit"
        ),
        "source": str(boundary_ledger.resolve()),
    }
    kept = [
        row for row in rows
        if not (
            str(row["transition"]) == "G/L"
            and isclose(float(row["chiN"]), 15.0, abs_tol=1.0e-12)
        )
    ]
    kept.append(override)
    kept.sort(key=lambda item: (
        str(item["transition"]), float(item["chiN"] or 0.0)
    ))
    return kept, [override]


def apply_ud_theta_lowchi_gl_status(
    rows: list[dict[str, object]],
    root: Path,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    """Apply direct low-chi G/L verification status to legacy locators."""
    if not root.is_dir():
        return rows, []
    audits: dict[float, tuple[Path, dict[str, str]]] = {}
    canonical = root / "boundary_candidates.csv"
    paths = [canonical] if canonical.is_file() else sorted(
        root.glob("chi*/boundary_point.csv")
    )
    for path in paths:
        source_rows = read_rows(path)
        if path == canonical:
            rows_to_read = source_rows
        else:
            if len(source_rows) != 1:
                raise ValueError(f"{path}: expected one low-chi boundary row")
            rows_to_read = source_rows
        for source in rows_to_read:
            if source.get("boundary", "").replace("_", "/") != "G/L":
                raise ValueError(f"{path}: expected G/L boundary")
            chi_n = finite_float(source.get("chiN"))
            if chi_n is None:
                raise ValueError(f"{path}: nonfinite chiN")
            if path == canonical:
                require_boundary_claim(source)
            elif source.get("schema") != "bvk2-ud-theta-lowchi-direct-v1":
                raise ValueError(f"{path}: unexpected low-chi schema")
            audits[chi_n] = (path, source)

    updated: list[dict[str, object]] = []
    overrides: list[dict[str, object]] = []
    matched: set[float] = set()
    for row in rows:
        key = float(row["chiN"])
        if row["transition"] != "G/L" or key not in audits:
            updated.append(row)
            continue
        matched.add(key)
        path, source = audits[key]
        direct_status = source.get("status", "")
        new = dict(row)
        if direct_status == CLAIM_ACCEPTED:
            f_a = finite_float(
                source.get("fA_secant", source.get("fA"))
            )
            width = finite_float(source.get("bracket_width"))
            delta_lower = finite_float(source.get("deltaF_lower"))
            delta_upper = finite_float(source.get("deltaF_upper"))
            if None in (f_a, width, delta_lower, delta_upper):
                raise ValueError(f"{path}: incomplete accepted bracket")
            if width > DF_TOL + 1.0e-12:
                raise ValueError(f"{path}: accepted bracket is too wide")
            if delta_lower * delta_upper > 0.0:
                raise ValueError(f"{path}: accepted endpoints do not bracket")
            new["fA"] = f_a
            new["bracket_width"] = width
            new["status"] = "accepted"
            new["acceptance_basis"] = (
                "BVK2 Uneyama-Doi theta-before-discretization; direct "
                "low-chi GYR/LAM factor-3 signed bracket"
            )
        elif direct_status == CLAIM_PROVISIONAL:
            f_a = finite_float(
                source.get("fA_secant", source.get("fA"))
            )
            width = finite_float(source.get("bracket_width"))
            if f_a is not None:
                new["fA"] = f_a
            if width is not None:
                new["bracket_width"] = width
            new["status"] = "provisional"
            new["acceptance_basis"] = (
                f"{row['acceptance_basis']}; direct low-chi theta claim "
                f"is provisional ({source.get('status_reason', '')})"
            )
        elif direct_status == CLAIM_REJECTED:
            new["status"] = "rejected"
            new["acceptance_basis"] = (
                f"rejected direct low-chi theta claim "
                f"({source.get('status_reason', '')})"
            )
        else:
            raise ValueError(
                f"{path}: unsupported low-chi status {direct_status!r}"
            )
        new["source"] = str(path.resolve())
        updated.append(new)
        overrides.append(new)

    # A legacy locator may have been rejected before this override stage and
    # therefore be absent from ``rows``.  A complete accepted direct bracket is
    # independent publication evidence and must be able to add that coordinate
    # back instead of requiring a surviving legacy row to overwrite.
    for key, (path, source) in audits.items():
        if key in matched or source.get("status", "") != CLAIM_ACCEPTED:
            continue
        f_a = finite_float(source.get("fA_secant", source.get("fA")))
        width = finite_float(source.get("bracket_width"))
        delta_lower = finite_float(source.get("deltaF_lower"))
        delta_upper = finite_float(source.get("deltaF_upper"))
        if None in (f_a, width, delta_lower, delta_upper):
            raise ValueError(f"{path}: incomplete accepted bracket")
        if width > DF_TOL + 1.0e-12:
            raise ValueError(f"{path}: accepted bracket is too wide")
        if delta_lower * delta_upper > 0.0:
            raise ValueError(f"{path}: accepted endpoints do not bracket")
        new = {
            "transition": "G/L",
            "chiN": key,
            "fA": f_a,
            "bracket_width": width,
            "claim_kind": BOUNDARY_BRACKET,
            "status": "accepted",
            "acceptance_basis": (
                "BVK2 Uneyama-Doi theta-before-discretization; direct "
                "low-chi GYR/LAM factor-3 signed bracket"
            ),
            "source": str(path.resolve()),
        }
        updated.append(new)
        overrides.append(new)

    updated.sort(key=lambda item: (
        str(item["transition"]), float(item["chiN"] or 0.0)
    ))
    return updated, overrides


def load_corrections(root: Path) -> dict[tuple[str, float], dict[str, str]]:
    """Richardson grid corrections, keyed by (transition, chiN).

    Written per-group under `<root>/g*/corrected_roots.csv` (separate dirs so
    concurrent writers cannot interleave appends).
    """
    found: dict[tuple[str, float], dict[str, str]] = {}
    if not root.is_dir():
        return found
    for path in sorted(root.glob("**/corrected_roots.csv")):
        for row in read_rows(path):
            chi = finite_float(row.get("chiN"))
            shift = finite_float(row.get("root_fA_corrected"))
            if chi is None or shift is None:
                continue
            found[(str(row.get("transition")), chi)] = row
    return found


def apply_corrections(
    rows: list[dict[str, object]],
    corrections: dict[tuple[str, float], dict[str, str]],
) -> tuple[list[dict[str, object]], dict[str, str]]:
    """Apply Richardson corrections PER POINT, split by display tier.

    A correction moves a root by up to ~0.04 -- far beyond the 0.002 gate --
    so corrected and raw points must never be CONNECTED on one curve.  But
    they can coexist on the plot: corrected points keep status "accepted"
    (solid curve, filled markers) while not-yet-corrected points are
    demoted to status "provisional", which the shared renderer draws as its
    own dashed, disconnected curve with distinct markers.  Corrected points
    therefore appear the moment their correction lands, raw points remain
    visible (honestly styled as provisional), and no artificial kink can
    arise because tiers are never joined.  Transitions whose corrections
    are all negligible (max |shift| within ~the f_A gate) render raw at
    FULL range as plain accepted points (S/C is the standing example).
    """
    by_transition: dict[str, list[dict[str, object]]] = {}
    for row in rows:
        by_transition.setdefault(str(row["transition"]), []).append(row)

    report: dict[str, str] = {}
    out: list[dict[str, object]] = []
    for transition, group in by_transition.items():
        have = {
            float(r["chiN"]) for r in group
            if (transition, float(r["chiN"])) in corrections
        }
        if not have:
            report[transition] = "raw (no corrections)"
            out.extend(group)
            continue
        max_shift = max(
            abs(finite_float(corrections[(transition, chi)].get("shift"))
                or 0.0)
            for chi in have
        )
        # 1.25x the 0.002 gate: hysteresis so a shift hovering at the gate
        # (S/C peaks at ~0.0020) cannot flip the policy between renders.
        if max_shift <= 0.0025:
            report[transition] = (
                f"raw (corrections negligible: max |shift| "
                f"{max_shift:.5f} <= 0.0025)"
            )
            out.extend(group)
            continue
        n_corrected = 0
        n_raw = 0
        for row in group:
            key = (transition, float(row["chiN"]))
            if key in corrections:
                corr = corrections[key]
                shift = finite_float(corr.get("shift")) or 0.0
                new = dict(row)
                new["fA"] = finite_float(corr.get("root_fA_corrected"))
                new["acceptance_basis"] = (
                    f"{row['acceptance_basis']}; Richardson grid-corrected "
                    f"(shift {shift:+.5f}; both competitors extrapolated to "
                    f"dx->0 via E=E_inf-C*dx^2)"
                )
                out.append(new)
                n_corrected += 1
            else:
                new = dict(row)
                new["status"] = "provisional"
                new["acceptance_basis"] = (
                    f"{row['acceptance_basis']}; RAW grid root, Richardson "
                    "correction pending -- rendered as a provisional "
                    "(dashed, disconnected) tier"
                )
                out.append(new)
                n_raw += 1
        report[transition] = (
            f"CORRECTED {n_corrected}/{len(group)}"
            + (f"; {n_raw} raw point(s) shown as provisional tier"
               if n_raw else "")
        )
    return out, report


def gate_pass(row: dict[str, str]) -> tuple[bool, str]:
    """Re-screen a unit's accepted row against the campaign gates."""
    if str(row.get("status")) != "accepted":
        return False, f"status={row.get('status')}"
    f_a = finite_float(row.get("root_fA"))
    chi_n = finite_float(row.get("root_chiN"))
    if f_a is None or chi_n is None:
        return False, "nonfinite_root"
    if not (0.0 < f_a < 1.0) or not (8.0 <= chi_n <= 64.0):
        return False, "root_outside_diagram_window"
    direction = str(row.get("direction", ""))
    width = finite_float(row.get("bracket_width"))
    if width is None:
        return False, "nonfinite_bracket"
    if direction == "fA":
        limit = DF_TOL_LOW_CHI if chi_n < LOW_CHI else DF_TOL
        if abs(width) > limit + 1.0e-12:
            return False, f"bracket_width={width:g}>{limit:g}"
    elif direction == "chiN":
        if abs(width) > CHIN_TOL + 1.0e-12:
            return False, f"bracket_width={width:g}>{CHIN_TOL:g}"
    else:
        return False, f"unknown_direction={direction}"
    # Order-order units carry a free-energy root residual; the ODT is a
    # bifurcation edge and reports a predicate bracket instead.
    if "final_delta" in row:
        residual = finite_float(row.get("final_delta"))
        if residual is None or abs(residual) > ROOT_RESIDUAL_TOL:
            return False, f"root_residual={row.get('final_delta')}"
    if "phases_accepted" in row and not truth(row.get("phases_accepted")):
        return False, "phase_gate_failed"
    if "ordered_endpoint_accepted" in row and not truth(
        row.get("ordered_endpoint_accepted")
    ):
        return False, "ordered_endpoint_gate_failed"
    return True, "accepted"


def collect(campaign: Path, sdis: Path) -> tuple[
    list[dict[str, object]], list[dict[str, object]],
    list[dict[str, object]],
]:
    accepted: list[dict[str, object]] = []
    rejected: list[dict[str, object]] = []
    diagnostics: list[dict[str, object]] = []
    for directory in (campaign, sdis):
        for row in read_points(directory):
            transition = str(row.get("transition", ""))
            if (
                transition not in PAIRWISE_TRANSITIONS
                and transition != "L/DIS"
            ):
                continue
            ok, reason = gate_pass(row)
            direction = str(row.get("direction", ""))
            delta = row.get("final_delta", "")
            basis = (
                f"BVK2 {direction}-direction root; "
                f"stress-free {row.get('gyr_n', '')}^3-GYR/"
                f"{row.get('bcc_n', '')}^3-BCC campaign grids; "
                f"c2={row.get('c2', '')} frozen; "
                f"reference prior={Path(str(row.get('reference_source', ''))).name}"
                f" (covered={row.get('reference_covered', '')})"
                if ok
                else f"REJECTED: {reason}"
            )
            record = {
                "transition": transition,
                "chiN": finite_float(row.get("root_chiN")),
                "fA": finite_float(row.get("root_fA")),
                "bracket_width": row.get("bracket_width", ""),
                "claim_kind": BOUNDARY_BRACKET,
                "status": "accepted" if ok else "rejected",
                "acceptance_basis": basis,
                "source": str(directory),
            }
            # Numeric diagnostics live alongside, outside the fixed
            # publication schema shared with the BURP-TI/BVK1 diagrams.
            diagnostic = {
                **record,
                "direction": direction,
                "final_delta": delta,
                "reference_source": row.get("reference_source", ""),
                "gate_reason": reason,
            }
            if ok:
                accepted.append(record)
                diagnostics.append(diagnostic)
            else:
                rejected.append(diagnostic)
    accepted.sort(key=lambda item: (item["transition"], item["chiN"] or 0.0))
    diagnostics.sort(key=lambda item: (item["transition"], item["chiN"] or 0.0))
    return accepted, rejected, diagnostics


def apply_sdis_fixed_fa_chi_refinement(
    rows: list[dict[str, object]],
    refinement: Path,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    """Promote completed high-chiN-resolution BCC/DIS roof slices."""
    if not refinement.is_dir():
        return rows, []
    overrides: list[dict[str, object]] = []
    seen: set[float] = set()
    for source in read_points(refinement):
        transition = str(source.get("transition", ""))
        direction = str(source.get("direction", ""))
        target_f_a = finite_float(source.get("target_fA"))
        root_f_a = finite_float(source.get("root_fA"))
        root_chi_n = finite_float(source.get("root_chiN"))
        lower = finite_float(source.get("bracket_lower"))
        upper = finite_float(source.get("bracket_upper"))
        width = finite_float(source.get("bracket_width"))
        resolution = finite_float(source.get("resolution_target"))
        lower_energy = finite_float(source.get("e_ordered_lower"))
        upper_energy = finite_float(source.get("e_ordered_upper"))
        ordered_tolerance = finite_float(source.get("ordered_energy_tolerance"))
        if transition != "S/DIS" or direction != "chiN":
            raise ValueError(f"{refinement}: unexpected refinement claim")
        if None in (
            target_f_a, root_f_a, root_chi_n, lower, upper, width,
            resolution, lower_energy, upper_energy, ordered_tolerance,
        ):
            raise ValueError(f"{refinement}: nonfinite S/DIS refinement field")
        declared = next((
            value for value in SDIS_FIXED_FA_REFINEMENT_SLICES
            if isclose(target_f_a, value, abs_tol=1.0e-12)
        ), None)
        if declared is None:
            raise ValueError(f"{refinement}: undeclared fixed-fA slice")
        if declared in seen:
            raise ValueError(f"{refinement}: duplicate fixed-fA slice {declared}")
        seen.add(declared)
        if (
            source.get("status") != CLAIM_ACCEPTED
            or not truth(source.get("bracketed"))
            or not truth(source.get("resolution_reached"))
            or not truth(source.get("ordered_endpoint_accepted"))
            or source.get("ordered_phase") != "BCC"
            or source.get("ordered_side") != "upper"
            or not truth(source.get("decayed_lower"))
            or truth(source.get("decayed_upper"))
            or not isclose(target_f_a, root_f_a, abs_tol=1.0e-12)
            or lower >= upper
            or not lower <= root_chi_n <= upper
            or not isclose(upper - lower, width, abs_tol=1.0e-12)
            or resolution > SDIS_FIXED_FA_CHIN_TARGET + 1.0e-12
            or width > SDIS_FIXED_FA_CHIN_TARGET + 1.0e-12
            or lower_energy < -ordered_tolerance
            or upper_energy >= -ordered_tolerance
        ):
            raise ValueError(
                f"{refinement}: fixed-fA={declared:g} refinement gates failed"
            )
        overrides.append({
            "transition": "S/DIS",
            "chiN": root_chi_n,
            "fA": root_f_a,
            "bracket_width": width,
            "claim_kind": BOUNDARY_BRACKET,
            "status": CLAIM_ACCEPTED,
            "acceptance_basis": (
                "fixed-fA BVK2 BCC/DIS ordered-branch survival bracket; "
                "BCC 48^3 production grid; Delta chiN <= 0.01"
            ),
            "source": str(refinement.resolve()),
        })

    replaced_f_a = {float(row["fA"]) for row in overrides}
    kept = [
        row for row in rows
        if not (
            str(row["transition"]) == "S/DIS"
            and any(
                isclose(float(row["fA"]), f_a, abs_tol=1.0e-12)
                for f_a in replaced_f_a
            )
        )
    ]
    kept.extend(overrides)
    kept.sort(key=lambda item: (
        str(item["transition"]), float(item["chiN"] or 0.0)
    ))
    return kept, overrides


def apply_sc_fixed_fa_chi_refinement(
    rows: list[dict[str, object]],
    refinement: Path,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    """Promote all accepted signed fixed-composition BCC/CYL brackets."""
    refinement_files = (
        sorted(refinement.rglob("S_C_target_fA*.csv"))
        if refinement.is_dir()
        else [refinement]
    )
    overrides: list[dict[str, object]] = []
    for refinement_file in refinement_files:
        sources = read_rows(refinement_file)
        if not sources:
            continue
        if len(sources) != 1:
            raise ValueError(
                f"{refinement_file}: expected one S/C refinement row"
            )

        source = sources[0]
        if refinement.is_dir() and source.get("status") != CLAIM_ACCEPTED:
            continue
        transition = str(source.get("transition", ""))
        direction = str(source.get("direction", ""))
        target_f_a = finite_float(source.get("target_fA"))
        root_f_a = finite_float(source.get("root_fA"))
        root_chi_n = finite_float(source.get("root_chiN"))
        lower = finite_float(source.get("bracket_lower"))
        upper = finite_float(source.get("bracket_upper"))
        width = finite_float(source.get("bracket_width"))
        delta_lower = finite_float(source.get("delta_lower"))
        delta_upper = finite_float(source.get("delta_upper"))
        final_delta = finite_float(source.get("final_delta"))
        chi_n_tolerance = finite_float(source.get("chiN_tolerance"))
        residual_tolerance = finite_float(source.get("root_residual_tolerance"))
        bcc_n = finite_float(source.get("bcc_n"))
        cyl_n = finite_float(source.get("cyl_n"))
        if None in (
            target_f_a, root_f_a, root_chi_n, lower, upper, width,
            delta_lower, delta_upper, final_delta, chi_n_tolerance,
            residual_tolerance, bcc_n, cyl_n,
        ):
            raise ValueError(
                f"{refinement_file}: nonfinite S/C refinement field"
            )
        if (
            transition != "S/C"
            or direction != "chiN"
            or source.get("status") != CLAIM_ACCEPTED
            or not truth(source.get("phases_accepted"))
            or not truth(source.get("reference_covered"))
            or not isclose(target_f_a, root_f_a, abs_tol=1.0e-12)
            or lower >= upper
            or not lower <= root_chi_n <= upper
            or not isclose(upper - lower, width, abs_tol=1.0e-12)
            or width > SC_FIXED_FA_CHIN_TARGET + 1.0e-12
            or width > chi_n_tolerance + 1.0e-12
            or delta_lower * delta_upper > 0.0
            or abs(final_delta) > residual_tolerance
            or int(bcc_n) != 48
            or int(cyl_n) != 48
            or source.get("cyl_dims") != "48x84"
        ):
            raise ValueError(
                f"{refinement_file}: fixed-fA S/C refinement gates failed"
            )

        overrides.append({
            "transition": "S/C",
            "chiN": root_chi_n,
            "fA": root_f_a,
            "bracket_width": width,
            "claim_kind": BOUNDARY_BRACKET,
            "status": CLAIM_ACCEPTED,
            "acceptance_basis": (
                "fixed-fA BVK2 BCC/CYL signed free-energy bracket; "
                "BCC 48^3 and CYL 48x84 production grids; Delta chiN <= 0.05"
            ),
            "source": str(refinement_file.resolve()),
        })

    if not overrides:
        return rows, []
    by_f_a: dict[float, dict[str, object]] = {}
    for override in overrides:
        f_a = float(override["fA"])
        previous = by_f_a.get(f_a)
        if previous is None or float(override["bracket_width"]) < float(
            previous["bracket_width"]
        ):
            by_f_a[f_a] = override
    overrides = list(by_f_a.values())
    kept = [
        row for row in rows
        if not (str(row["transition"]) == "S/C" and any(
            isclose(float(row["fA"]), float(override["fA"]), abs_tol=1.0e-12)
            for override in overrides
        ))
    ]
    kept.extend(overrides)
    kept.sort(key=lambda item: (
        str(item["transition"]), float(item["chiN"] or 0.0)
    ))
    return kept, overrides


def write_wide_csv(path: Path, rows: list[dict[str, object]]) -> None:
    """Write rows carrying fields outside the fixed publication schema."""
    if not rows:
        return
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def critical_endpoint(rows: list[dict[str, object]]) -> dict[str, object] | None:
    """The f_A = 0.5 endpoint computed by the BVK2 L/DIS unit."""
    for row in rows:
        if row["transition"] == "L/DIS" and row["fA"] is not None:
            if abs(float(row["fA"]) - 0.5) <= 1.0e-9:
                return {
                    "transition": "S/DIS",
                    "chiN": row["chiN"],
                    "fA": 0.5,
                    "direction": row["direction"],
                    "bracket_width": row["bracket_width"],
                    "final_delta": "",
                    "claim_kind": BOUNDARY_BRACKET,
                    "status": "accepted",
                    "acceptance_basis": "bvk2_critical_endpoint",
                    "source": row["source"],
                }
    return None


def write_summary(path: Path, rows, rejected, endpoint, fcc_triple) -> None:
    counts: dict[str, dict[str, int]] = {}
    for row in rows:
        transition = str(row["transition"])
        status = str(row["status"])
        counts.setdefault(transition, {})
        counts[transition][status] = (
            counts[transition].get(status, 0) + 1
        )
    lines = [
        "# BVK2 publication phase diagram",
        "",
        "Accepted BVK2 boundary roots assembled from the campaign units. The",
        "computed `f_A <= 0.5` half is reflected about `f_A = 0.5` using A/B",
        "diblock symmetry, with symbols only on computed roots. Dashed",
        "neutral-gray SCFT literature curves are a comparison overlay only:",
        "they carry no symbols and are not BVK2 calculations.",
        "The neutral-gray dotted `RPA stability limit` is the AB-diblock spinodal",
        "computed directly with `Polyorder.RPA.compute_stability_limit`; it is",
        "a linear-stability approximation to the ODT, not a BVK2 boundary.",
        "Its `chiN*`, `D*`, and `k*` data are in `rpa_stability_limit.csv`.",
        "`boundary_candidates.csv` contains every accepted/provisional claim;",
        "`accepted_phase_boundaries.csv` contains accepted boundary claims only.",
        "`phase_diagram_plot_data.csv` is the canonical stable-topology plot table;",
        "it excludes pairwise continuations outside the displayed phase pockets.",
        "Rebuild scripts are stored beside the figures: `phase_diagram.py`,",
        "`fcc_region_zoom.py`, and `o70_region_zoom.py`. Each launcher resolves",
        "the repository root from its own location and invokes the canonical",
        "source under `scripts/`. Full assembly writes the canonical plot table,",
        "then the shared ACS/mpltex snapshot renderer produces Figure 6 without",
        "duplicating or reconstructing the scientific boundary selection.",
        "The low-segregation O70 pocket is shown separately in",
        "`o70_region_zoom.svg` (`0.4 <= f_A <= 0.5`, `10 <= chiN <= 15`).",
        "Its SCFT reference paths use display-only monotone smoothing to remove",
        "digitization stair steps; source vertices and curve endpoints are unchanged.",
        "The accepted chiN=12 free-energy evidence for this pocket is in",
        "`../bvk2_fcc_o70_boundary_campaign/chi12_phase_stability/`.",
        "The FCC/BCC/DIS pocket closes at the display-only estimate",
        f"`f_A={fcc_triple[0]:.9f}, chiN={fcc_triple[1]:.9f}` derived in",
        "`fcc_region_zoom_triple_point.csv`. The full plot uses only strict",
        "matched-dx FCC edges, terminates BCC/DIS at this point, and omits",
        "FCC/DIS below it; excluded pairwise continuations remain diagnostic.",
        "",
        "| boundary | accepted roots | provisional roots |",
        "| --- | ---: | ---: |",
    ]
    for transition in TRANSITIONS:
        lines.append(
            f"| {transition} | "
            f"{counts.get(transition, {}).get('accepted', 0)} | "
            f"{counts.get(transition, {}).get('provisional', 0)} |"
        )
    lines += [
        "",
        f"Critical endpoint: "
        + (
            f"`f_A = 0.5, chiN = {endpoint['chiN']:.6f}` (BVK2 L/DIS unit)"
            if endpoint
            else "not yet computed"
        ),
        "",
        f"Gate-rejected point rows: {len(rejected)} "
        "(recorded in `rejected_points.csv`, never plotted).",
        "",
        "Constitutive scalar `c2 = 0.16` is frozen by",
        "`results/bvk2_c2_calibration_nx128/selected_model.toml`.",
        "",
        "For `chiN >= 40`, S/C roots are replaced by the validated",
        "Uneyama--Doi theta-before-discretization ledger at",
        "`results/bvk2_ud_theta_sc_highchi/accepted_boundaries.csv`.",
        "Those points use 3x optimization and signed 4x endpoint audits.",
        "Near the symmetric ODT roof, completed fixed-composition BCC/DIS",
        "refinements at `f_A=0.40,0.41,0.43...0.49` replace their legacy roots as each",
        "accepted artifact lands. They use BCC 48^3 and `Delta chiN <= 0.01`;",
        "unfinished slices retain their previous accepted publication points.",
        "",
        "High-chiN G/C and G/L roots are replaced by the validated",
        "Uneyama--Doi theta-before-discretization G-pocket ledgers at",
        "`results/bvk2_ud_theta_g_pocket_highchi_ledger/`.",
        "At `chiN=54`, the direct G/C verification in",
        "`results/bvk2_ud_theta_gc_chi54_direct/boundary_point.csv`",
        "takes precedence and uses direct GYR 112^3 and CYL 96x168",
        "factor-3 energies.",
        "At `chiN=30,32,34,37,40,42,45,48,50`, the completed global-filter",
        "G/C production ledger takes final precedence. It uses one fixed",
        "adaptive-curvature sensor-filter ratio `r=6`, direct factor-3",
        "GYR 96^3 and CYL 128x224 energies, stress-free cells, phase-identity",
        "checks, and signed fixed-chiN brackets no wider than 0.002 in f_A.",
        "The canonical source is",
        "`results/bvk2_ud_theta_gc_filtered_r6/boundary_points.csv`.",
        "Legacy G/L points in that ledger use",
        "GYR 96^3/112^3 Richardson correction and direct LAM N=1024.",
        "Matching rows in",
        "`results/bvk2_ud_theta_gl_direct_fine/boundary_candidates.csv`",
        "take precedence and use direct GYR 112^3 and LAM N=1024 energies.",
        "These roots select the lower-free-energy one-period LAM basin, use",
        "signed factor-3 brackets, and apply no Richardson extrapolation.",
        "Accepted matched-spacing promotion rows normally take final precedence at",
        "matching chiN slices; the canonical input is `matched_boundary_promotions.csv`",
        "under `results/bvk2_ud_theta_gl_matched_spacing_audit/`.",
        "Each promotion uses hashed accepted endpoint records,",
        "one fixed GYR/LAM grid pair across its bracket, and <=3% physical",
        "spacing mismatch.",
        "At `chiN=15`, the direct GYR64/LAM1024 G/L bracket instead takes final",
        "precedence because the adjacent O70/L campaign consistently uses LAM1024",
        "and its audited O70 edges have crossed. The closure record is",
        "`results/bvk2_o70_gl_triple_chi15_matched/thermodynamic_audit.csv`.",
        "",
        "For `chiN <= 42`, the direct theta GYR/LAM verification campaign",
        "takes precedence over legacy Richardson roots. Legacy coordinates",
        "remain provisional locators until a direct narrow signed bracket",
        "passes both endpoint gates.",
        "The common claim/status contract is documented in",
        "`docs/phase_boundary_acceptance_contract.md`: accepted phase",
        "endpoints do not become solid boundary points until the deterministic",
        "boundary assembler emits `claim_kind=boundary_bracket,status=accepted`.",
        "",
        "Pairwise C/L roots are retained in `boundary_diagnostics.csv` as",
        "metastable evidence only. They are excluded from the thermodynamic",
        "publication diagram because the stable GYR pocket is bounded by",
        "G/C and G/L.",
    ]
    path.write_text("\n".join(lines) + "\n")


def write_output_rebuild_scripts(outdir: Path) -> None:
    """Write self-locating figure rebuild entry points beside the outputs."""
    launchers = {
        "phase_diagram.py": ("assemble_bvk2_publication_phase_diagram.py", False),
        "fcc_region_zoom.py": ("render_bvk2_fcc_zoom.py", True),
        "o70_region_zoom.py": ("render_bvk2_o70_zoom.py", True),
    }
    for filename, (source_name, use_uv) in launchers.items():
        command = (
            '["uv", "run", "--frozen", "python", str(source)]'
            if use_uv
            else '[sys.executable, str(source)]'
        )
        text = f'''#!/usr/bin/env python3
"""Rebuild the colocated BVK2 publication figure from canonical sources."""

from pathlib import Path
import subprocess
import sys


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
source = PROJECT / "scripts" / "{source_name}"
if not source.is_file():
    raise SystemExit(f"canonical renderer not found: {{source}}")
subprocess.run({command}, cwd=PROJECT, check=True)
'''
        launcher = outdir / filename
        launcher.write_text(text)
        launcher.chmod(0o755)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign-dir", type=Path, default=DEFAULT_CAMPAIGN)
    parser.add_argument("--sdis-dir", type=Path, default=DEFAULT_SDIS)
    parser.add_argument(
        "--sdis-fixed-fa-refinement",
        type=Path,
        default=DEFAULT_SDIS_FIXED_FA_REFINEMENT,
        help="accepted fixed-fA BCC/DIS chiN-refinement campaign",
    )
    parser.add_argument(
        "--sc-fixed-fa-refinement",
        type=Path,
        default=DEFAULT_SC_FIXED_FA_REFINEMENT,
        help="accepted fixed-fA BCC/CYL signed chiN-refinement file or tree",
    )
    parser.add_argument(
        "--scft-reference", type=Path, default=DEFAULT_SCFT_REFERENCE
    )
    parser.add_argument(
        "--scft-source-svg", type=Path, default=DEFAULT_SCFT_SOURCE_SVG
    )
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument(
        "--corrections", type=Path, default=DEFAULT_CORRECTIONS,
        help="root holding Richardson corrected_roots.csv files; applied "
        "per transition only when every accepted point of that transition "
        "has one",
    )
    parser.add_argument(
        "--ud-theta-sc", type=Path, default=DEFAULT_UD_THETA_SC,
        help="validated high-chiN UD-theta S/C ledger",
    )
    parser.add_argument(
        "--ud-theta-gc", type=Path, default=DEFAULT_UD_THETA_GC,
        help="validated high-chiN UD-theta G/C boundary ledger",
    )
    parser.add_argument(
        "--ud-theta-gc-direct",
        type=Path,
        default=DEFAULT_UD_THETA_GC_DIRECT,
        help="direct GYR112/CYL96x168 G/C boundary ledger",
    )
    parser.add_argument(
        "--ud-theta-gc-filtered-r6",
        type=Path,
        default=DEFAULT_UD_THETA_GC_FILTERED_R6,
        help="complete global-filter r=6 GYR96/CYL128x224 G/C ledger",
    )
    parser.add_argument(
        "--ud-theta-gl", type=Path, default=DEFAULT_UD_THETA_GL,
        help="validated high-chiN UD-theta G/L boundary ledger",
    )
    parser.add_argument(
        "--ud-theta-gl-direct",
        type=Path,
        default=DEFAULT_UD_THETA_GL_DIRECT,
        help="direct GYR112/LAM1024 G/L boundary ledger",
    )
    parser.add_argument(
        "--ud-theta-gl-matched",
        type=Path,
        default=DEFAULT_UD_THETA_GL_MATCHED,
        help="accepted matched-physical-spacing GYR/LAM G/L promotion ledger",
    )
    parser.add_argument(
        "--o70-gl-chi15-audit",
        type=Path,
        default=DEFAULT_O70_GL_CHI15_AUDIT,
        help="LAM1024-consistent chiN=15 O70 closure audit",
    )
    parser.add_argument(
        "--ud-theta-gl-lowchi",
        type=Path,
        default=DEFAULT_UD_THETA_GL_LOWCHI,
        help="direct low-chi GYR/LAM G/L verification root",
    )
    parser.add_argument(
        "--fcc-o70", type=Path, default=DEFAULT_FCC_O70,
        help="accepted signed FCC and O70 boundary-bracket ledger",
    )
    parser.add_argument(
        "--strict-fcc", type=Path, default=DEFAULT_STRICT_FCC,
        help="accepted strict matched-dx FCC boundary ledger",
    )
    parser.add_argument(
        "--strict-fcc-points", type=Path, default=DEFAULT_STRICT_FCC_POINTS,
        help="matched-dx FCC energy samples used for triple-point closure",
    )
    parser.add_argument(
        "--rpa-stability", type=Path,
        help="existing Polyorder.RPA CSV; omit to recompute it in the outdir",
    )
    parser.add_argument(
        "--polyorder-project", type=Path, default=DEFAULT_POLYORDER_PROJECT,
        help="Julia environment containing Polyorder, Optim, and Roots",
    )
    parser.add_argument(
        "--rpa-script", type=Path, default=DEFAULT_RPA_SCRIPT,
        help="Julia script that computes the Polyorder.RPA stability limit",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    rpa_path = (
        args.rpa_stability.resolve()
        if args.rpa_stability is not None
        else (args.outdir / RPA_FILENAME).resolve()
    )
    if args.rpa_stability is None:
        compute_rpa_stability_limit(
            rpa_path,
            project=args.polyorder_project,
            script=args.rpa_script,
        )
    rpa_stability_limit = load_rpa_stability_limit(rpa_path)

    rows, rejected, diagnostics = collect(
        args.campaign_dir.resolve(), args.sdis_dir.resolve()
    )
    endpoint = critical_endpoint(rows)
    # L/DIS exists only to supply the shared critical endpoint; it is not a
    # separately plotted boundary (BCC is the first ordered phase above DIS).
    rows = [row for row in rows if row["transition"] != "L/DIS"]
    if endpoint is not None:
        rows.append(endpoint)

    # Richardson grid corrections: per point; corrected rows stay "accepted"
    # (solid) while uncorrected rows are demoted to the "provisional" display
    # tier (dashed, disconnected) so the two are never joined on one curve.
    corrections = load_corrections(args.corrections.resolve())
    raw_by_key = {(str(r["transition"]), float(r["chiN"])): r["fA"]
                  for r in rows if r["chiN"] is not None}
    rows, correction_report = apply_corrections(rows, corrections)
    for transition, note in sorted(correction_report.items()):
        print(f"  grid correction | {transition:<6} {note}", flush=True)
    # Keep the uncorrected value alongside every corrected row.
    for diag in diagnostics:
        key = (str(diag["transition"]), float(diag["chiN"]))
        if key in corrections:
            diag["fA_raw"] = raw_by_key.get(key, "")
            diag["fA_corrected"] = corrections[key].get("root_fA_corrected", "")
            diag["grid_shift"] = corrections[key].get("shift", "")
            diag["dDelta_correction"] = corrections[key].get(
                "dDelta_correction", "")

    rows, sdis_fixed_fa_refinement = apply_sdis_fixed_fa_chi_refinement(
        rows, args.sdis_fixed_fa_refinement.resolve()
    )
    refined_sdis_f_a = {float(row["fA"]) for row in sdis_fixed_fa_refinement}
    diagnostics = [
        diag for diag in diagnostics
        if not (
            str(diag["transition"]) == "S/DIS"
            and any(
                isclose(float(diag["fA"]), f_a, abs_tol=1.0e-12)
                for f_a in refined_sdis_f_a
            )
        )
    ]
    diagnostics.extend({
        **row,
        "direction": "chiN",
        "final_delta": "",
        "reference_source": "",
        "gate_reason": "accepted_fixed_fa_bcc_dis_chin_refinement",
    } for row in sdis_fixed_fa_refinement)

    rows, sc_fixed_fa_refinement = apply_sc_fixed_fa_chi_refinement(
        rows, args.sc_fixed_fa_refinement.resolve()
    )
    refined_sc_f_a = {float(row["fA"]) for row in sc_fixed_fa_refinement}
    diagnostics = [
        diag for diag in diagnostics
        if not (
            str(diag["transition"]) == "S/C"
            and any(
                isclose(float(diag["fA"]), f_a, abs_tol=1.0e-12)
                for f_a in refined_sc_f_a
            )
        )
    ]
    diagnostics.extend({
        **row,
        "direction": "chiN",
        "final_delta": "",
        "reference_source": "",
        "gate_reason": "accepted_fixed_fa_bcc_cyl_signed_chin_refinement",
    } for row in sc_fixed_fa_refinement)

    rows, ud_theta_sc = apply_ud_theta_sc_overrides(
        rows, args.ud_theta_sc.resolve()
    )
    diagnostics = [
        diag for diag in diagnostics
        if not (
            str(diag["transition"]) == "S/C"
            and any(
                float(diag["chiN"]) == float(row["chiN"])
                for row in ud_theta_sc
            )
        )
    ]
    diagnostics.extend({
        **row,
        "direction": "fA",
        "final_delta": "",
        "reference_source": "",
        "gate_reason": "accepted_ud_theta_signed_4x_bracket",
    } for row in ud_theta_sc)

    rows, ud_theta_gc = apply_ud_theta_order_order_override(
        rows,
        args.ud_theta_gc.resolve(),
        source_boundary="G_C",
        transition="G/C",
    )
    diagnostics = [
        diag for diag in diagnostics
        if not (
            str(diag["transition"]) == "G/C"
            and any(
                float(diag["chiN"]) == float(row["chiN"])
                for row in ud_theta_gc
            )
        )
    ]
    diagnostics.extend({
        **row,
        "direction": "fA",
        "final_delta": "",
        "reference_source": "",
        "gate_reason": "accepted_ud_theta_signed_audit_bracket",
    } for row in ud_theta_gc)

    rows, ud_theta_gc_direct = apply_ud_theta_direct_gc_override(
        rows, args.ud_theta_gc_direct.resolve()
    )
    diagnostics = [
        diag for diag in diagnostics
        if not (
            str(diag["transition"]) == "G/C"
            and any(
                float(diag["chiN"]) == float(row["chiN"])
                for row in ud_theta_gc_direct
            )
        )
    ]
    diagnostics.extend({
        **row,
        "direction": "fA",
        "final_delta": "",
        "reference_source": "",
        "gate_reason": "accepted_ud_theta_direct_gyr112_cyl96x168",
    } for row in ud_theta_gc_direct)

    # The completed global-filter campaign is the final G/C authority at its
    # nine declared chiN slices.  Apply it after legacy and direct overrides.
    rows, ud_theta_gc_filtered_r6 = apply_ud_theta_filtered_r6_gc_override(
        rows, args.ud_theta_gc_filtered_r6.resolve()
    )
    filtered_r6_keys = {
        (str(row["transition"]), float(row["chiN"]))
        for row in ud_theta_gc_filtered_r6
    }
    diagnostics = [
        diag for diag in diagnostics
        if (str(diag["transition"]), float(diag["chiN"]))
        not in filtered_r6_keys
    ]
    diagnostics.extend({
        **row,
        "direction": "fA",
        "final_delta": "",
        "reference_source": "",
        "gate_reason": "accepted_ud_theta_global_filter_r6_gyr96_cyl128x224",
    } for row in ud_theta_gc_filtered_r6)

    rows, ud_theta_gl = apply_ud_theta_corrected_gl_override(
        rows, args.ud_theta_gl.resolve()
    )
    diagnostics = [
        diag for diag in diagnostics
        if not (
            str(diag["transition"]) == "G/L"
            and any(
                float(diag["chiN"]) == float(row["chiN"])
                for row in ud_theta_gl
            )
        )
    ]
    diagnostics.extend({
        **row,
        "direction": "fA",
        "final_delta": "",
        "reference_source": "",
        "gate_reason": "accepted_ud_theta_signed_corrected_bracket",
    } for row in ud_theta_gl)

    rows, ud_theta_gl_lowchi = apply_ud_theta_lowchi_gl_status(
        rows, args.ud_theta_gl_lowchi.resolve()
    )
    diagnostics = [
        diag for diag in diagnostics
        if not (
            str(diag["transition"]) == "G/L"
            and any(
                float(diag["chiN"]) == float(row["chiN"])
                for row in ud_theta_gl_lowchi
            )
        )
    ]
    diagnostics.extend({
        **row,
        "direction": "fA",
        "final_delta": "",
        "reference_source": "",
        "gate_reason": (
            "accepted_ud_theta_lowchi_direct"
            if row["status"] == "accepted"
            else "provisional_ud_theta_lowchi_revalidation"
        ),
    } for row in ud_theta_gl_lowchi)

    rows, ud_theta_gl_direct = apply_ud_theta_direct_gl_override(
        rows, args.ud_theta_gl_direct.resolve()
    )
    diagnostics = [
        diag for diag in diagnostics
        if not (
            str(diag["transition"]) == "G/L"
            and any(
                float(diag["chiN"]) == float(row["chiN"])
                for row in ud_theta_gl_direct
            )
        )
    ]
    diagnostics.extend({
        **row,
        "direction": "fA",
        "final_delta": "",
        "reference_source": "",
        "gate_reason": (
            "accepted_ud_theta_direct_gyr112_lam1024"
        ),
    } for row in ud_theta_gl_direct)

    fcc_o70 = [
        *load_fcc_o70_boundaries(args.fcc_o70.resolve()),
        *load_strict_fcc_boundaries(args.strict_fcc.resolve()),
    ]
    rows.extend(fcc_o70)

    # Matched-spacing promotions are the normal final G/L authority.  Apply
    # them after direct and FCC/O70 inputs so legacy same-chi roots cannot be
    # reintroduced.  The LAM1024-consistent chiN=15 O70 closure is the one
    # deliberate exception and is applied immediately afterward.
    rows, ud_theta_gl_matched = apply_ud_theta_matched_gl_override(
        rows, args.ud_theta_gl_matched.resolve()
    )
    matched_keys = {
        (str(row["transition"]), float(row["chiN"]))
        for row in ud_theta_gl_matched
    }
    fcc_o70 = [
        row for row in fcc_o70
        if (str(row["transition"]), float(row["chiN"])) not in matched_keys
    ]
    diagnostics = [
        diag for diag in diagnostics
        if (str(diag["transition"]), float(diag["chiN"])) not in matched_keys
    ]
    diagnostics.extend({
        **row,
        "direction": "fA",
        "final_delta": "",
        "reference_source": "",
        "gate_reason": "accepted_ud_theta_matched_spacing_gyr_lam",
    } for row in ud_theta_gl_matched)

    rows, o70_gl_chi15 = apply_o70_chi15_lam1024_gl_override(
        rows,
        args.o70_gl_chi15_audit.resolve(),
        args.fcc_o70.resolve(),
    )
    chi15_keys = {
        (str(row["transition"]), float(row["chiN"]))
        for row in o70_gl_chi15
    }
    fcc_o70 = [
        row for row in fcc_o70
        if (str(row["transition"]), float(row["chiN"])) not in chi15_keys
    ]
    fcc_o70.extend(o70_gl_chi15)

    rows, metastable_o70 = filter_metastable_o70_pairwise(
        rows, open_o70_slices(args.fcc_o70.resolve())
    )
    stable_fcc_o70 = [row for row in fcc_o70 if row not in metastable_o70]
    stable_fcc_o70_keys = {
        (str(row["transition"]), float(row["chiN"]))
        for row in stable_fcc_o70
    }
    metastable_keys = {
        (str(row["transition"]), float(row["chiN"]))
        for row in metastable_o70
    }
    diagnostics = [
        row for row in diagnostics
        if (str(row["transition"]), float(row["chiN"]))
        not in metastable_keys | stable_fcc_o70_keys
    ]
    diagnostics.extend({
        **row,
        "direction": "fA",
        "final_delta": "",
        "reference_source": (
            str(args.o70_gl_chi15_audit.resolve())
            if (str(row["transition"]), float(row["chiN"])) in chi15_keys
            else ""
        ),
        "gate_reason": (
            "accepted_o70_closure_direct_gyr64_lam1024"
            if (str(row["transition"]), float(row["chiN"])) in chi15_keys
            else "accepted_signed_fcc_o70_bracket"
        ),
    } for row in stable_fcc_o70)
    diagnostics.extend({
        **row,
        "status": "metastable",
        "direction": "fA",
        "final_delta": "",
        "reference_source": "",
        "gate_reason": (
            "metastable_pairwise_g_l_inside_stable_o70_interval"
            if row["transition"] == "G/L"
            else "metastable_pairwise_o70_edge_not_thermodynamic"
        ),
    } for row in metastable_o70)

    # C/L is a useful pairwise crossing but not a stable thermodynamic
    # boundary while the GYR pocket lies between CYL and LAM.
    rows = [row for row in rows if row["transition"] != "C/L"]
    for diag in diagnostics:
        if diag["transition"] == "C/L":
            diag["status"] = "metastable"
            diag["gate_reason"] = "metastable_pairwise_only"

    write_csv(args.outdir / "boundary_candidates.csv", rows)
    write_csv(
        args.outdir / "accepted_phase_boundaries.csv",
        [row for row in rows if row["status"] == "accepted"],
    )
    write_wide_csv(args.outdir / "boundary_diagnostics.csv", diagnostics)
    write_wide_csv(args.outdir / "rejected_points.csv", rejected)

    scft_rows = load_scft_reference(args.scft_reference)
    write_scft_reference_csv(
        args.outdir / "scft_reference_boundaries.csv",
        scft_rows,
        args.scft_reference,
        args.scft_source_svg,
        model_label="BVK2",
    )

    subprocess.run(
        [
            "uv",
            "run",
            "--frozen",
            "python",
            str(HERE / "render_bvk2_fcc_zoom.py"),
            "--boundaries",
            str(args.outdir / "accepted_phase_boundaries.csv"),
            "--legacy-fcc",
            str(args.fcc_o70.resolve()),
            "--strict-fcc",
            str(args.strict_fcc.resolve()),
            "--strict-points",
            str(args.strict_fcc_points.resolve()),
            "--scft",
            str(args.outdir / "scft_reference_boundaries.csv"),
            "--rpa-stability",
            str(rpa_path),
            "--output-prefix",
            str(args.outdir / "fcc_region_zoom"),
        ],
        cwd=PROJECT,
        check=True,
    )

    fcc_triple = load_fcc_triple_point(args.outdir / FCC_TRIPLE_FILENAME)
    plot_rows = stable_fcc_plot_rows(rows, fcc_triple[1])
    write_csv(args.outdir / "phase_diagram_plot_data.csv", plot_rows)

    estimated_triple_points = (
        fcc_triple,
        load_estimated_triple_point(
            args.outdir / LOWER_O70_TRIPLE_FILENAME,
            ("G/C", "G/O70", "C/O70"),
        ),
        load_estimated_triple_point(
            args.outdir / UPPER_O70_TRIPLE_FILENAME,
            ("G/L", "G/O70", "O70/L"),
        ),
    )

    render_png(
        args.outdir / "phase_diagram.png",
        plot_rows,
        [],
        scft_rows,
        title=TITLE,
        transitions=TRANSITIONS,
        critical_endpoint_exclusions=EXTENDED_TRANSITIONS,
        estimated_triple_points=estimated_triple_points,
        chi_n_max=50.0,
        rpa_stability_limit=rpa_stability_limit,
    )
    # Figure 6 uses the shared ACS/mpltex publication renderer. Scientific
    # assembly is complete before this call; the renderer consumes only the
    # canonical stable plot table and topology-junction records written above.
    subprocess.run(
        [
            "uv",
            "run",
            "--frozen",
            "python",
            str(HERE / "render_bvk2_phase_diagram_snapshot.py"),
            "--plot-data",
            str(args.outdir / "phase_diagram_plot_data.csv"),
            "--scft",
            str(args.outdir / "scft_reference_boundaries.csv"),
            "--rpa",
            str(rpa_path),
            "--fcc-triple",
            str(args.outdir / FCC_TRIPLE_FILENAME),
            "--g-c-o70-triple",
            str(args.outdir / LOWER_O70_TRIPLE_FILENAME),
            "--g-o70-l-triple",
            str(args.outdir / UPPER_O70_TRIPLE_FILENAME),
            "--output",
            str(args.outdir / "phase_diagram.svg"),
        ],
        cwd=PROJECT,
        check=True,
    )

    subprocess.run(
        [
            "uv",
            "run",
            "--frozen",
            "python",
            str(HERE / "render_bvk2_o70_zoom.py"),
            "--boundaries",
            str(args.outdir / "phase_diagram_plot_data.csv"),
            "--brackets",
            str(args.fcc_o70.resolve()),
            "--scft",
            str(args.outdir / "scft_reference_boundaries.csv"),
            "--rpa-stability",
            str(rpa_path),
            "--triple-point",
            str(args.outdir / LOWER_O70_TRIPLE_FILENAME),
            "--upper-triple-point",
            str(args.outdir / UPPER_O70_TRIPLE_FILENAME),
            "--output-prefix",
            str(args.outdir / "o70_region_zoom"),
        ],
        cwd=PROJECT,
        check=True,
    )
    subprocess.run(
        [
            "uv",
            "run",
            "--frozen",
            "python",
            str(HERE / "render_bvk2_o70_chi12_stability.py"),
        ],
        cwd=PROJECT,
        check=True,
    )

    write_summary(
        args.outdir / "phase_diagram.md",
        rows,
        rejected,
        endpoint,
        fcc_triple,
    )
    write_output_rebuild_scripts(args.outdir)
    accepted_count = sum(row["status"] == "accepted" for row in rows)
    provisional_count = sum(
        row["status"] == "provisional" for row in rows
    )
    print(
        f"assembled {accepted_count} accepted and "
        f"{provisional_count} provisional BVK2 roots "
        f"({len(rejected)} rejected) -> {args.outdir}",
        flush=True,
    )


if __name__ == "__main__":
    main()
