#!/usr/bin/env python3
"""Render every Macromolecules main-text and SI figure with Matplotlib.

The numerical generators remain responsible for producing the machine-readable
CSV tables.  This script is the publication rendering layer: it reads only
frozen tables/source figures, preserves the scientific data, writes
deterministic SVG, and uses caption-led journal styling
without figure-level titles or explanatory subtitles.

Selected quantized SCFT overlay paths use bounded, endpoint-preserving
display smoothing. Their source tables and all computed AQCE roots stay intact.
AQCE boundary points are joined by straight-line segments.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams.update(
    {
        "svg.fonttype": "none",
        "svg.hashsalt": "macromolecules-publication-figures-v1",
    }
)
import matplotlib.pyplot as plt
import mpltex
import mpltex.acs as mpltex_acs
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import FormatStrFormatter, MultipleLocator

try:
    from .scft_display import smooth_scft_display_curve
    from .assemble_burp_ti_publication_phase_diagram import (
        boundary_geometry as phase_boundary_geometry,
        rows_with_estimated_triple_points,
    )
except ImportError:  # Direct execution: python scripts/render_macromolecules_figures.py
    from scft_display import smooth_scft_display_curve
    from assemble_burp_ti_publication_phase_diagram import (
        boundary_geometry as phase_boundary_geometry,
        rows_with_estimated_triple_points,
    )


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
SVG_METADATA = {"Date": None, "Creator": "Matplotlib + mpltex ACS via uv"}

FIGURE6_DIR = RESULTS / "bvk2_publication_phase_diagram"
FIGURE6_SOURCE_SVG = FIGURE6_DIR / "phase_diagram.svg"
FIGURE6_PLOT_DATA = FIGURE6_DIR / "phase_diagram_plot_data.csv"
FIGURE6_SCFT_DATA = FIGURE6_DIR / "scft_reference_boundaries.csv"
FIGURE6_RPA_DATA = FIGURE6_DIR / "rpa_stability_limit.csv"
FIGURE6_TITLE = "AQCE continuous phase diagram"
FIGURE6_TRANSITIONS = (
    "S/C",
    "S/DIS",
    "G/C",
    "G/L",
    "FCC/BCC",
    "FCC/DIS",
    "C/O70",
    "G/O70",
    "O70/L",
)
FIGURE6_CRITICAL_ENDPOINT_EXCLUSIONS = (
    "FCC/BCC",
    "FCC/DIS",
    "C/O70",
    "G/O70",
    "O70/L",
)
FIGURE6_COLORS = {
    "S/DIS": "#1f77b4",
    "S/C": "#d62728",
    "G/C": "#9467bd",
    "G/L": "#ff7f0e",
    "FCC/BCC": "#8c564b",
    "FCC/DIS": "#e377c2",
    "C/O70": "#17becf",
    "G/O70": "#bcbd22",
    "O70/L": "#393b79",
}
FIGURE6_TRIPLE_POINT_SOURCES = (
    (
        FIGURE6_DIR / "fcc_region_zoom_triple_point.csv",
        "fA_estimate",
        "chiN_estimate",
        ("S/DIS", "FCC/DIS", "FCC/BCC"),
    ),
    (
        FIGURE6_DIR / "estimated_g_c_o70_triple_point.csv",
        "fA",
        "chiN",
        ("G/C", "G/O70", "C/O70"),
    ),
    (
        FIGURE6_DIR / "estimated_g_o70_l_triple_point.csv",
        "fA",
        "chiN",
        ("G/L", "G/O70", "O70/L"),
    ),
)

ACS_COLORS = mpltex.tableau_10
# One Tableau palette: orange QCE/red AQCE form the warm proposed-model pair,
# contrasted with blue OPF, purple BURP, green UD, and neutral-gray SCFT.
COLORS = {
    "scft": ACS_COLORS[7],
    "exact_rpa": mpltex.almost_black,
    "ohta_kawasaki": ACS_COLORS[7],
    "ok": ACS_COLORS[1],
    "uneyama_doi": ACS_COLORS[2],
    "liu2019_opf": ACS_COLORS[0],
    "burp_ti": ACS_COLORS[3],
    "bvk1": ACS_COLORS[9],
    "bvk2": ACS_COLORS[1],
    "bvk2_fixed": ACS_COLORS[4],
    "bvk1_bvk2": ACS_COLORS[9],
    "S/DIS": ACS_COLORS[0],
    "S/C": ACS_COLORS[1],
    "G/C": ACS_COLORS[3],
    "G/L": ACS_COLORS[2],
}

SCFT_CRITICAL_CHI = 10.537542662116

LAMELLAR_AGGREGATE_MODELS = (
    "ohta_kawasaki",
    "uneyama_doi",
    "liu2019_opf",
    "burp_ti",
    "bvk2",
)
LAMELLAR_AGGREGATE_PLOT_MODELS = (
    "ohta_kawasaki",
    "uneyama_doi",
    "liu2019_opf",
    "burp_ti",
    "bvk2_fixed",
    "bvk2",
)
LAMELLAR_AGGREGATE_SCHEMA = "liu2019-stress-free-period-map-v2"
LAMELLAR_AGGREGATE_CLAIM_KIND = "stress_free_lamellar_observable"
LAMELLAR_AGGREGATE_SOURCE_FINGERPRINTS = {
    "ohta_kawasaki": {
        "model_label": frozenset({"Ohta-Kawasaki (recomputed)"}),
        "calibration_role": frozenset(
            {"liu2019_ok_asymptotic_quadratic_mapped_c3_c4"}
        ),
        "period_optimizer": frozenset({"StressRootBisection"}),
        "reduced_field_protocol": frozenset(
            {"mean-constrained Fourier nonlinear conjugate gradient"}
        ),
        "grid_count": frozenset({"128"}),
        "quadrature_order": frozenset({""}),
        "discretization_schema": frozenset({""}),
    },
    "uneyama_doi": {
        "model_label": frozenset({"Uneyama-Doi"}),
        "calibration_role": frozenset({"not_applicable"}),
        "period_optimizer": frozenset({"BootstrapLocal"}),
        "reduced_field_protocol": frozenset(
            {"mean-constrained fused-spectral L-BFGS with profile continuation"}
        ),
        "grid_count": frozenset({"128"}),
        "quadrature_order": frozenset({""}),
        "discretization_schema": frozenset({""}),
    },
    "liu2019_opf": {
        "model_label": frozenset({"OPF (force/stress mapped)"}),
        "calibration_role": frozenset({"per_state_force_stress_mapping"}),
        "period_optimizer": frozenset({"StressRootBisection"}),
        "reduced_field_protocol": frozenset(
            {"mean-constrained Fourier nonlinear conjugate gradient"}
        ),
        "grid_count": frozenset({"128"}),
        "quadrature_order": frozenset({""}),
        "discretization_schema": frozenset({""}),
        "opf_c2": frozenset({"-1.0"}),
    },
    "burp_ti": {
        "model_label": frozenset({"BURP-TI"}),
        "calibration_role": frozenset({"not_applicable"}),
        "period_optimizer": frozenset(
            {"BoundedKKTNewtonGMRES+AnalyticStressRoot"}
        ),
        "reduced_field_protocol": frozenset({"bounded-KKT Newton-GMRES"}),
        "grid_count": frozenset({"128"}),
        "quadrature_order": frozenset({"64"}),
        "discretization_schema": frozenset({""}),
    },
    "bvk2": {
        "model_label": frozenset({"BVK2"}),
        "calibration_role": frozenset(
            {"training", "holdout", "frozen_c2_evaluation"}
        ),
        "period_optimizer": frozenset({"FilteredUDThetaAnalyticCellStressRoot"}),
        "reduced_field_protocol": frozenset(
            {"oversampled filtered UD-theta L-BFGS with cosine-Newton polish"}
        ),
        "grid_count": frozenset({"256"}),
        "quadrature_order": frozenset({""}),
        "discretization_schema": frozenset(
            {"ud-theta-spectral-adaptive-k-physical-sensor-filter-oversampled-v2"}
        ),
        "c2": frozenset({"0.16"}),
        "oversample_factor": frozenset({"2"}),
        "sensor_filter_ratio": frozenset({"12.0"}),
    },
}
LAMELLAR_AGGREGATE_EXPECTED_STATES = frozenset(
    (f_a, float(chi_n))
    for f_a, first_chi_n in (
        (0.20, 30),
        (0.25, 21),
        (0.30, 17),
        (0.35, 14),
        (0.40, 13),
        (0.45, 12),
        (0.50, 12),
    )
    for chi_n in range(first_chi_n, 36)
)
assert len(LAMELLAR_AGGREGATE_EXPECTED_STATES) == 133

STRONG_EXTENSION_DIR = RESULTS / "bvk_strong_segregation_extension"
STRONG_EXTENSION_SCHEMA = "bvk-strong-segregation-extension-v1"
STRONG_EXTENSION_STATES = frozenset(
    (f_a, chi_n)
    for f_a in (0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50)
    for chi_n in map(float, range(36, 51))
)
STRONG_EXTENSION_MAIN_MODELS = (
    "ohta_kawasaki",
    "uneyama_doi",
    "liu2019_opf",
    "burp_ti",
    "bvk2",
)
STRONG_EXTENSION_MODELS = (
    "ohta_kawasaki",
    "uneyama_doi",
    "liu2019_opf",
    "burp_ti",
    "bvk2_fixed",
    "bvk2",
)
STRONG_EXTENSION_REPORT_MODELS = STRONG_EXTENSION_MODELS


def rows(relative: str | Path) -> list[dict[str, str]]:
    path = Path(relative)
    if not path.is_absolute():
        path = ROOT / path
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def strong_extension_plot_rows(
    result_dir: str | Path = STRONG_EXTENSION_DIR,
) -> list[dict[str, str]]:
    """Load the accepted 105-state extension and fail closed on inventory drift."""
    result_dir = Path(result_dir).resolve()
    report = json.loads(
        (result_dir / "validation_report.json").read_text(encoding="utf-8")
    )
    if report != {
        "schema": STRONG_EXTENSION_SCHEMA,
        "status": "passed",
        "cohort": "integer_chiN_36_50",
        "state_count": 105,
        "main_row_count": 630,
        "fixed_row_count": 105,
        "models": list(STRONG_EXTENSION_REPORT_MODELS),
    }:
        raise ValueError("strong-segregation extension validation report drifted")

    main = rows(result_dir / "summary.csv")
    fixed = rows(result_dir / "fixed_stiffness_summary.csv")
    selected = [
        row
        for row in main
        if (float(row["f"]), float(row["chiN"])) in STRONG_EXTENSION_STATES
        and row["model"] in {"scft", *STRONG_EXTENSION_MAIN_MODELS}
    ]
    expected_main = {
        (f_a, chi_n, model)
        for f_a, chi_n in STRONG_EXTENSION_STATES
        for model in ("scft", *STRONG_EXTENSION_MAIN_MODELS)
    }
    main_keys = {
        (float(row["f"]), float(row["chiN"]), row["model"])
        for row in selected
    }
    if (
        len(selected) != len(expected_main)
        or main_keys != expected_main
        or any(row["status"] != "accepted" for row in selected)
    ):
        raise ValueError("strong-segregation main-model inventory is incomplete")

    fixed_keys = {
        (float(row["f"]), float(row["chiN"]), row["model"])
        for row in fixed
    }
    expected_fixed = {
        (f_a, chi_n, "bvk2_fixed") for f_a, chi_n in STRONG_EXTENSION_STATES
    }
    if (
        len(fixed) != len(expected_fixed)
        or fixed_keys != expected_fixed
        or any(row["status"] != "accepted" for row in fixed)
    ):
        raise ValueError("strong-segregation fixed-stiffness inventory is incomplete")
    for row in [*selected, *fixed]:
        for metric in ("period_rg", "scft_period_rg", "signed_period_error",
                       "abs_period_error", "profile_rms"):
            if not math.isfinite(float(row[metric])):
                raise ValueError(
                    f"strong-segregation {metric} is nonfinite for "
                    f"model={row['model']}, f={row['f']}, chiN={row['chiN']}"
                )
    return [*selected, *fixed]


def strong_extension_aggregate_rows(
    result_dir: str | Path = STRONG_EXTENSION_DIR,
) -> list[dict[str, str]]:
    """Load the separately reported strong-segregation aggregate."""
    strong_extension_plot_rows(result_dir)
    source_rows = rows(Path(result_dir).resolve() / "aggregate.csv")
    if (
        len(source_rows) != len(STRONG_EXTENSION_REPORT_MODELS)
        or tuple(row["model"] for row in source_rows)
        != STRONG_EXTENSION_REPORT_MODELS
    ):
        raise ValueError("strong-segregation aggregate model inventory drifted")
    for row in source_rows:
        if (
            row["schema"] != STRONG_EXTENSION_SCHEMA
            or row["cohort"] != "integer_chiN_36_50"
            or row["state_count"] != "105"
            or not math.isfinite(float(row["mean_abs_period_error_percent"]))
            or not math.isfinite(float(row["mean_profile_rms"]))
        ):
            raise ValueError("strong-segregation aggregate contains an invalid row")
    by_model = {row["model"]: row for row in source_rows}
    return [by_model[model] for model in STRONG_EXTENSION_MODELS]


def figure6_toc_sources() -> tuple[
    list[dict[str, str]],
    list[dict[str, str]],
    tuple[tuple[float, float, tuple[str, ...]], ...],
]:
    """Load Figure 6 boundaries, SCFT overlay, and topology junctions for the TOC."""
    boundaries = rows(FIGURE6_PLOT_DATA)
    scft = rows(FIGURE6_SCFT_DATA)
    triple_points = []
    for path, f_a_field, chi_n_field, transitions in FIGURE6_TRIPLE_POINT_SOURCES:
        source_rows = rows(path)
        if len(source_rows) != 1:
            raise ValueError(f"Figure 6 requires exactly one row in {path}")
        row = source_rows[0]
        triple_points.append(
            (float(row[f_a_field]), float(row[chi_n_field]), transitions)
        )
    if not boundaries or any(row["status"] != "accepted" for row in boundaries):
        raise ValueError("Figure 6 plotted boundary data must contain accepted rows only")
    if {row["transition"] for row in boundaries} != set(FIGURE6_TRANSITIONS):
        raise ValueError("Figure 6 plotted boundary data has an unexpected transition set")
    return boundaries, scft, tuple(triple_points)


def _boundary_f_a_at_chi(
    boundary_rows: list[dict[str, str]], transition: str, chi_n: float
) -> float:
    points = sorted(
        (float(row["chiN"]), float(row["fA"]))
        for row in boundary_rows
        if row["transition"] == transition
    )
    if len(points) < 2 or not points[0][0] <= chi_n <= points[-1][0]:
        raise ValueError(f"cannot locate {transition} at chiN={chi_n:g}")
    return float(
        np.interp(chi_n, [point[0] for point in points], [point[1] for point in points])
    )


def figure6_region_anchors(
    boundary_rows: list[dict[str, str]],
) -> dict[str, tuple[float, float]]:
    """Place phase labels at midpoints of their displayed AQCE boundaries."""
    specifications = {
        "FCC": (37.0, "FCC/DIS", "FCC/BCC"),
        "BCC": (37.0, "FCC/BCC", "S/C"),
        "CYL": (33.0, "S/C", "G/C"),
        "GYR": (29.0, "G/C", "G/L"),
        "O70": (13.5, "G/O70", "O70/L"),
    }
    anchors = {}
    for label, (chi_n, left_transition, right_transition) in specifications.items():
        left = _boundary_f_a_at_chi(boundary_rows, left_transition, chi_n)
        right = _boundary_f_a_at_chi(boundary_rows, right_transition, chi_n)
        anchors[label] = (0.5 * (left + right), chi_n)
    lam_chi_n = 31.0
    lam_left = _boundary_f_a_at_chi(boundary_rows, "G/L", lam_chi_n)
    anchors["LAM"] = (0.5 * (lam_left + (1.0 - lam_left)), lam_chi_n)
    anchors["DIS"] = (0.5, 9.25)
    return anchors


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_model_hierarchy_source_manifest() -> None:
    sources = (
        (RESULTS / "diblock_kernel_comparison/kernel_samples.csv", "Figure 1 panel (a)"),
        (
            RESULTS / "diblock_weak_response_map/weak_response_map.csv",
            "Figure 1 panels (b) and (c)",
        ),
        (
            RESULTS / "bvk2_aba_comprehensive_validation/benchmark_summary.csv",
            "TOC ABA UD, BURP, and AQCE period errors",
        ),
        (
            RESULTS / "bvk2_aba_fixed_stiffness_validation/summary.csv",
            "TOC ABA QCE period errors",
        ),
        (FIGURE6_PLOT_DATA, "TOC Figure 6 plotted AQCE phase boundaries"),
        (FIGURE6_SCFT_DATA, "TOC Figure 6 SCFT reference boundaries"),
        *(
            (path, "TOC Figure 6 estimated triple-point topology")
            for path, _f_a, _chi_n, _transitions in FIGURE6_TRIPLE_POINT_SOURCES
        ),
    )
    output = RESULTS / "bvk2_model_hierarchy_graphics/source_manifest.csv"
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("source", "role", "data_rows", "sha256"),
            lineterminator="\n",
        )
        writer.writeheader()
        for path, role in sources:
            writer.writerow(
                {
                    "source": path.relative_to(ROOT).as_posix(),
                    "role": role,
                    "data_rows": len(rows(path)),
                    "sha256": _sha256(path),
                }
            )


def save(fig, relative: str | Path, *, dpi: int = 180) -> None:
    path = Path(relative)
    if not path.is_absolute():
        path = ROOT / path
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        path,
        format="svg",
        dpi=dpi,
        metadata=SVG_METADATA,
    )
    plt.close(fig)
    path.write_text(
        "\n".join(line.rstrip() for line in path.read_text(encoding="utf-8").splitlines())
        + "\n",
        encoding="utf-8",
    )


def panel(ax, label: str) -> None:
    ax.text(
        -0.12,
        1.03,
        label,
        transform=ax.transAxes,
        fontsize=10,
        fontweight="bold",
        va="bottom",
    )


def model_label(model: str) -> str:
    return {
        "scft": "SCFT",
        "exact_rpa": "Ideal-chain vertex",
        "ohta_kawasaki": "OK",
        "ok": "OK",
        "uneyama_doi": "UD",
        "liu2019_opf": "OPF",
        "burp_ti": "BURP",
        "bvk1": "nonsmooth precursor",
        "bvk2_fixed": "QCE",
        "bvk2": "AQCE",
        "bvk1_bvk2": "Reduced QCE/AQCE vertex",
    }.get(model, model)


def fixed_stiffness_period_rows(
    source: str | Path = RESULTS / "bvk_fixed_stiffness_period_map_nx256/summary.csv",
    expected_states: set[tuple[float, float]] | None = None,
) -> list[dict[str, str]]:
    """Load the accepted production fixed-stiffness lamellar map."""
    source = Path(source).resolve()
    data = rows(source)
    required = {
        "schema",
        "claim_kind",
        "f",
        "chiN",
        "c2",
        "adaptive",
        "grid_count",
        "oversample_factor",
        "sensor_filter_ratio",
        "discretization_schema",
        "scft_period_rg",
        "period_rg",
        "signed_period_error",
        "abs_period_error",
        "profile_rms",
        "field_gate_pass",
        "composition_gate_pass",
        "morphology_gate_pass",
        "cell_gate_pass",
        "local_minimum_check_pass",
        "stress_orientation_pass",
        "accepted",
        "status",
        "source_root_sha256",
        "source_profile_sha256",
    }
    selected: dict[tuple[float, float], dict[str, str]] = {}
    for row_index, row in enumerate(data, start=2):
        missing = sorted(required - row.keys())
        if missing:
            raise ValueError(
                f"fixed-stiffness source row {row_index} lacks fields {missing}"
            )
        if (
            row["schema"] != "bvk2-fixed-stiffness-period-map-v1"
            or row["claim_kind"] != "stress_free_lamellar_observable"
            or row["status"] != "accepted"
        ):
            raise ValueError(
                "fixed-stiffness source has an invalid schema, claim kind, "
                f"or status at row {row_index}"
            )
        true_gates = (
            "field_gate_pass",
            "composition_gate_pass",
            "morphology_gate_pass",
            "cell_gate_pass",
            "local_minimum_check_pass",
            "stress_orientation_pass",
            "accepted",
        )
        if any(row[field] != "true" for field in true_gates):
            raise ValueError(
                f"fixed-stiffness source row {row_index} fails a required gate"
            )
        if (
            row["adaptive"] != "false"
            or abs(float(row["c2"]) - 0.16) > 1.0e-14
            or row["grid_count"] != "256"
            or row["oversample_factor"] != "2"
            or abs(float(row["sensor_filter_ratio"]) - 12.0) > 1.0e-14
            or row["discretization_schema"]
            != "ud-theta-spectral-fixed-k-physical-sensor-filter-oversampled-v2"
        ):
            raise ValueError(
                f"fixed-stiffness source row {row_index} has wrong production settings"
            )
        for fingerprint in ("source_root_sha256", "source_profile_sha256"):
            if len(row[fingerprint]) != 64 or any(
                char not in "0123456789abcdef" for char in row[fingerprint]
            ):
                raise ValueError(
                    f"fixed-stiffness source row {row_index} has invalid provenance"
                )
        try:
            f_value = float(row["f"])
            chi_n = float(row["chiN"])
            scft_period = float(row["scft_period_rg"])
            model_period = float(row["period_rg"])
            signed_error = float(row["signed_period_error"])
            absolute_error = float(row["abs_period_error"])
            profile_rms = float(row["profile_rms"])
        except (TypeError, ValueError) as error:
            raise ValueError(
                f"fixed-stiffness source row {row_index} has nonnumeric data"
            ) from error
        values = (
            f_value,
            chi_n,
            scft_period,
            model_period,
            signed_error,
            absolute_error,
            profile_rms,
        )
        if not all(math.isfinite(value) for value in values):
            raise ValueError(
                f"fixed-stiffness source row {row_index} has nonfinite data"
            )
        if scft_period <= 0.0 or model_period <= 0.0 or profile_rms < 0.0:
            raise ValueError(
                f"fixed-stiffness source row {row_index} has invalid observables"
            )
        calculated_error = (model_period - scft_period) / scft_period
        if (
            abs(calculated_error - signed_error) > 5.0e-13
            or abs(abs(signed_error) - absolute_error) > 5.0e-13
        ):
            raise ValueError(
                f"fixed-stiffness source row {row_index} has inconsistent errors"
            )
        state = (f_value, chi_n)
        if state in selected:
            raise ValueError(f"duplicate fixed-stiffness state {state}")
        selected[state] = {
            "model": "bvk2_fixed",
            "model_label": "BVK, fixed stiffness",
            "f": row["f"],
            "chiN": row["chiN"],
            "period_rg": row["period_rg"],
            "scft_period_rg": row["scft_period_rg"],
            "signed_period_error": row["signed_period_error"],
            "abs_period_error": row["abs_period_error"],
            "profile_rms": row["profile_rms"],
            "status": "accepted",
            "grid_count": row["grid_count"],
            "source_schema": row["schema"],
        }
    if expected_states is not None and set(selected) != expected_states:
        missing = sorted(expected_states - set(selected))
        unexpected = sorted(set(selected) - expected_states)
        raise ValueError(
            "fixed-stiffness state inventory mismatch: "
            f"missing={missing}, unexpected={unexpected}"
        )
    return [selected[state] for state in sorted(selected)]


def _pdf_panel(pdf: Path, output: str) -> None:
    with tempfile.TemporaryDirectory() as directory:
        prefix = Path(directory) / "panel"
        subprocess.run(
            ["pdftocairo", "-singlefile", "-png", "-r", "300", str(pdf), str(prefix)],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        image = plt.imread(prefix.with_suffix(".png"))
    fig, ax = plt.subplots(
        figsize=(mpltex_acs.width_double_column, 4.8), constrained_layout=True
    )
    ax.imshow(image)
    ax.axis("off")
    save(fig, output, dpi=300)


@mpltex.acs_decorator
def render_homopolymer() -> None:
    source = RESULTS / "homopolymer_urp_vk1_validation/source"
    _pdf_panel(
        source / "accuracy_saddle_approx/accuracy_saddle_approx-eps-converted-to.pdf",
        "results/homopolymer_urp_vk1_validation/urp_vk1_error_map.svg",
    )
    _pdf_panel(
        source / "iw_saddle_approx/iw_saddle_approx-eps-converted-to.pdf",
        "results/homopolymer_urp_vk1_validation/urp_vk1_saddle_profile.svg",
    )


def _kernel_panel(ax) -> None:
    data = rows("results/diblock_kernel_comparison/kernel_samples.csv")
    compositions = (0.1, 0.3, 0.5)
    models = ("uneyama_doi", "bvk1_bvk2", "exact_rpa")
    for model in models:
        for f_a in compositions:
            subset = sorted(
                [
                    row
                    for row in data
                    if row["model"] == model
                    and math.isclose(float(row["f"]), f_a)
                ],
                key=lambda row: float(row["qRg"]),
            )
            ax.plot(
                [float(row["qRg"]) ** 2 for row in subset],
                [float(row["kernel"]) for row in subset],
                color=COLORS[model],
            )
    ax.legend(
        handles=[
            Line2D(
                [0],
                [0],
                color=COLORS[model],
                label={
                    "uneyama_doi": "UD vertex",
                    "bvk1_bvk2": (
                        "QCE/AQCE vertex "
                        r"($\Gamma_\psi$)"
                    ),
                }.get(model, model_label(model)),
            )
            for model in models
        ],
        loc="upper right",
        ncol=1,
        fontsize=6.5,
    )
    annotation_qrg = {0.1: 4.55, 0.3: 4.15, 0.5: 4.65}
    annotation_offset = {0.1: (-2, 7), 0.3: (-2, 7), 0.5: (-2, -14)}
    for f_a in compositions:
        exact_rows = [
            row
            for row in data
            if row["model"] == "exact_rpa"
            and math.isclose(float(row["f"]), f_a)
        ]
        anchor = min(
            exact_rows,
            key=lambda row: abs(float(row["qRg"]) - annotation_qrg[f_a]),
        )
        ax.annotate(
            rf"$f_A={f_a:.2f}$",
            xy=(float(anchor["qRg"]) ** 2, float(anchor["kernel"])),
            xytext=annotation_offset[f_a],
            textcoords="offset points",
            ha="right",
            va="bottom" if annotation_offset[f_a][1] >= 0 else "top",
            color=mpltex.almost_black,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.8, "pad": 0.1},
        )
    ax.set(
        xlabel=r"$k^2R_g^2$",
        ylabel=r"$\Gamma_\psi(k)$",
        yscale="log",
        xlim=(0, 25.5),
        ylim=(1.0e1, 1.0e3),
    )


def _weak_panels(ax_q, ax_spinodal) -> None:
    data = rows("results/diblock_weak_response_map/weak_response_map.csv")
    for ax in (ax_q, ax_spinodal):
        ax.axhline(0, color=mpltex.almost_black, linewidth=0.7, zorder=0)
    for model in ("uneyama_doi", "bvk1_bvk2"):
        subset = sorted([r for r in data if r["model"] == model], key=lambda r: float(r["f"]))
        x = [float(r["f"]) for r in subset]
        ax_q.plot(
            x,
            [100 * float(r["relative_qRg_error"]) for r in subset],
            color=COLORS[model],
            label=model_label(model),
        )
        ax_spinodal.plot(
            x,
            [100 * float(r["relative_spinodal_error"]) for r in subset],
            color=COLORS[model],
        )
    for ax in (ax_q, ax_spinodal):
        ax.set_xlabel(r"$f_A$")
    ax_q.set_ylabel(r"$k^*$ error (\%)")
    ax_spinodal.set_ylabel(r"spinodal error (\%)")


def _weak_error_panel(ax, field: str, ylabel: str) -> None:
    data = rows("results/diblock_weak_response_map/weak_response_map.csv")
    ax.axhline(0, color=COLORS["exact_rpa"], linewidth=0.7, zorder=0)
    for model in ("uneyama_doi", "bvk1_bvk2"):
        subset = sorted(
            [row for row in data if row["model"] == model],
            key=lambda row: float(row["f"]),
        )
        ax.plot(
            [float(row["f"]) for row in subset],
            [100 * float(row[field]) for row in subset],
            color=COLORS[model],
        )
    ax.set(xlabel=r"$f_A$", ylabel=ylabel, xlim=(0.1, 0.5))


def _shared_weak_period_panel(ax) -> None:
    data = rows("results/diblock_weak_response_map/weak_response_map.csv")
    subset = sorted(
        [row for row in data if row["model"] == "bvk1_bvk2"],
        key=lambda row: float(row["f"]),
    )
    ax.plot(
        [float(row["f"]) for row in subset],
        [100 * float(row["relative_period_error"]) for row in subset],
        color=COLORS["bvk1_bvk2"],
        label="UD; QCE/AQCE vertices",
    )
    ax.axhline(
        0,
        color=COLORS["exact_rpa"],
        linewidth=0.7,
        label=model_label("exact_rpa"),
        zorder=0,
    )
    ax.set(
        xlabel=r"$f_A$",
        ylabel=r"weak-period error (\%)",
        xlim=(0.1, 0.5),
    )
    ax.legend(loc="lower right", fontsize=5.7, handlelength=1.6)


@mpltex.acs_decorator
def render_kernel_and_hierarchy(
    *,
    weak_output_dir: Path | None = None,
    include_toc: bool = True,
) -> None:
    def weak_output(default: str) -> str | Path:
        if weak_output_dir is None:
            return default
        return weak_output_dir / Path(default).name

    fig, ax = plt.subplots(
        figsize=(mpltex_acs.width_single_column, 2.35), constrained_layout=True
    )
    _kernel_panel(ax)
    save(
        fig,
        weak_output(
            "results/diblock_kernel_comparison/diblock_kernel_comparison.svg"
        ),
    )

    fig, axes = plt.subplots(
        1, 2, figsize=(mpltex_acs.width_double_column, 2.8), constrained_layout=True
    )
    _weak_panels(*axes)
    axes[0].legend()
    panel(axes[0], "(a)")
    panel(axes[1], "(b)")
    save(
        fig,
        weak_output(
            "results/diblock_weak_response_map/diblock_weak_response_map.svg"
        ),
    )

    fig = plt.figure(
        figsize=(mpltex_acs.width_double_column, 4.2), constrained_layout=True
    )
    grid = fig.add_gridspec(2, 2, width_ratios=(1.12, 1.0))
    kernel = fig.add_subplot(grid[:, 0])
    spinodal = fig.add_subplot(grid[0, 1])
    weak_period = fig.add_subplot(grid[1, 1], sharex=spinodal)
    _kernel_panel(kernel)
    panel(kernel, "(a)")
    _weak_error_panel(
        spinodal,
        "relative_spinodal_error",
        r"spinodal error (\%)",
    )
    panel(spinodal, "(b)")
    _shared_weak_period_panel(weak_period)
    panel(weak_period, "(c)")
    save(
        fig,
        weak_output(
            "results/bvk2_model_hierarchy_graphics/"
            "bvk2_model_hierarchy_weak_response.svg"
        ),
    )

    if not include_toc:
        return

    render_toc()


def _toc_aba_period_errors() -> dict[str, np.ndarray]:
    """Use the same signed SCFT-relative ABA period errors as Figure 7d."""
    summary = rows("results/bvk2_aba_comprehensive_validation/benchmark_summary.csv")
    summary += rows("results/bvk2_aba_fixed_stiffness_validation/summary.csv")
    errors = {}
    for model in ("uneyama_doi", "burp_ti", "bvk2_fixed", "bvk2"):
        subset = sorted(
            (row for row in summary if row["model"] == model),
            key=lambda row: float(row["chiN"]),
        )
        data = np.array([(float(row["chiN"]), 100 * float(row["period_relative_error"]))
                         for row in subset])
        if data.shape != (7, 2) or not np.isfinite(data).all():
            raise ValueError(f"incomplete TOC ABA period errors for {model}")
        if not np.array_equal(data[:, 0], [20, 22, 25, 30, 35, 40, 45]):
            raise ValueError(f"unexpected TOC ABA conditions for {model}")
        if any(row["accepted"].lower() != "true" or float(row["fA"]) != 0.5
               for row in subset):
            raise ValueError(f"invalid TOC ABA state for {model}")
        errors[model] = data
    return errors


@mpltex.acs_decorator
def render_toc() -> None:
    """Architectural transfer TOC, independent of the main-figure renderers."""

    # mpltex supplies the ACS visual style, but Macromolecules requires the TOC
    # graphic itself to be exactly 3.25 x 1.75 in at final size.
    fig = plt.figure(figsize=(3.25, 1.75), constrained_layout=False)
    grid = fig.add_gridspec(
        1,
        2,
        width_ratios=(1.5, 1.0),
        left=0.085,
        right=0.985,
        bottom=0.14,
        top=0.77,
        wspace=0.32,
    )
    phase_ax = fig.add_subplot(grid[0, 0])
    period_ax = fig.add_subplot(grid[0, 1])
    # Continuous, color-segmented curves depict chain architecture only.
    for ax, blocks, label in (
        (phase_ax, (("A", 0.5), ("B", 0.5)), "AB"),
        (period_ax, (("A", 0.25), ("B", 0.5), ("A", 0.25)),
         "ABA"),
    ):
        center = (ax.get_position().x0 + ax.get_position().x1) / 2
        start = 0.0
        for chemistry, fraction in blocks:
            contour = np.linspace(start, start + fraction, 100)
            color = "#3b82a0" if chemistry == "A" else "#df9351"
            fig.add_artist(Line2D(
                center + 0.25 * (contour - 0.5),
                0.905 + 0.022 * np.sin(6 * np.pi * contour),
                transform=fig.transFigure, color=color, linewidth=2.0,
                solid_capstyle="round", gid=f"chain-{label}-{start:g}",
            ))
            start += fraction
        fig.text(center, 0.955, label, ha="center", fontsize=6.0)

    period_data = _toc_aba_period_errors()
    period_ax.axhline(0, color=COLORS["scft"], linewidth=0.55, zorder=1,
                      gid="scft-zero-reference")
    for model, color, marker in (
        ("uneyama_doi", COLORS["uneyama_doi"], "o"),
        ("burp_ti", COLORS["burp_ti"], "^"),
        ("bvk2_fixed", COLORS["bvk2_fixed"], "s"),
        ("bvk2", COLORS["bvk2"], "D"),
    ):
        subset = period_data[model]
        period_ax.plot(
            subset[:, 0],
            subset[:, 1],
            color=color,
            linestyle="-",
            label=model_label(model),
            linewidth=0.85, marker=marker, markersize=2.6,
            markeredgewidth=0.3, clip_on=False, zorder=3,
        )
    period_ax.set(xlim=(20, 45), ylim=(-20, 1), xticks=(20, 30, 45), yticks=(-20, -10, 0))
    period_ax.set_ylabel(r"Period error (\%)", labelpad=0, fontsize=6)
    period_ax.set_xlabel(r"$\chi N$", labelpad=0, fontsize=6)
    period_ax.tick_params(labelsize=5.0, pad=1)
    period_ax.legend(
        loc="center right",
        ncol=1,
        frameon=True,
        facecolor="white",
        framealpha=0.78,
        edgecolor="none",
        fontsize=5.0,
        handlelength=1.3,
        handletextpad=0.25,
        columnspacing=0.55,
        borderaxespad=0.25,
        bbox_to_anchor=(0.98, 0.48),
    )
    fig.text((period_ax.get_position().x0 + period_ax.get_position().x1) / 2,
             0.795, r"$f_A=0.5$", ha="center", fontsize=4.8)

    figure6_rows, scft, triple_points = figure6_toc_sources()
    for (_transition, _side, _curve), group in _groups(
        scft, ("transition", "side", "curve_id")
    ).items():
        scft_curve = sorted(group, key=lambda row: int(row["point_index"]))
        if len(scft_curve) < 2:
            continue
        phase_ax.plot(
            [float(row["fA"]) for row in scft_curve],
            [float(row["chiN"]) for row in scft_curve],
            color=COLORS["scft"],
            linestyle="--",
            linewidth=0.6,
            zorder=1,
        )
    plotted_rows = [dict(row) for row in figure6_rows]
    for f_a, chi_n, transitions in triple_points:
        plotted_rows.extend(
            {
                "transition": transition,
                "fA": f_a,
                "chiN": chi_n,
                "status": "accepted",
                "estimated_display_only": True,
            }
            for transition in transitions
        )
    critical = next(
        (
            (float(row["fA"]), float(row["chiN"]))
            for row in plotted_rows
            if row["transition"] == "S/DIS"
            and abs(float(row["fA"]) - 0.5) <= 1.0e-9
        ),
        None,
    )
    for transition in FIGURE6_TRANSITIONS:
        points = sorted(
            (
                (float(row["fA"]), float(row["chiN"]))
                for row in plotted_rows
                if row["transition"] == transition
            ),
            key=lambda point: point[1],
        )
        if (
            critical is not None
            and transition != "S/DIS"
            and transition not in FIGURE6_CRITICAL_ENDPOINT_EXCLUSIONS
        ):
            points.insert(0, critical)
        for mirrored in (False, True):
            phase_ax.plot(
                [1.0 - point[0] if mirrored else point[0] for point in points],
                [point[1] for point in points],
                color=FIGURE6_COLORS[transition],
                linewidth=0.85,
                zorder=3,
            )
    for label, x, y in (
        ("BCC", 0.173, 43),
        ("C", 0.255, 36),
        ("L", 0.500, 34),
        ("DIS", 0.080, 28.0),
    ):
        phase_ax.text(x, y, label, ha="center", va="center",
                      fontsize=4.4 if label == "BCC" else 5.2)
    gyroid_anchor_chi = 34.0
    gyroid_anchor_f_a = 0.5 * (
        _boundary_f_a_at_chi(figure6_rows, "G/C", gyroid_anchor_chi)
        + _boundary_f_a_at_chi(figure6_rows, "G/L", gyroid_anchor_chi)
    )
    phase_ax.annotate(
        "G", xy=(gyroid_anchor_f_a, gyroid_anchor_chi), xytext=(0.285, 28.5),
        ha="center", va="center", fontsize=5.2,
        arrowprops={"arrowstyle": "-", "color": "#4b5563", "linewidth": 0.45},
    )
    fcc_anchor_chi = 40.0
    fcc_anchor_f_a = 0.5 * (
        _boundary_f_a_at_chi(figure6_rows, "FCC/DIS", fcc_anchor_chi)
        + _boundary_f_a_at_chi(figure6_rows, "FCC/BCC", fcc_anchor_chi)
    )
    phase_ax.annotate(
        "FCC",
        xy=(fcc_anchor_f_a, fcc_anchor_chi),
        xytext=(0.078, 45.0),
        ha="center",
        va="center",
        fontsize=4.8,
        arrowprops={"arrowstyle": "-", "color": "#4b5563", "linewidth": 0.45},
    )
    o70_anchor_chi = 13.0
    o70_anchor_f_a = 0.5 * (
        _boundary_f_a_at_chi(figure6_rows, "G/O70", o70_anchor_chi)
        + _boundary_f_a_at_chi(figure6_rows, "O70/L", o70_anchor_chi)
    )
    phase_ax.annotate(
        r"$O^{70}$",
        xy=(o70_anchor_f_a, o70_anchor_chi),
        xytext=(0.500, 20.5),
        ha="center",
        va="center",
        fontsize=4.8,
        arrowprops={"arrowstyle": "-", "color": "#4b5563", "linewidth": 0.45},
    )
    phase_ax.set(
        xlim=(0.0, 1.0),
        ylim=(10, 50),
        xticks=(0.0, 0.5, 1.0),
        yticks=(10, 50),
        xlabel=r"$f_A$",
        ylabel=r"$\chi N$",
    )
    phase_ax.tick_params(axis="both", labelsize=4.8, pad=1.0, length=2.0)
    phase_ax.xaxis.label.set_size(5.5)
    phase_ax.yaxis.label.set_size(5.5)
    phase_ax.xaxis.labelpad = 1.0
    phase_ax.yaxis.labelpad = 1.0
    phase_ax.legend(
        handles=[
            # Boundary families use their own colors; this neutral key encodes style.
            Line2D([], [], color=mpltex.almost_black, linestyle="-", linewidth=0.85, label="AQCE"),
            Line2D([], [], color=COLORS["scft"], linestyle="--", linewidth=0.6, label="SCFT"),
        ],
        loc="upper center",
        ncol=1,
        frameon=True,
        facecolor="white",
        framealpha=0.78,
        edgecolor="none",
        fontsize=5.0,
        handlelength=1.3,
        handletextpad=0.25,
        columnspacing=0.55,
        borderaxespad=0.25,
        bbox_to_anchor=(0.50, 0.98),
    )

    toc_png = ROOT / "results/bvk2_model_hierarchy_graphics/bvk2_macromolecules_toc.png"
    fig.savefig(
        toc_png,
        format="png",
        dpi=300,
        metadata={"Software": "Matplotlib + mpltex ACS via uv"},
    )
    review_pdf = ROOT / "docs/manuscript/macromolecules/latex_review/figures/toc.pdf"
    review_pdf.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(review_pdf, metadata={"CreationDate": None, "ModDate": None})
    save(fig, "results/bvk2_model_hierarchy_graphics/bvk2_macromolecules_toc.svg", dpi=300)
    _write_model_hierarchy_source_manifest()


def _profile_axes(ax, profile_rows, case_id: str, models: list[str], *, state_text: bool = True) -> None:
    for model in models:
        subset = sorted([r for r in profile_rows if r["case_id"] == case_id and r["model"] == model], key=lambda r: float(r["s"]))
        if subset:
            ax.plot([float(r["s"]) for r in subset], [float(r["phi_a"]) for r in subset], color=COLORS[model], label=model_label(model))
    ax.set(xlabel=r"$x/D$", ylabel=r"$\phi_A$", xlim=(0, 1), ylim=(-0.03, 1.03))
    if state_text:
        summary = next(r for r in rows("results/unified_lamellar_benchmark/summary.csv") if r["case_id"] == case_id)
        ax.text(0.03, 0.93, rf"$f_A={float(summary['f']):.2f}$, $\chi N={float(summary['chiN']):g}$", transform=ax.transAxes, va="top")


def _lamellar_error_axes(period_ax, profile_ax, summary) -> None:
    models = ["ohta_kawasaki", "uneyama_doi", "liu2019_opf", "burp_ti", "bvk2"]
    x = np.arange(len(models))
    period = []
    profile = []
    for model in models:
        subset = [r for r in summary if r["model"] == model]
        period.append(100 * np.mean([float(r["abs_period_error"]) for r in subset]))
        profile.append(np.mean([float(r["profile_rms"]) for r in subset]))
    colors = [COLORS[m] for m in models]
    period_ax.bar(x, period, color=colors, width=0.72)
    profile_ax.bar(x, profile, color=colors, width=0.72)
    for ax in (period_ax, profile_ax):
        ax.set_xticks(x, [model_label(m) for m in models], rotation=30, ha="right")
    period_ax.set_ylabel(r"mean period error (\%)")
    profile_ax.set_ylabel("mean profile RMS")


@mpltex.acs_decorator
def render_lamellar() -> None:
    profile_rows = rows("results/unified_lamellar_benchmark/profiles.csv")
    summary = rows("results/unified_lamellar_benchmark/summary.csv")
    cases = sorted({r["case_id"] for r in summary})
    models = ["scft", "ohta_kawasaki", "uneyama_doi", "liu2019_opf", "burp_ti", "bvk2"]

    fig, axes = plt.subplots(
        2,
        3,
        figsize=(mpltex_acs.width_double_column, 4.6),
        sharex=True,
        sharey=True,
        constrained_layout=True,
    )
    for index, (ax, case) in enumerate(zip(axes.flat, cases)):
        _profile_axes(ax, profile_rows, case, models)
        panel(ax, f"({chr(97 + index)})")
    axes.flat[0].legend(ncol=2, loc="lower center")
    save(fig, "results/unified_lamellar_benchmark/stress_free_lamella_profiles.svg")

    fig, axes = plt.subplots(
        1, 2, figsize=(mpltex_acs.width_double_column, 2.8), constrained_layout=True
    )
    _lamellar_error_axes(*axes, summary)
    panel(axes[0], "(a)")
    panel(axes[1], "(b)")
    save(fig, "results/unified_lamellar_benchmark/stress_free_lamella_errors.svg")


@mpltex.acs_decorator
def render_liu2019_stress_free_period_map(
    result_dir: str | Path = RESULTS / "liu2019_stress_free_period_map",
) -> None:
    """Reproduce Liu et al. Figure 2 and add the QCE/AQCE model hierarchy.

    The numerical campaign owns point selection, convergence gates, and
    profile alignment. The independently recomputed OPF points use the same
    per-state Eq. 18 mapping. Digitized literature data are reserved for the
    numerical validator and are not rendered. This function renders accepted
    calculated rows only and performs no interpolation beyond line segments
    between calculated integer-chiN states.
    """
    result_dir = Path(result_dir).resolve()
    source = result_dir / "summary.csv"
    data = rows(source)
    if not data:
        raise ValueError(f"empty period-map table: {source}")
    accepted = [row for row in data if row["status"] == "accepted"]
    weak_response = rows(RESULTS / "diblock_weak_response_map/weak_response_map.csv")
    ud_literal_spinodals = {
        round(float(row["f"]), 8): float(row["chiN_spinodal"])
        for row in weak_response
        if row["model"] == "uneyama_doi"
    }

    def classified_nonaccepted_ud(row):
        if row["model"] != "uneyama_doi":
            return False
        f_value = round(float(row["f"]), 8)
        spinodal = ud_literal_spinodals.get(f_value)
        if row["status"] not in {"provisional", "rejected"}:
            return False
        reason = row.get("status_reason")
        if reason in {
            "no_bracketed_primitive_period_minimum",
            "no_two_sided_primitive_period_minimum",
        }:
            return (
                row.get("field_gate_pass") == "true"
                and row.get("resolution_gate_pass") == "true"
                and row.get("composition_gate_pass") == "true"
                and row.get("morphology_gate_pass") == "true"
                and row.get("period_local_minimum_check_pass") == "false"
            )
        return (
            spinodal is not None
            and float(row["chiN"]) <= spinodal + 1.0e-10
            and reason
            in {
                "one_or_more_model_gates_failed",
                "solver_failure",
                "outside_ordered_ud_domain",
            }
        )

    classified_nonaccepted = {
        (float(row["f"]), float(row["chiN"]), row["model"])
        for row in data
        if row["status"] != "accepted" and classified_nonaccepted_ud(row)
    }
    summary_models = ["liu2019_opf", "uneyama_doi", "burp_ti", "bvk2"]
    spinodals = {
        0.20: 24.613,
        0.25: 18.172,
        0.30: 14.635,
        0.35: 12.562,
        0.40: 11.344,
        0.45: 10.698,
        0.50: 10.495,
    }
    requested_states = {
        (f_value, float(chi_n))
        for f_value, spinodal in spinodals.items()
        for chi_n in range(math.floor(spinodal) + 1, 36)
    }
    for f_value, chi_n in sorted(requested_states):
        present = {
            row["model"]
            for row in accepted
            if abs(float(row["f"]) - f_value) < 1.0e-12
            and abs(float(row["chiN"]) - chi_n) < 1.0e-12
        }
        missing = set(["scft", *summary_models]) - present
        missing = {
            model
            for model in missing
            if (f_value, chi_n, model) not in classified_nonaccepted
        }
        if missing:
            raise ValueError(
                f"accepted state f={f_value:g}, chiN={chi_n:g} is incomplete: "
                f"{sorted(missing)}"
            )

    actual_provisional = {
        (float(row["f"]), float(row["chiN"]), row["model"])
        for row in data
        if row["status"] != "accepted"
    }
    if actual_provisional != classified_nonaccepted:
        raise ValueError(
            "period-map nonaccepted rows differ from the classified exclusion: "
            f"{sorted(actual_provisional)}"
        )

    fixed_stiffness = fixed_stiffness_period_rows(
        expected_states=requested_states
    )
    extension_rows = strong_extension_plot_rows()
    plot_rows = [*accepted, *fixed_stiffness, *extension_rows]
    period_models = [
        "liu2019_opf",
        "uneyama_doi",
        "burp_ti",
        "bvk2_fixed",
        "bvk2",
    ]
    error_models = [
        "uneyama_doi",
        "liu2019_opf",
        "burp_ti",
        "bvk2_fixed",
        "bvk2",
    ]
    profile_models = [
        "uneyama_doi",
        "liu2019_opf",
        "burp_ti",
        "bvk2_fixed",
        "bvk2",
    ]

    def integer_grid_with_gaps(subset, value_key):
        x_values = []
        y_values = []
        previous = None
        for row in subset:
            chi_n = float(row["chiN"])
            if previous is not None and chi_n - previous > 1.0 + 1.0e-12:
                x_values.append(math.nan)
                y_values.append(math.nan)
            x_values.append(chi_n)
            y_values.append(float(row[value_key]))
            previous = chi_n
        return x_values, y_values

    markers = {
        "liu2019_opf": "o",
        "uneyama_doi": "v",
        "burp_ti": "s",
        "bvk2_fixed": "P",
        "bvk2": "D",
    }
    composition_values = sorted({float(row["f"]) for row in accepted})
    composition_colors = {
        f_value: plt.get_cmap("viridis")(
            index / max(1, len(composition_values) - 1)
        )
        for index, f_value in enumerate(composition_values)
    }

    fig = plt.figure(
        figsize=(mpltex_acs.width_double_column, 5.0), constrained_layout=True
    )
    grid = fig.add_gridspec(2, 1, height_ratios=(1.2, 1.0))
    period_grid = grid[0].subgridspec(1, 2, wspace=0.12)
    error_grid = grid[1].subgridspec(1, 5, wspace=0.06)
    period_axes = [
        fig.add_subplot(period_grid[0, 0]),
        fig.add_subplot(period_grid[0, 1]),
    ]
    error_axes = [
        fig.add_subplot(error_grid[0, index])
        for index in range(5)
    ]

    for index, (period_ax, f_value) in enumerate(
        zip(period_axes, (0.50, 0.30))
    ):
        scft = sorted(
            [
                row
                for row in accepted
                if row["model"] == "scft"
                and abs(float(row["f"]) - f_value) < 1.0e-12
            ],
            key=lambda row: float(row["chiN"]),
        )
        if scft:
            period_ax.plot(
                [float(row["chiN"]) for row in scft],
                [float(row["period_rg"]) for row in scft],
                color=COLORS["scft"],
                linestyle="-",
                linewidth=1.1,
            )
        extension_scft = sorted(
            [
                row
                for row in extension_rows
                if row["model"] == "scft"
                and abs(float(row["f"]) - f_value) < 1.0e-12
            ],
            key=lambda row: float(row["chiN"]),
        )
        if extension_scft:
            period_ax.plot(
                [float(row["chiN"]) for row in extension_scft],
                [float(row["period_rg"]) for row in extension_scft],
                color=COLORS["scft"],
                linestyle="-",
                linewidth=1.1,
            )
        for model in period_models:
            subset = sorted(
                [
                    row
                    for row in plot_rows
                    if row["model"] == model
                    and abs(float(row["f"]) - f_value) < 1.0e-12
                ],
                key=lambda row: float(row["chiN"]),
            )
            if not subset:
                continue
            period_ax.plot(
                [float(row["chiN"]) for row in subset],
                [float(row["period_rg"]) for row in subset],
                linestyle="none",
                marker=markers[model],
                markersize=3.2,
                markerfacecolor=COLORS[model],
                markeredgecolor=COLORS[model],
                markeredgewidth=0.65,
            )
        period_ax.axvspan(35.5, 50.5, color="#eeeeee", zorder=-10)
        period_ax.set(
            xlim=(10, 50.5) if index == 0 else (14.5, 50.5),
            ylim=(3.0, 5.30) if index == 0 else (3.0, 5.15),
            xlabel=r"$\chi N$",
            title=rf"$f_A={f_value:.2f}$",
        )
        if index == 0:
            period_ax.set_ylabel(r"$D_0/R_g$")
        panel(period_ax, f"({chr(97 + index)})")
    legend_models = [*error_models, "scft"]
    model_handles = [
        Line2D(
            [0],
            [0],
            color=COLORS["scft"] if model == "scft" else COLORS[model],
            linestyle="-" if model == "scft" else "none",
            marker=None if model == "scft" else markers[model],
            markersize=3.5,
            label=(
                "OPF"
                if model == "liu2019_opf"
                else "QCE"
                if model == "bvk2_fixed"
                else model_label(model)
            ),
        )
        for model in legend_models
    ]
    fig.legend(
        handles=model_handles,
        ncol=6,
        loc="outside upper center",
        frameon=False,
    )

    finite_errors = [
        float(row["signed_period_error"])
        for row in plot_rows
        if row["model"] in period_models
        and math.isfinite(float(row["signed_period_error"]))
    ]
    error_min = min(-0.02, min(finite_errors) - 0.015)
    error_max = max(0.02, max(finite_errors) + 0.015)
    for index, (ax, model) in enumerate(zip(error_axes, error_models)):
        for f_value in composition_values:
            subset = sorted(
                [
                    row
                    for row in plot_rows
                    if row["model"] == model
                    and abs(float(row["f"]) - f_value) < 1.0e-12
                ],
                key=lambda row: float(row["chiN"]),
            )
            if not subset:
                continue
            chi_values, error_values = integer_grid_with_gaps(
                subset, "signed_period_error"
            )
            ax.plot(
                chi_values,
                error_values,
                color=composition_colors[f_value],
                marker="o",
                markersize=2.4,
                markeredgecolor=mpltex.almost_black,
                markeredgewidth=0.3,
                linewidth=0.8,
                label=rf"$f_A={f_value:.2f}$",
            )
        ax.axhline(0, color=mpltex.almost_black, linestyle="--", linewidth=0.6)
        ax.axvspan(35.5, 50.5, color="#eeeeee", zorder=-10)
        ax.set(
            xlim=(10, 50.5),
            ylim=(error_min, error_max),
            xlabel=r"$\chi N$",
        )
        if index == 0:
            ax.set_ylabel(r"$\epsilon=(D_0-D_{0,\mathrm{SCFT}})/D_{0,\mathrm{SCFT}}$")
        else:
            ax.tick_params(labelleft=False)
        ax.text(
            0.05,
            0.93,
            (
                "OPF"
                if model == "liu2019_opf"
                else "QCE"
                if model == "bvk2_fixed"
                else model_label(model)
            ),
            transform=ax.transAxes,
            va="top",
            ha="left",
        )
        panel(ax, f"({chr(99 + index)})")
    composition_legend = [
        Line2D(
            [0],
            [0],
            color=composition_colors[f_value],
            marker="o",
            markersize=2.8,
            linewidth=0.8,
            label=rf"$f_A={f_value:.2f}$",
        )
        for f_value in composition_values
    ]
    fig.legend(
        handles=composition_legend,
        ncol=7,
        loc="outside lower center",
        frameon=False,
    )
    save(
        fig,
        result_dir / "liu2019_stress_free_period_comparison.svg",
    )

    representative_cases = (
        ("f0.5_chiN12", 0.50, 12.0),
        ("f0.5_chiN30", 0.50, 30.0),
        ("f0.3_chiN17", 0.30, 17.0),
        ("f0.3_chiN30", 0.30, 30.0),
    )
    representative_models = (
        "scft",
        "uneyama_doi",
        "liu2019_opf",
        "burp_ti",
        "bvk2_fixed",
        "bvk2",
    )
    canonical_profiles = rows(result_dir / "profiles.csv")
    fixed_profiles = rows(
        RESULTS / "bvk_fixed_stiffness_period_map_nx256/profiles.csv"
    )
    representative_profiles = [
        *[
            row
            for row in canonical_profiles
            if row["case_id"] in {case[0] for case in representative_cases}
            and row["model"] in representative_models
        ],
        *[
            row
            for row in fixed_profiles
            if row["case_id"] in {case[0] for case in representative_cases}
        ],
    ]
    expected_samples = {
        "scft": 256,
        "uneyama_doi": 256,
        "liu2019_opf": 256,
        "burp_ti": 256,
        "bvk2_fixed": 1024,
        "bvk2": 1024,
    }
    profile_inventory = {
        (row["case_id"], row["model"])
        for row in representative_profiles
    }
    expected_inventory = {
        (case_id, model)
        for case_id, _, _ in representative_cases
        for model in representative_models
    }
    if profile_inventory != expected_inventory:
        raise ValueError(
            "representative profile inventory is incomplete: "
            f"{sorted(expected_inventory - profile_inventory)}"
        )
    for case_id, f_value, chi_n in representative_cases:
        for model in representative_models:
            subset = sorted(
                [
                    row
                    for row in representative_profiles
                    if row["case_id"] == case_id and row["model"] == model
                ],
                key=lambda row: float(row["s"]),
            )
            count = expected_samples[model]
            if len(subset) != count:
                raise ValueError(
                    f"{case_id}/{model} has {len(subset)} profile rows; "
                    f"expected {count}"
                )
            coordinates = [float(row["s"]) for row in subset]
            values = [float(row["phi_a"]) for row in subset]
            if (
                any(not math.isfinite(value) for value in values)
                or any(
                    not math.isclose(coordinate, index / count, abs_tol=1.0e-12)
                    for index, coordinate in enumerate(coordinates)
                )
                or any(
                    not math.isclose(float(row["f"]), f_value, abs_tol=1.0e-12)
                    or not math.isclose(float(row["chiN"]), chi_n, abs_tol=1.0e-12)
                    for row in subset
                )
            ):
                raise ValueError(f"{case_id}/{model} profile data are malformed")

    fig = plt.figure(
        figsize=(mpltex_acs.width_double_column, 5.6), constrained_layout=False
    )
    fig.subplots_adjust(
        left=0.085, right=0.975, top=0.895, bottom=0.155,
        hspace=0.43, wspace=0.25,
    )
    grid = fig.add_gridspec(2, 20, height_ratios=(1.05, 1.0))
    representative_axes = []
    for index in range(4):
        representative_axes.append(
            fig.add_subplot(
                grid[0, 5 * index : 5 * (index + 1)],
                sharey=representative_axes[0] if representative_axes else None,
            )
        )
    rms_axes = []
    for index in range(5):
        rms_axes.append(
            fig.add_subplot(
                grid[1, 4 * index : 4 * (index + 1)],
                sharey=rms_axes[0] if rms_axes else None,
            )
        )

    for index, (ax, (case_id, f_value, chi_n)) in enumerate(
        zip(representative_axes, representative_cases)
    ):
        for model in representative_models:
            subset = sorted(
                [
                    row
                    for row in representative_profiles
                    if row["case_id"] == case_id and row["model"] == model
                ],
                key=lambda row: float(row["s"]),
            )
            ax.plot(
                [float(row["s"]) for row in subset],
                [float(row["phi_a"]) for row in subset],
                color=COLORS[model],
                linewidth=1.0,
            )
        ax.set(
            xlim=(0.0, 1.0),
            ylim=(-0.10, 1.03),
            xlabel=r"$s=x/D_0$",
        )
        ax.set_xticks(
            (0.0, 0.5, 1.0),
            ("0.0", "0.5", "1.0") if index == 0 else ("", "0.5", "1.0"),
        )
        ax.text(
            0.50, 0.97,
            rf"$f_A={f_value:.2f}$, $\chi N={chi_n:g}$",
            transform=ax.transAxes, va="top", ha="center",
        )
        if index == 0:
            ax.set_ylabel(r"$\phi_A(s)$")
        else:
            ax.tick_params(labelleft=False)
        panel(ax, f"({chr(97 + index)})")

    representative_handles = [
        Line2D(
            [0], [0], color=COLORS[model], linewidth=1.1,
            label=(
                "OPF"
                if model == "liu2019_opf"
                else "QCE"
                if model == "bvk2_fixed"
                else model_label(model)
            ),
        )
        for model in representative_models
    ]
    fig.legend(
        handles=representative_handles,
        ncol=6,
        loc="upper center",
        bbox_to_anchor=(0.54, 0.995),
        frameon=False,
    )

    for index, (ax, model) in enumerate(zip(rms_axes, profile_models)):
        for f_value in composition_values:
            subset = sorted(
                [
                    row
                    for row in plot_rows
                    if row["model"] == model
                    and abs(float(row["f"]) - f_value) < 1.0e-12
                ],
                key=lambda row: float(row["chiN"]),
            )
            chi_values, rms_values = integer_grid_with_gaps(subset, "profile_rms")
            ax.plot(
                chi_values,
                rms_values,
                color=composition_colors[f_value],
                marker="o",
                markersize=2.2,
                markeredgecolor=mpltex.almost_black,
                markeredgewidth=0.25,
                linewidth=0.75,
            )
        ax.axvspan(35.5, 50.5, color="#eeeeee", zorder=-10)
        ax.set(xlim=(10, 50.5), xlabel=r"$\chi N$")
        if index == 0:
            ax.set_ylabel(r"profile RMS relative to SCFT")
        else:
            ax.tick_params(labelleft=False)
        label = (
            "OPF"
            if model == "liu2019_opf"
            else "QCE"
            if model == "bvk2_fixed"
            else model_label(model)
        )
        ax.text(
            0.05, 0.93, label, transform=ax.transAxes,
            va="top", ha="left"
        )
        panel(ax, f"({chr(101 + index)})")
    fig.legend(
        handles=composition_legend,
        ncol=7,
        loc="lower center",
        bbox_to_anchor=(0.54, 0.01),
        frameon=False,
    )
    save(
        fig,
        result_dir / "liu2019_stress_free_profile_comparison.svg",
    )


def matched_lamellar_aggregate_subsets(
    accepted: list[dict[str, str]],
) -> dict[str, list[dict[str, str]]]:
    """Return the locked cohort after validating metrics and provenance."""
    accepted_by_model: dict[str, dict[tuple[float, float], dict[str, str]]] = {
        model: {} for model in LAMELLAR_AGGREGATE_MODELS
    }
    for row in accepted:
        required_fields = {
            "schema",
            "claim_kind",
            "model",
            "model_label",
            "calibration_role",
            "f",
            "chiN",
            "status",
            "abs_period_error",
            "profile_rms",
            "period_optimizer",
            "reduced_field_protocol",
            "grid_count",
            "quadrature_order",
            "discretization_schema",
        }
        missing_fields = sorted(required_fields - row.keys())
        if missing_fields:
            raise ValueError(
                f"aggregate-error source lacks required fields {missing_fields}"
            )
        model = row["model"]
        if model not in accepted_by_model:
            continue
        missing_fingerprint_fields = sorted(
            LAMELLAR_AGGREGATE_SOURCE_FINGERPRINTS[model].keys() - row.keys()
        )
        if missing_fingerprint_fields:
            raise ValueError(
                "aggregate-error source lacks fingerprint fields "
                f"{missing_fingerprint_fields} for model={model}"
            )
        if row["status"] != "accepted":
            raise ValueError(
                f"aggregate-error row for model={model} is not accepted"
            )
        if (
            row["schema"] != LAMELLAR_AGGREGATE_SCHEMA
            or row["claim_kind"] != LAMELLAR_AGGREGATE_CLAIM_KIND
        ):
            raise ValueError(
                "aggregate-error source schema/claim fingerprint mismatch for "
                f"model={model}: schema={row['schema']!r}, "
                f"claim_kind={row['claim_kind']!r}"
            )
        state = (float(row["f"]), float(row["chiN"]))
        for metric in ("abs_period_error", "profile_rms"):
            try:
                value = float(row[metric])
            except (TypeError, ValueError) as error:
                raise ValueError(
                    f"aggregate-error {metric} is not numeric for model={model}, "
                    f"f={state[0]:g}, chiN={state[1]:g}"
                ) from error
            if not math.isfinite(value):
                raise ValueError(
                    f"aggregate-error {metric} is not finite for model={model}, "
                    f"f={state[0]:g}, chiN={state[1]:g}"
                )
        if state in accepted_by_model[model]:
            raise ValueError(
                "aggregate-error source contains duplicate accepted row for "
                f"model={model}, f={state[0]:g}, chiN={state[1]:g}"
            )
        accepted_by_model[model][state] = row

    missing_models = [
        model for model, states in accepted_by_model.items() if not states
    ]
    if missing_models:
        raise ValueError(
            f"aggregate-error source lacks accepted rows for {missing_models}"
        )
    for model, expected_fingerprint in (
        LAMELLAR_AGGREGATE_SOURCE_FINGERPRINTS.items()
    ):
        model_rows = accepted_by_model[model].values()
        for field, expected_values in expected_fingerprint.items():
            actual_values = frozenset(row[field] for row in model_rows)
            if actual_values != expected_values:
                raise ValueError(
                    "aggregate-error source fingerprint mismatch for "
                    f"model={model}, field={field}: expected="
                    f"{sorted(expected_values)}, actual={sorted(actual_values)}"
                )
    common_states = set.intersection(
        *(set(states) for states in accepted_by_model.values())
    )
    if common_states != LAMELLAR_AGGREGATE_EXPECTED_STATES:
        missing_states = sorted(LAMELLAR_AGGREGATE_EXPECTED_STATES - common_states)
        unexpected_states = sorted(
            common_states - LAMELLAR_AGGREGATE_EXPECTED_STATES
        )
        raise ValueError(
            "aggregate-error common cohort differs from the locked 133-state "
            f"manifest: missing={missing_states}, unexpected={unexpected_states}"
        )
    return {
        model: [accepted_by_model[model][state] for state in sorted(common_states)]
        for model in LAMELLAR_AGGREGATE_MODELS
    }


def matched_lamellar_aggregate_subsets_with_fixed(
    accepted: list[dict[str, str]],
    fixed_source: str | Path = RESULTS / "bvk_fixed_stiffness_period_map_nx256/summary.csv",
) -> dict[str, list[dict[str, str]]]:
    """Add the fixed-stiffness production control on the locked cohort."""
    production = matched_lamellar_aggregate_subsets(accepted)
    fixed_rows = fixed_stiffness_period_rows(fixed_source)
    fixed_by_state = {
        (float(row["f"]), float(row["chiN"])): row for row in fixed_rows
    }
    if len(fixed_by_state) != len(fixed_rows):
        raise ValueError("fixed-stiffness aggregate source contains duplicate states")
    missing = sorted(LAMELLAR_AGGREGATE_EXPECTED_STATES - set(fixed_by_state))
    if missing:
        raise ValueError(
            "fixed-stiffness aggregate source lacks locked 133-state rows: "
            f"{missing}"
        )
    fixed_subset = [
        fixed_by_state[state] for state in sorted(LAMELLAR_AGGREGATE_EXPECTED_STATES)
    ]
    combined = {
        **production,
        "bvk2_fixed": fixed_subset,
    }
    return {model: combined[model] for model in LAMELLAR_AGGREGATE_PLOT_MODELS}


@mpltex.acs_decorator
def render_lamellar_aggregate_errors(
    result_dir: str | Path = RESULTS / "liu2019_stress_free_period_map",
) -> None:
    """Render both lamellar cohorts in one grouped comparison."""
    result_dir = Path(result_dir).resolve()
    source = result_dir / "summary.csv"
    accepted = [row for row in rows(source) if row["status"] == "accepted"]
    subsets = matched_lamellar_aggregate_subsets_with_fixed(accepted)
    strong = strong_extension_aggregate_rows()
    models = LAMELLAR_AGGREGATE_PLOT_MODELS
    strong_models = tuple(row["model"] for row in strong)
    cohorts = (
        (
            "Liu-domain common cohort",
            models,
            [
                100.0
                * float(
                    np.mean(
                        [float(row["abs_period_error"]) for row in subsets[model]]
                    )
                )
                for model in models
            ],
            [
                100.0
                * float(np.mean([float(row["profile_rms"]) for row in subsets[model]]))
                for model in models
            ],
        ),
        (
            r"Strong segregation, integer $\chi N=36$--$50$",
            strong_models,
            [float(row["mean_abs_period_error_percent"]) for row in strong],
            [100.0 * float(row["mean_profile_rms"]) for row in strong],
        ),
    )

    fig, ax = plt.subplots(
        figsize=(mpltex_acs.width_double_column, 0.47 * mpltex_acs.width_double_column),
        constrained_layout=True,
    )
    bar_width = 0.17
    upper = 1.12 * max(
        value
        for _title, _models, period_errors, profile_errors in cohorts
        for value in (*period_errors, *profile_errors)
    )
    x_positions = np.arange(len(models))
    common = cohorts[0]
    strong_cohort = cohorts[1]
    bar_specs = (
        (-1.65, common[2], ACS_COLORS[0], 0.48),
        (-0.65, strong_cohort[2], ACS_COLORS[0], 0.96),
        (0.65, common[3], ACS_COLORS[1], 0.48),
        (1.65, strong_cohort[3], ACS_COLORS[1], 0.96),
    )
    for offset, values, color, alpha in bar_specs:
        ax.bar(
            x_positions + offset * bar_width,
            values,
            width=bar_width,
            color=color,
            alpha=alpha,
            linewidth=0,
        )
    labels = [
        "QCE" if model == "bvk2_fixed" else model_label(model)
        for model in models
    ]
    ax.set(
        ylim=(0.0, upper),
        ylabel=r"mean error (\%), RMS ($\times 100$)",
    )
    ax.set_xticks(x_positions, labels)
    ax.yaxis.labelpad = 5
    ax.grid(axis="y", color="0.88", linewidth=0.5)
    ax.set_axisbelow(True)
    ax.legend(
        handles=(
            Patch(facecolor=ACS_COLORS[0], label="period"),
            Patch(facecolor=ACS_COLORS[1], label="profile"),
            Patch(
                facecolor="0.35",
                alpha=0.48,
                label="weak to intermediate segregation",
            ),
            Patch(facecolor="0.35", alpha=0.96, label="strong segregation"),
        ),
        frameon=False,
        loc="upper right",
        ncol=2,
    )
    save(fig, result_dir / "aggregate_error_bars.svg")


@mpltex.acs_decorator
def render_cell_stress() -> None:
    result_dir = RESULTS / "lamellar_cell_stress_mechanism"
    marker = result_dir / ".promotion-in-progress"
    if marker.exists():
        raise ValueError("cell-stress data promotion is incomplete")
    manifest_path = result_dir / "generation_manifest.csv"
    if not manifest_path.is_file():
        raise ValueError("cell-stress generation manifest is missing")
    with manifest_path.open(newline="", encoding="utf-8") as handle:
        manifest_rows = list(csv.DictReader(handle))
    if len(manifest_rows) != 1:
        raise ValueError("cell-stress generation manifest must contain one row")
    manifest = manifest_rows[0]
    if (
        manifest.get("schema") != "lamellar-cell-stress-generation-v1"
        or manifest.get("data_schema") != "lamellar-cell-stress-mechanism-v2"
    ):
        raise ValueError("cell-stress generation manifest has wrong schema")
    artifact_fields = {
        "scan.csv": "scan_sha256",
        "summary.csv": "summary_sha256",
        "README.md": "readme_sha256",
    }
    fingerprints: list[str] = []
    for name, field in artifact_fields.items():
        artifact = result_dir / name
        if not artifact.is_file():
            raise ValueError(f"cell-stress generation lacks {name}")
        fingerprint = hashlib.sha256(artifact.read_bytes()).hexdigest()
        if fingerprint != manifest.get(field):
            raise ValueError(
                f"cell-stress generation fingerprint mismatch for {name}"
            )
        fingerprints.append(fingerprint)
    generation_id = hashlib.sha256("|".join(fingerprints).encode()).hexdigest()
    if generation_id != manifest.get("generation_id"):
        raise ValueError("cell-stress generation id does not match artifacts")

    data = rows(result_dir / "scan.csv")
    cases = ["f0.5_chiN20", "f0.35_chiN30"]
    source_models = [
        "uneyama_doi",
        "liu2019_opf",
        "burp_ti",
        "bvk2_fixed",
        "bvk2",
    ]
    plot_models = [
        "uneyama_doi",
        "burp_ti",
        "bvk2_fixed",
        "bvk2",
    ]
    if not data or {row["schema"] for row in data} != {
        "lamellar-cell-stress-mechanism-v2"
    }:
        raise ValueError("cell-stress source must use one v2 schema")
    expected_pairs = {(case, model) for case in cases for model in source_models}
    observed_pairs = {(row["case_id"], row["model"]) for row in data}
    if observed_pairs != expected_pairs:
        raise ValueError("cell-stress source has incomplete model/case inventory")
    row_keys = [
        (row["case_id"], row["model"], round(float(row["period_ratio"]), 12))
        for row in data
    ]
    if len(row_keys) != len(set(row_keys)):
        raise ValueError("cell-stress source contains duplicate branch rows")
    for row in data:
        if row["force_gate_pass"] != "true" or row["morphology_gate_pass"] != "true":
            raise ValueError("cell-stress source contains an unaccepted branch point")
        source = Path(row["source_artifact"])
        if source.is_absolute() or not row["source_fingerprint"]:
            raise ValueError("cell-stress source provenance is invalid")
    fig, axes = plt.subplots(
        2, 2, figsize=(mpltex_acs.width_double_column, 5.0), constrained_layout=True
    )
    for ax in axes.flat:
        ax.grid(False)
    for col, case in enumerate(cases):
        case_rows = [row for row in data if row["case_id"] == case]
        if not case_rows:
            raise ValueError(f"cell-stress source lacks case {case}")
        scft_period = float(case_rows[0]["period_rg"]) / float(
            case_rows[0]["period_ratio"]
        )
        for model in plot_models:
            subset = sorted(
                [r for r in case_rows if r["model"] == model],
                key=lambda r: float(r["period_rg"]),
            )
            if not subset:
                raise ValueError(f"cell-stress source lacks {case} {model}")
            roots = [row for row in subset if row["is_source_root"] == "true"]
            if len(roots) != 1:
                raise ValueError(f"cell-stress source requires one root for {case} {model}")
            x = [float(r["period_rg"]) for r in subset]
            axes[0, col].plot(x, [100.0 * float(r["normalized_excess_density"]) for r in subset], color=COLORS[model], label=model_label(model))
            axes[1, col].plot(x, [float(r["normalized_log_period_stress"]) for r in subset], color=COLORS[model])
            root = roots[0]
            root_x = float(root["period_rg"])
            axes[0, col].plot(
                root_x,
                100.0 * float(root["normalized_excess_density"]),
                marker="o",
                markersize=4.0,
                color=COLORS[model],
                linestyle="none",
            )
            axes[1, col].plot(
                root_x,
                float(root["normalized_log_period_stress"]),
                marker="o",
                markersize=4.0,
                color=COLORS[model],
                linestyle="none",
            )
        for ax in axes[:, col]:
            ax.axvline(
                scft_period,
                color=COLORS["scft"],
                linewidth=0.8,
                linestyle=(0, (4, 3)),
            )
        axes[1, col].set_xlabel(r"$D/R_g$")
        axes[1, col].axhline(0, color=mpltex.almost_black, linewidth=0.7)
        for ax in axes[:, col]:
            ax.yaxis.set_major_formatter(FormatStrFormatter("%.1f"))
        axes[1, col].yaxis.set_major_locator(MultipleLocator(0.5))
        label = r"$f_A=0.50$, $\chi N=20$" if col == 0 else r"$f_A=0.35$, $\chi N=30$"
        axes[0, col].text(0.04, 0.92, label, transform=axes[0, col].transAxes, va="top")
    axes[0, 0].set_ylabel(r"$10^2[F/V-(F/V)_{\min}]/s_m$")
    axes[1, 0].set_ylabel(r"$\sigma/s_m$")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    axes[1, 1].legend(handles, labels, ncol=2, loc="upper left")
    for index, ax in enumerate(axes.flat):
        panel(ax, f"({chr(97 + index)})")
    output = result_dir / "lamellar_cell_stress_mechanism.svg"
    staged = output.with_name(".lamellar_cell_stress_mechanism.svg.tmp")
    save(fig, staged)
    staged.replace(output)




@mpltex.acs_decorator
def render_phase_diagram(
    plot_data_path: str | Path = FIGURE6_PLOT_DATA,
    scft_path: str | Path = FIGURE6_SCFT_DATA,
    output_path: str | Path = FIGURE6_SOURCE_SVG,
    rpa_path: str | Path = FIGURE6_RPA_DATA,
    triple_point_sources: tuple[
        tuple[Path, str, str, tuple[str, ...]], ...
    ] = FIGURE6_TRIPLE_POINT_SOURCES,
) -> None:
    """Render Figure 6 from the canonical assembled publication snapshot.

    Scientific assembly belongs exclusively to
    ``assemble_bvk2_publication_phase_diagram.py``. This renderer consumes its
    stable plot-data table and topology-junction records, so rebuilding the
    other manuscript figures cannot fall back to the legacy four-boundary
    diagram.
    """
    plotted = rows(plot_data_path)
    scft = rows(scft_path)
    # Retain rpa_path for compatibility with older snapshot commands, but
    # Figure 6 no longer reads or displays the RPA stability-limit data.
    if not plotted or any(row["status"] != "accepted" for row in plotted):
        raise ValueError("Figure 6 plot data must contain accepted roots only")
    transitions = {row["transition"] for row in plotted}
    if transitions != set(FIGURE6_TRANSITIONS):
        raise ValueError(
            "Figure 6 plot data has an unexpected transition set: "
            f"{sorted(transitions)}"
        )
    triple_points = []
    for path, f_a_field, chi_n_field, connected_transitions in triple_point_sources:
        source_rows = rows(path)
        if len(source_rows) != 1:
            raise ValueError(f"Figure 6 requires exactly one row in {path}")
        row = source_rows[0]
        triple_points.append(
            (
                float(row[f_a_field]),
                float(row[chi_n_field]),
                connected_transitions,
            )
        )
    plotted_with_junctions = rows_with_estimated_triple_points(
        plotted, tuple(triple_points)
    )
    with mpltex.acs_decorator:
        fig = plt.figure(
            figsize=(mpltex_acs.width_double_column, 4.8),
            constrained_layout=True,
        )
        grid = fig.add_gridspec(2, 2, width_ratios=(4.0, 1.35))
        ax = fig.add_subplot(grid[:, 0])
        fcc_ax = fig.add_subplot(grid[0, 1])
        o70_ax = fig.add_subplot(grid[1, 1])

        scft_curves = []
        for (transition, _side, _curve), group in _groups(
            scft, ("transition", "side", "curve_id")
        ).items():
            curve = sorted(group, key=lambda row: int(row["point_index"]))
            if len(curve) >= 2:
                points = tuple(
                    (float(row["fA"]), float(row["chiN"])) for row in curve
                )
                if transition in {"C/O70", "G/O70", "O70/L", "G/L"}:
                    # Display-only dequantization: retain raw CSV coordinates,
                    # endpoints, mirror symmetry, and shared junctions.
                    smoothed = tuple(smooth_scft_display_curve(
                        list(points), smooth_chi=transition != "G/L"
                    ))
                    displacement = np.abs(np.array(smoothed[::5]) - np.array(points))
                    if (displacement[:, 0].max() > 1e-3
                            or displacement[:, 1].max() > 0.03):
                        raise ValueError("SCFT display smoothing exceeds displacement limits")
                    points = smoothed
                scft_curves.append(points)
        boundary_curves = {}
        for transition in FIGURE6_TRANSITIONS:
            curves, _computed = phase_boundary_geometry(
                plotted_with_junctions,
                transition,
                [],
                critical_endpoint_exclusions=FIGURE6_CRITICAL_ENDPOINT_EXCLUSIONS,
            )
            boundary_curves[transition] = tuple(
                tuple(curve) for curve in curves if len(curve) >= 2
            )

        def plot_layers(target_ax, *, legend_labels: bool) -> None:
            for curve in scft_curves:
                target_ax.plot(
                    [point[0] for point in curve],
                    [point[1] for point in curve],
                    color=COLORS["scft"],
                    linewidth=0.75,
                    linestyle="--",
                    zorder=1,
                )
            for transition in FIGURE6_TRANSITIONS:
                for curve in boundary_curves[transition]:
                    target_ax.plot(
                        [point[0] for point in curve],
                        [point[1] for point in curve],
                        color=FIGURE6_COLORS[transition],
                        linewidth=1.05,
                        zorder=3,
                    )
            for index, (f_a, chi_n, _connected) in enumerate(triple_points):
                target_ax.plot(
                    (f_a, 1.0 - f_a),
                    (chi_n, chi_n),
                    linestyle="none",
                    marker="o",
                    markersize=3.5,
                    markeredgewidth=0,
                    color="#808080",
                    label=(
                        "triple point"
                        if legend_labels and index == 0 else None
                    ),
                    zorder=5,
                )
            if legend_labels:
                target_ax.plot(
                    [], [], color=COLORS["scft"], linewidth=0.75,
                    linestyle="--", label="SCFT"
                )

        for target_ax in (ax, fcc_ax, o70_ax):
            plot_layers(target_ax, legend_labels=target_ax is ax)

        ax.axvline(0.5, color="#9ca3af", linestyle=":", linewidth=0.7)
        ax.set(
            xlabel=r"$f_A$",
            ylabel=r"$\chi N$",
            xlim=(0.0, 1.0),
            ylim=(8.0, 50.0),
        )
        region_anchors = figure6_region_anchors(plotted)
        for label in ("BCC", "CYL", "LAM", "DIS"):
            f_a, chi_n = region_anchors[label]
            ax.text(
                f_a, chi_n, label, ha="center", va="center", fontsize=7.0,
                fontweight="bold", color="#111827", zorder=6,
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.78,
                      "pad": 0.8},
            )
        for label, text_position in (
            ("FCC", (0.105, 39.0)),
            ("GYR", (0.395, 27.5)),
        ):
            ax.annotate(
                label, xy=region_anchors[label], xytext=text_position,
                ha="center", va="center", fontsize=7.0, fontweight="bold",
                color="#111827", zorder=6,
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.78,
                      "pad": 0.8},
                arrowprops={"arrowstyle": "-", "color": "#4b5563",
                            "linewidth": 0.55},
            )

        fcc_ax.set(
            xlim=(0.095, 0.205), ylim=(21.0, 50.0),
            title="FCC pocket", ylabel=r"$\chi N$",
        )
        fcc_ax.annotate(
            "FCC", xy=region_anchors["FCC"], xytext=(0.112, 39.0),
            ha="center", va="center", fontsize=6.5, fontweight="bold", zorder=6,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.78,
                  "pad": 0.6},
            arrowprops={"arrowstyle": "-", "color": "#4b5563",
                        "linewidth": 0.5},
        )
        o70_ax.set(
            xlim=(0.405, 0.475), ylim=(10.6, 15.3),
            title=r"$O^{70}$ pocket", xlabel=r"$f_A$", ylabel=r"$\chi N$",
        )
        o70_ax.annotate(
            r"$O^{70}$", xy=region_anchors["O70"], xytext=(0.45, 14.0),
            ha="center", va="center", fontsize=6.5,
            fontweight="bold", zorder=6,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.78,
                  "pad": 0.6},
            arrowprops={"arrowstyle": "-", "color": "#4b5563",
                        "linewidth": 0.5},
        )
        for inset_ax in (fcc_ax, o70_ax):
            inset_ax.tick_params(axis="both", labelsize=6.0)
            inset_ax.title.set_size(7.0)
            inset_ax.xaxis.label.set_size(7.0)
            inset_ax.yaxis.label.set_size(7.0)
        ax.legend(
            ncol=1,
            loc="lower right",
            fontsize=7.0,
        )
        save(fig, output_path)


def _groups(data, keys):
    grouped = defaultdict(list)
    for row in data:
        grouped[tuple(row[key] for key in keys)].append(row)
    return grouped


def _scft_display_curve(group, transition):
    """Return one calibrated JCP2021 vector branch without curve fitting.

    PostScript integer coordinates occasionally place adjacent vertices on
    the same chiN row.  Averaging those repeated rows removes sampling-order
    duplicates while retaining an explicitly stored critical endpoint at
    f_A=0.5.  No interpolation, smoothing, or extrapolation is applied.
    """
    by_chi = defaultdict(list)
    for row in group:
        chi = float(row["chiN"])
        f_a = float(row["fA"])
        if math.isfinite(chi) and math.isfinite(f_a):
            by_chi[chi].append(f_a)
    if len(by_chi) < 2:
        return np.array([]), np.array([])

    chi = np.array(sorted(by_chi))
    f_a = np.array([float(np.mean(by_chi[value])) for value in chi])
    if any(abs(value - 0.5) < 1.0e-12 for value in by_chi[chi[0]]):
        f_a[0] = 0.5
    return f_a, chi


@mpltex.acs_decorator
def render_crossings() -> None:
    summaries = rows("results/bvk2_boundary_crossing_panels/crossing_summary.csv")
    endpoints = rows("results/bvk2_boundary_crossing_panels/crossing_endpoints.csv")
    fig, axes = plt.subplots(
        2, 2, figsize=(mpltex_acs.width_double_column, 5.2), constrained_layout=True
    )
    for index, (ax, summary) in enumerate(zip(axes.flat, summaries)):
        transition = summary["transition"]
        group = sorted([r for r in endpoints if r["transition"] == transition], key=lambda r: float(r["coordinate"]))
        x = [float(r["coordinate"]) for r in group]
        y = [1000 * float(r["deltaF_phase_a_minus_phase_b"]) for r in group]
        ax.plot(x, y, color=COLORS[transition], marker="o")
        ax.axhline(0, color=mpltex.almost_black, linewidth=0.7)
        ax.axvline(float(summary["root_fA"]), color=COLORS[transition], linestyle="--", linewidth=0.9)
        ax.axvspan(float(summary["bracket_lower"]), float(summary["bracket_upper"]), color=COLORS[transition], alpha=0.10)
        ax.text(0.04, 0.92, rf"{transition}, $\chi N={float(summary['chiN']):g}$", transform=ax.transAxes, va="top")
        ax.set(xlabel=r"$f_A$", ylabel=r"$10^3\Delta F/V$")
        panel(ax, f"({chr(97 + index)})")
    save(fig, "results/bvk2_boundary_crossing_panels/bvk2_boundary_crossings.svg")


def _aba_profile(ax, profile_rows, case, models) -> None:
    for model in models:
        subset = sorted([r for r in profile_rows if r["case_id"] == case and r["model"] == model], key=lambda r: float(r["s"]))
        ax.plot([float(r["s"]) for r in subset], [float(r["phi_a"]) for r in subset], color=COLORS[model], label=model_label(model))
    chi = float(next(r for r in profile_rows if r["case_id"] == case)["chiN"])
    ax.text(0.04, 0.92, rf"$\chi N={chi:g}$", transform=ax.transAxes, va="top")
    ax.set(xlabel=r"$x/D$", ylabel=r"$\phi_A$", xlim=(0, 1), ylim=(-0.03, 1.03))


@mpltex.acs_decorator
def render_aba() -> None:
    profile_rows = rows("results/bvk2_aba_comprehensive_validation/benchmark_profiles.csv")
    summary = rows("results/bvk2_aba_comprehensive_validation/benchmark_summary.csv")
    models = ["scft", "uneyama_doi", "burp_ti", "bvk1", "bvk2"]
    cases = sorted({r["case_id"] for r in summary}, key=lambda c: float(c.split("chiN")[1]))
    fig, axes = plt.subplots(
        2, 4, figsize=(mpltex_acs.width_double_column, 4.5), constrained_layout=True
    )
    for index, case in enumerate(cases):
        _aba_profile(axes.flat[index], profile_rows, case, models)
        panel(axes.flat[index], f"({chr(97 + index)})")
    for model in ("uneyama_doi", "burp_ti", "bvk1", "bvk2"):
        subset = sorted(
            [r for r in summary if r["model"] == model],
            key=lambda r: float(r["chiN"]),
        )
        axes.flat[7].plot(
            [float(r["chiN"]) for r in subset],
            [100 * float(r["period_relative_error"]) for r in subset],
            color=COLORS[model],
            marker="o",
            label=model_label(model),
        )
    axes.flat[7].axhline(0, color=mpltex.almost_black, linewidth=0.7)
    axes.flat[7].set(xlabel=r"$\chi N$", ylabel=r"period error (\%)")
    panel(axes.flat[7], "(h)")
    axes.flat[0].legend(ncol=2, loc="lower center")
    save(fig, "results/bvk2_aba_comprehensive_validation/aba_unified_benchmark.svg")

    stress = rows("results/bvk2_aba_comprehensive_validation/cell_stress_scan.csv")
    stress_cases = sorted({r["case_id"] for r in stress})
    fig, axes = plt.subplots(
        2, 2, figsize=(mpltex_acs.width_double_column, 5.0), constrained_layout=True
    )
    for col, case in enumerate(stress_cases):
        for model in models:
            subset = sorted([r for r in stress if r["case_id"] == case and r["model"] == model], key=lambda r: float(r["period_rg"]))
            x = [float(r["period_rg"]) for r in subset]
            axes[0, col].plot(x, [float(r["delta_energy"]) for r in subset], color=COLORS[model], label=model_label(model))
            axes[1, col].plot(x, [float(r["cell_stress"]) for r in subset], color=COLORS[model])
        axes[0, col].set_ylabel(r"$F/V-F_{\min}/V$")
        axes[1, col].set(xlabel=r"$D/R_g$", ylabel="cell stress")
        axes[1, col].axhline(0, color=mpltex.almost_black, linewidth=0.7)
        axes[0, col].text(0.04, 0.92, rf"$\chi N={float(case.split('chiN')[1]):g}$", transform=axes[0, col].transAxes, va="top")
    axes[0, 0].legend(ncol=2)
    for index, ax in enumerate(axes.flat):
        panel(ax, f"({chr(97 + index)})")
    save(fig, "results/bvk2_aba_comprehensive_validation/aba_cell_stress_mechanism.svg")

    crossing = rows("results/bvk2_aba_comprehensive_validation/ldis_crossing.csv")
    fig, axes = plt.subplots(
        1, 2, figsize=(mpltex_acs.width_double_column, 2.8), constrained_layout=True
    )
    for model in ("scft", "bvk2"):
        subset = sorted([r for r in crossing if r["model"] == model], key=lambda r: float(r["chiN"]))
        x = [float(r["chiN"]) for r in subset]
        axes[0].plot(x, [float(r["delta_free_energy"]) for r in subset], color=COLORS[model], marker="o", label=model_label(model))
        axes[1].plot(x, [float(r["amplitude"]) for r in subset], color=COLORS[model], marker="o")
    axes[0].axhline(0, color=mpltex.almost_black, linewidth=0.7)
    axes[0].set(xlabel=r"$\chi N$", ylabel=r"$F_{LAM}/V-F_{DIS}/V$")
    axes[1].set(xlabel=r"$\chi N$", ylabel="LAM amplitude")
    for index, ax in enumerate(axes):
        panel(ax, f"({chr(97 + index)})")
    axes[0].legend()
    save(fig, "results/bvk2_aba_comprehensive_validation/aba_ldis_crossing.svg")

@mpltex.acs_decorator
def render_convergence_morphology() -> None:
    traces = rows("results/bvk2_si_convergence_morphology/convergence_traces.csv")
    sections = rows("results/bvk2_si_convergence_morphology/morphology_sections.csv")
    phases = ["LAM", "CYL", "BCC", "GYR"]
    fig, axes = plt.subplots(
        2, 4, figsize=(mpltex_acs.width_double_column, 4.5), constrained_layout=True
    )
    for col, phase in enumerate(phases):
        trace = sorted([r for r in traces if r["phase"] == phase and r["record_kind"] == "trace"], key=lambda r: int(r["iteration"]))
        axes[0, col].semilogy([int(r["iteration"]) for r in trace], [max(float(r["abs_gap_to_selected_energy_density"]), 1e-16) for r in trace], color=ACS_COLORS[0], label=r"$|F-F^*|/V$")
        axes[0, col].semilogy([int(r["iteration"]) for r in trace], [max(float(r["projected_theta_r2"]), 1e-16) for r in trace], color=ACS_COLORS[1], linestyle="--", label=r"$R_2$")
        axes[0, col].text(0.05, 0.90, phase, transform=axes[0, col].transAxes, va="top")
        axes[0, col].set_xlabel("iteration")
        if col == 0:
            axes[0, col].set_ylabel("diagnostic magnitude")
            axes[0, col].legend()
        group = [r for r in sections if r["phase"] == phase]
        u = sorted({float(r["u"]) for r in group})
        v = sorted({float(r["v"]) for r in group if r["v"]})
        if not v:
            ordered = sorted(group, key=lambda r: float(r["u"]))
            axes[1, col].plot([float(r["u"]) for r in ordered], [float(r["phi"]) for r in ordered], color=ACS_COLORS[2])
            axes[1, col].set(xlabel=r"$x/D$", ylabel=r"$\phi_A$" if col == 0 else "")
        else:
            ui = {value: index for index, value in enumerate(u)}
            vi = {value: index for index, value in enumerate(v)}
            image = np.full((len(v), len(u)), np.nan)
            for row in group:
                image[vi[float(row["v"])]][ui[float(row["u"])] ] = float(row["phi"])
            axes[1, col].imshow(image, origin="lower", cmap="viridis", vmin=0, vmax=1, aspect="equal", interpolation="nearest")
            axes[1, col].set(xticks=[], yticks=[])
        panel(axes[0, col], f"({chr(97 + col)})")
        panel(axes[1, col], f"({chr(101 + col)})")
    save(fig, "results/bvk2_si_convergence_morphology/bvk2_si_convergence_morphology.svg")


@mpltex.acs_decorator
def render_performance() -> None:
    cost = rows("results/bvk2_scft_cost_benchmark/summary.csv")
    cases = sorted({r["case_id"] for r in cost})
    fig, axes = plt.subplots(
        1, 2, figsize=(mpltex_acs.width_double_column, 2.8), constrained_layout=True
    )
    x = np.arange(len(cases))
    width = 0.36
    for offset, method in ((-width / 2, "bvk2"), (width / 2, "scft")):
        subset = {r["case_id"]: r for r in cost if r["method"] == method}
        axes[0].bar(x + offset, [float(subset[c]["median_solver_wall_s"]) for c in cases], width, color=COLORS[method], label=model_label(method))
        axes[1].bar(x + offset, [float(subset[c]["peak_rss_mib"]) for c in cases], width, color=COLORS[method])
    axes[0].set_ylabel("solver time (s)")
    axes[1].set_ylabel("peak RSS (MiB)")
    for index, ax in enumerate(axes):
        ax.set_xticks(x, [c.replace("_", " ") for c in cases], rotation=25, ha="right")
        panel(ax, f"({chr(97 + index)})")
    axes[0].legend()
    save(fig, "results/bvk2_scft_cost_benchmark/cost_memory_benchmark.svg")

    kernel = rows("results/bvk2_scft_fixedcell_performance/kernel_summary.csv")
    fig, axes = plt.subplots(
        1, 2, figsize=(mpltex_acs.width_double_column, 2.8), constrained_layout=True
    )
    for dimension, ax in (("1", axes[0]), ("3", axes[1])):
        for method in ("bvk2", "scft"):
            subset = sorted([r for r in kernel if r["dimension"] == dimension and r["method"] == method], key=lambda r: int(r["grid_points"]))
            ax.loglog([int(r["grid_points"]) for r in subset], [float(r["median_wall_s"]) for r in subset], marker="o", color=COLORS[method], label=model_label(method))
        ax.set(xlabel="grid points", ylabel="kernel time (s)")
        ax.text(0.05, 0.90, "1D" if dimension == "1" else "3D", transform=ax.transAxes, va="top")
    axes[0].legend()
    panel(axes[0], "(a)")
    panel(axes[1], "(b)")
    save(fig, "results/bvk2_scft_fixedcell_performance/kernel_benchmark.svg")

    continuation = rows("results/bvk2_scft_fixedcell_performance/continuation_summary.csv")
    fig, axes = plt.subplots(
        1, 2, figsize=(mpltex_acs.width_double_column, 2.8), constrained_layout=True
    )
    for architecture, ax in (("AB", axes[0]), ("ABA", axes[1])):
        for method in ("bvk2", "scft"):
            subset = sorted([r for r in continuation if r["architecture"] == architecture and r["method"] == method], key=lambda r: float(r["chiN"]))
            ax.plot([float(r["chiN"]) for r in subset], [float(r["median_total_wall_s"]) for r in subset], marker="o", color=COLORS[method], label=model_label(method))
        ax.set(xlabel=r"$\chi N$", ylabel="fixed-cell time (s)")
        ax.text(0.05, 0.90, architecture, transform=ax.transAxes, va="top")
    axes[0].legend()
    panel(axes[0], "(a)")
    panel(axes[1], "(b)")
    save(fig, "results/bvk2_scft_fixedcell_performance/warm_continuation_benchmark.svg")


def main() -> None:
    render_homopolymer()
    render_kernel_and_hierarchy()
    render_lamellar()
    render_lamellar_aggregate_errors()
    render_cell_stress()
    render_phase_diagram()
    render_crossings()
    render_aba()
    render_convergence_morphology()
    render_performance()
    print("rendered all Macromolecules main-text, SI, and TOC figures with Matplotlib")


if __name__ == "__main__":
    main()
