#!/usr/bin/env python3
"""Validate solver gates, profile provenance, and published-marker agreement.

Run after assembling a Liu et al. period-map campaign. The validator writes
``validation_report.json`` into the selected results directory and exits
nonzero if any model-state row, profile table, OPF workflow invariant, or
published-marker comparison fails.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "liu2019_stress_free_period_map"
MODELS = (
    "scft",
    "liu2019_opf",
    "ohta_kawasaki",
    "uneyama_doi",
    "burp_ti",
    "bvk2",
)
SPINODALS = {
    0.20: 24.613260805718557,
    0.25: 18.17192404431113,
    0.30: 14.634878125696673,
    0.35: 12.562004349286187,
    0.40: 11.343968614645185,
    0.45: 10.697794367172573,
    0.50: 10.494868245363874,
}
SAMPLE_COUNT = 256
BVK2_PROFILE_SAMPLE_COUNT = 1024
SCHEMA = "liu2019-stress-free-period-map-v2"
SCFT_RESIDUAL_TOL = 1.0e-6
SCFT_REQUESTED_SPACING_RG = 0.1
SCFT_CONTOUR_STEP = 0.01
SCFT_REFINED_PROTOCOLS = {
    (0.20, 25.0): (0.05, 0.005),
    (0.20, 26.0): (0.05, 0.005),
    (0.20, 27.0): (0.05, 0.005),
}
SCFT_FIELD_UPDATER = "Anderson(SD(0.1); m=10, warmup=100, αw=0.2)"
SCFT_CELL_UPDATER = "VariableCell(BB(1.0), field_updater)"
REDUCED_FIELD_PROTOCOLS = {
    "liu2019_opf": "mean-constrained Fourier nonlinear conjugate gradient",
    "ohta_kawasaki": "mean-constrained Fourier nonlinear conjugate gradient",
    "uneyama_doi": "mean-constrained fused-spectral L-BFGS with profile continuation",
    "burp_ti": "bounded-KKT Newton-GMRES",
    "bvk2": "oversampled filtered UD-theta L-BFGS with cosine-Newton polish",
}
OPF_STRESS_TOL = 1.0e-4
OPF_FIELD_RMS_TOL = 1.0e-5
OPF_FIELD_MAX_TOL = 1.0e-4
POLYNOMIAL_PROFILE_MODELS = ("liu2019_opf", "ohta_kawasaki")
POLYNOMIAL_SPECTRAL_TAIL_FRACTION = 0.25
# Require the relative RMS amplitude in the highest resolved quarter of the
# native spectrum to be no larger than the accepted OPF/OK field-residual RMS.
# This is a truncation check on the archived grid, not a replacement for a
# matched-grid convergence study.
POLYNOMIAL_SPECTRAL_TAIL_RMS_TOL = OPF_FIELD_RMS_TOL
UD_FIELD_RMS_TOL = 1.0e-6
UD_FIELD_MAX_TOL = 1.0e-5
UD_LITERAL_SPINODALS = {
    0.20: 27.65484092457152,
    0.25: 20.257673359966226,
    0.30: 16.228288742704713,
    0.35: 13.881174016398957,
    0.40: 12.50800212915063,
    0.45: 11.78158581873333,
    0.50: 11.553777520530048,
}
OPF_PUBLISHED_MODEL = "liu2019_opf_mapped_published"
OPF_PUBLISHED_COUNT = 126
OPF_PUBLISHED_RMSE_TOL = 8.0e-3
OPF_PUBLISHED_BIAS_TOL = 5.0e-3
OPF_PUBLISHED_MAXABS_TOL = 2.5e-2
OK_PUBLISHED_MODEL = "ohta_kawasaki_published"
OK_PUBLISHED_COUNT = 77
# Liu et al. used state-specific mapped c3 and c4 values that were not tabulated.
# The independent reproduction therefore uses their published Eqs. 21--22
# regression for c3 and c4.  This envelope tests the resulting Figure 2(c)
# trace without fitting the implementation to digitized marker coordinates.
OK_PUBLISHED_RMSE_TOL = 1.2e-2
OK_PUBLISHED_BIAS_TOL = 5.0e-3
OK_PUBLISHED_MAXABS_TOL = 2.5e-2
BURP_STRESS_TOL = 1.0e-4
BURP_NQUAD = 64
BURP_NEWTON_CAP = 256
BURP_OPTIMIZER = "BoundedKKTNewtonGMRES+AnalyticStressRoot"
BVK2_FIELD_RMS_TOL = 1.0e-7
BVK2_FIELD_MAX_TOL = 5.0e-7
BVK2_STRESS_TOL = 1.0e-8
BVK2_OPTIMIZER = "FilteredUDThetaAnalyticCellStressRoot"
BVK2_OVERSAMPLE = 2
BVK2_SENSOR_FILTER_RATIO = 12.0
BVK2_DISCRETIZATION_SCHEMA = (
    "ud-theta-spectral-adaptive-k-physical-sensor-filter-oversampled-v2"
)


def dominant_nonzero_mode(values: list[float], mean_value: float) -> int:
    count = len(values)
    amplitudes: list[float] = []
    for mode in range(1, min(8, count // 2) + 1):
        coefficient = sum(
            (value - mean_value)
            * complex(
                math.cos(-2.0 * math.pi * mode * index / count),
                math.sin(-2.0 * math.pi * mode * index / count),
            )
            for index, value in enumerate(values)
        )
        amplitudes.append(abs(coefficient))
    return 1 + max(range(len(amplitudes)), key=amplitudes.__getitem__)


def periodic_domain_count(values: list[float], contrast_tolerance: float = 1.0e-4) -> int:
    contrast = max(values) - min(values)
    if not all(math.isfinite(value) for value in values) or contrast <= contrast_tolerance:
        return 0
    threshold = sum(values) / len(values)
    return sum(
        values[index - 1] <= threshold < values[index]
        for index in range(len(values))
    )


def spectral_tail_rms_fraction(
    values: list[float], tail_fraction: float = POLYNOMIAL_SPECTRAL_TAIL_FRACTION
) -> float:
    """Return the relative RMS amplitude in the highest resolved modes.

    The mean mode is removed.  The numerator contains the final
    ``tail_fraction`` of the positive-frequency spectrum through Nyquist, and
    the denominator contains all nonzero positive-frequency modes.  The ratio
    is therefore independent of FFT normalization.
    """
    if len(values) < 8 or not 0.0 < tail_fraction < 1.0:
        return math.inf
    mean_value = sum(values) / len(values)
    mode_powers: list[float] = []
    for mode in range(1, len(values) // 2 + 1):
        coefficient = sum(
            (value - mean_value)
            * complex(
                math.cos(-2.0 * math.pi * mode * index / len(values)),
                math.sin(-2.0 * math.pi * mode * index / len(values)),
            )
            for index, value in enumerate(values)
        )
        mode_powers.append(abs(coefficient) ** 2)
    total_power = sum(mode_powers)
    if not math.isfinite(total_power) or total_power <= 0.0:
        return math.inf
    tail_count = max(1, math.ceil(tail_fraction * len(mode_powers)))
    return math.sqrt(sum(mode_powers[-tail_count:]) / total_power)


def native_profile_inventory_errors(
    rows: list[dict[str, str]], expected_count: int
) -> list[str]:
    """Validate the native-grid row inventory and coordinate provenance."""
    errors: list[str] = []
    if len(rows) != expected_count:
        errors.append(f"has {len(rows)} rows; expected {expected_count}")
    indices: list[int] = []
    parsed_rows: list[tuple[int, float]] = []
    for row_number, row in enumerate(rows, start=1):
        try:
            native_index = int(row["native_index"])
        except (KeyError, TypeError, ValueError):
            errors.append(f"row {row_number} has an invalid native_index")
            continue
        indices.append(native_index)
        try:
            native_grid_count = int(row["native_grid_count"])
        except (KeyError, TypeError, ValueError):
            errors.append(
                f"native_index {native_index} has an invalid native_grid_count"
            )
        else:
            if native_grid_count != expected_count:
                errors.append(
                    f"native_index {native_index} records native_grid_count "
                    f"{native_grid_count}; expected {expected_count}"
                )
        try:
            coordinate = float(row["s"])
        except (KeyError, TypeError, ValueError):
            errors.append(f"native_index {native_index} has an invalid coordinate")
            continue
        if not math.isfinite(coordinate):
            errors.append(f"native_index {native_index} has a nonfinite coordinate")
        else:
            parsed_rows.append((native_index, coordinate))
        for field in ("phi_a", "phi_b"):
            try:
                value = float(row[field])
            except (KeyError, TypeError, ValueError):
                errors.append(f"native_index {native_index} has an invalid {field}")
            else:
                if not math.isfinite(value):
                    errors.append(
                        f"native_index {native_index} has a nonfinite {field}"
                    )

    duplicate_indices = sorted(
        index for index, count in Counter(indices).items() if count > 1
    )
    if duplicate_indices:
        errors.append(f"duplicate native_index values: {duplicate_indices}")
    expected_indices = list(range(1, expected_count + 1))
    if sorted(indices) != expected_indices:
        errors.append(
            "native_index inventory is not the contiguous range "
            f"1:{expected_count}"
        )
    if sorted(indices) == expected_indices:
        for native_index, coordinate in sorted(parsed_rows):
            expected_coordinate = (native_index - 1) / expected_count
            if not math.isclose(
                coordinate, expected_coordinate, rel_tol=0.0, abs_tol=1.0e-12
            ):
                errors.append(
                    f"native_index {native_index} coordinate is {coordinate:.12g}; "
                    f"expected {expected_coordinate:.12g}"
                )
    return errors


def periodic_resample(values: list[float], count: int) -> list[float]:
    output: list[float] = []
    for index in range(count):
        scaled = index * len(values) / count
        left_index = math.floor(scaled)
        fraction = scaled - left_index
        left = values[left_index % len(values)]
        right = values[(left_index + 1) % len(values)]
        output.append((1.0 - fraction) * left + fraction * right)
    return output


def read_csv(results: Path, name: str) -> list[dict[str, str]]:
    with (results / name).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def case_key(row: dict[str, str]) -> tuple[float, float]:
    return round(float(row["f"]), 8), round(float(row["chiN"]), 8)


def is_classified_nonaccepted_ud(row: dict[str, str]) -> bool:
    if row.get("model") != "uneyama_doi":
        return False
    try:
        f_value = round(float(row["f"]), 8)
        chi_n = float(row["chiN"])
        spinodal = UD_LITERAL_SPINODALS[f_value]
    except (KeyError, TypeError, ValueError):
        return False
    if row.get("status") not in {"provisional", "rejected"}:
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
    if chi_n > spinodal + 1.0e-10:
        return False
    return reason in {
        "one_or_more_model_gates_failed",
        "solver_failure",
        "outside_ordered_ud_domain",
    }


def error_metrics(errors: list[float]) -> dict[str, float | int]:
    count = len(errors)
    return {
        "count": count,
        "bias": sum(errors) / count,
        "rmse": math.sqrt(sum(error * error for error in errors) / count),
        "mae": sum(abs(error) for error in errors) / count,
        "maximum_absolute_error": max(abs(error) for error in errors),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=RESULTS,
        help="assembled campaign directory (default: %(default)s)",
    )
    return parser.parse_args()


def main(results: Path = RESULTS) -> None:
    results = results.resolve()
    summary = read_csv(results, "summary.csv")
    profiles = read_csv(results, "profiles.csv")
    native = read_csv(results, "native_profiles.csv")
    digitized = read_csv(results, "liu2019_figure2_digitized.csv")
    expected_states = {
        (f_value, float(chi_n))
        for f_value, spinodal in SPINODALS.items()
        for chi_n in range(math.floor(spinodal) + 1, 36)
    }
    failures: list[str] = []

    wrong_schema = [row for row in summary if row.get("schema") != SCHEMA]
    current_summary = [row for row in summary if row.get("schema") == SCHEMA]
    summary_by_key = {
        (case_key(row), row["model"]): row for row in current_summary
    }
    expected_summary_keys = {
        (state, model) for state in expected_states for model in MODELS
    }
    missing = expected_summary_keys - set(summary_by_key)
    extra = set(summary_by_key) - expected_summary_keys
    if missing:
        failures.append(f"missing summary rows: {len(missing)}")
    if extra:
        failures.append(f"unexpected summary rows: {len(extra)}")
    nonaccepted = [row for row in current_summary if row["status"] != "accepted"]
    classified_nonaccepted = [
        row for row in nonaccepted if is_classified_nonaccepted_ud(row)
    ]
    unexpected_nonaccepted = [
        row for row in nonaccepted if not is_classified_nonaccepted_ud(row)
    ]
    if unexpected_nonaccepted:
        failures.append(
            f"unexpected nonaccepted summary rows: {len(unexpected_nonaccepted)}"
        )
    if wrong_schema:
        failures.append(f"rows with wrong schema: {len(wrong_schema)}")

    for state in sorted(expected_states):
        scft_row = summary_by_key.get((state, "scft"))
        if scft_row is None:
            continue
        residual = float(scft_row["scft_residual_norm"])
        requested_spacing = float(scft_row["scft_requested_spacing_rg"])
        actual_spacing = float(scft_row["scft_actual_spacing_rg"])
        contour_step = float(scft_row["scft_contour_step"])
        expected_spacing, expected_contour_step = SCFT_REFINED_PROTOCOLS.get(
            state, (SCFT_REQUESTED_SPACING_RG, SCFT_CONTOUR_STEP)
        )
        if not residual < SCFT_RESIDUAL_TOL:
            failures.append(f"{state} SCFT residual is {residual:.12g}")
        if not math.isclose(
            requested_spacing, expected_spacing, abs_tol=1.0e-14
        ):
            failures.append(
                f"{state} requested SCFT spacing is {requested_spacing:.12g}"
            )
        if not 0.8 * expected_spacing <= actual_spacing <= 1.2 * expected_spacing:
            failures.append(f"{state} actual SCFT spacing is {actual_spacing:.12g}")
        if not math.isclose(
            contour_step, expected_contour_step, abs_tol=1.0e-14
        ):
            failures.append(f"{state} SCFT contour step is {contour_step:.12g}")
        if scft_row["scft_field_updater"] != SCFT_FIELD_UPDATER:
            failures.append(f"{state} has the wrong SCFT field updater")
        if scft_row["scft_cell_updater"] != SCFT_CELL_UPDATER:
            failures.append(f"{state} has the wrong SCFT cell updater")
        for model, expected_protocol in REDUCED_FIELD_PROTOCOLS.items():
            model_row = summary_by_key.get((state, model))
            if model_row is None:
                continue
            if model_row.get("status") != "accepted":
                continue
            try:
                residual = float(model_row["reduced_field_residual_norm"])
                maxabs = float(model_row["reduced_field_residual_maxabs"])
            except (KeyError, ValueError):
                failures.append(f"{state} {model} lacks a field residual")
                continue
            polynomial_phase_field = model in POLYNOMIAL_PROFILE_MODELS
            residual_tol = (
                UD_FIELD_RMS_TOL
                if model == "uneyama_doi"
                else BVK2_FIELD_RMS_TOL
                if model == "bvk2"
                else OPF_FIELD_RMS_TOL
                if polynomial_phase_field
                else 5.0e-4
            )
            maxabs_tol = (
                UD_FIELD_MAX_TOL
                if model == "uneyama_doi"
                else BVK2_FIELD_MAX_TOL
                if model == "bvk2"
                else OPF_FIELD_MAX_TOL
                if polynomial_phase_field
                else 4.0e-3
            )
            if residual > residual_tol or maxabs > maxabs_tol:
                failures.append(
                    f"{state} {model} field residuals are "
                    f"{residual:.6g}, {maxabs:.6g}"
                )
            if model_row.get("reduced_field_protocol") != expected_protocol:
                failures.append(f"{state} {model} has the wrong field protocol")
            if model == "liu2019_opf":
                if model_row.get("calibration_role") != "per_state_force_stress_mapping":
                    failures.append(f"{state} OPF does not use Eq. 18 mapping")
                if model_row.get("period_optimizer") != "StressRootBisection":
                    failures.append(f"{state} OPF does not use a stress-root cell solve")
                try:
                    mapped_c2 = float(model_row["opf_c2"])
                    mapping_residual = float(model_row["opf_mapping_residual"])
                except (KeyError, ValueError):
                    failures.append(f"{state} OPF lacks mapped coefficients")
                else:
                    if not math.isclose(mapped_c2, -1.0, abs_tol=1.0e-14):
                        failures.append(f"{state} OPF mapped c2 is {mapped_c2:.6g}")
                    if not math.isfinite(mapping_residual):
                        failures.append(f"{state} OPF mapping residual is nonfinite")
                try:
                    stress = abs(float(model_row["stress_norm"]))
                except (KeyError, ValueError):
                    failures.append(f"{state} OPF lacks an isotropic stress")
                    continue
                if stress > OPF_STRESS_TOL:
                    failures.append(f"{state} OPF stress is {stress:.6g}")
            elif model == "ohta_kawasaki":
                if (
                    model_row.get("calibration_role")
                    != "liu2019_ok_asymptotic_quadratic_mapped_c3_c4"
                ):
                    failures.append(f"{state} OK has the wrong coefficient route")
                if model_row.get("period_optimizer") != "StressRootBisection":
                    failures.append(f"{state} OK does not use a stress-root cell solve")
                try:
                    c2 = float(model_row["opf_c2"])
                    c3 = float(model_row["opf_c3"])
                    c4 = float(model_row["opf_c4"])
                    c5 = float(model_row["opf_c5"])
                    c6 = float(model_row["opf_c6"])
                    f_value = float(model_row["f"])
                    chi_n = float(model_row["chiN"])
                except (KeyError, ValueError):
                    failures.append(f"{state} OK lacks its five coefficients")
                else:
                    spinodal = SPINODALS[state[0]]
                    expected_c5 = 1.0 / (4.0 * f_value * (1.0 - f_value))
                    expected_c6 = 3.0 / (
                        4.0 * f_value**2 * (1.0 - f_value) ** 2
                    )
                    k_star = (expected_c6 / expected_c5) ** 0.25
                    expected_c2 = (
                        -chi_n
                        + spinodal
                        - expected_c5 * k_star**2
                        - expected_c6 / k_star**2
                    )
                    for name, actual, expected in (
                        ("c2", c2, expected_c2),
                        ("c5", c5, expected_c5),
                        ("c6", c6, expected_c6),
                    ):
                        if not math.isclose(actual, expected, rel_tol=2.0e-13):
                            failures.append(
                                f"{state} OK {name} is {actual:.12g}; "
                                f"expected {expected:.12g}"
                            )
                    if not math.isfinite(c3) or not math.isfinite(c4):
                        failures.append(f"{state} OK mapped c3 or c4 is nonfinite")
                try:
                    stress = abs(float(model_row["stress_norm"]))
                except (KeyError, ValueError):
                    failures.append(f"{state} OK lacks an isotropic stress")
                    continue
                if stress > OPF_STRESS_TOL:
                    failures.append(f"{state} OK stress is {stress:.6g}")
            elif model == "burp_ti":
                if model_row.get("period_optimizer") != BURP_OPTIMIZER:
                    failures.append(
                        f"{state} BURP-TI does not use the analytic stress root"
                    )
                if model_row.get("period_local_minimum_check_pass") != "true":
                    failures.append(
                        f"{state} BURP-TI lacks the two-sided energy minimum"
                    )
                if model_row.get("period_boundary_limited") != "false":
                    failures.append(f"{state} BURP-TI is boundary limited")
                try:
                    stress = abs(float(model_row["stress_norm"]))
                    nquad = int(float(model_row["quadrature_order"]))
                    newton_cap = int(float(model_row["field_iteration_cap"]))
                    dominant_mode = int(float(model_row["branch_dominant_mode"]))
                except (KeyError, ValueError):
                    failures.append(
                        f"{state} BURP-TI lacks production provenance"
                    )
                else:
                    if stress > BURP_STRESS_TOL:
                        failures.append(f"{state} BURP-TI stress is {stress:.6g}")
                    if nquad != BURP_NQUAD:
                        failures.append(f"{state} BURP-TI nquad is {nquad}")
                    if newton_cap != BURP_NEWTON_CAP:
                        failures.append(
                            f"{state} BURP-TI Newton cap is {newton_cap}"
                        )
                    if dominant_mode != 1:
                        failures.append(
                            f"{state} BURP-TI dominant mode is {dominant_mode}"
                        )
            elif model == "bvk2":
                if model_row.get("period_optimizer") != BVK2_OPTIMIZER:
                    failures.append(
                        f"{state} BVK2 does not use the analytic stress root"
                    )
                if model_row.get("period_local_minimum_check_pass") != "true":
                    failures.append(
                        f"{state} BVK2 lacks the two-sided energy minimum"
                    )
                if model_row.get("period_boundary_limited") != "false":
                    failures.append(f"{state} BVK2 is boundary limited")
                try:
                    stress = abs(float(model_row["stress_norm"]))
                    oversample = int(float(model_row["oversample_factor"]))
                    sensor_filter_ratio = float(model_row["sensor_filter_ratio"])
                except (KeyError, ValueError):
                    failures.append(
                        f"{state} BVK2 lacks discretization or stress provenance"
                    )
                else:
                    if stress > BVK2_STRESS_TOL:
                        failures.append(f"{state} BVK2 stress is {stress:.6g}")
                    if oversample != BVK2_OVERSAMPLE:
                        failures.append(
                            f"{state} BVK2 oversample factor is {oversample}"
                        )
                    if sensor_filter_ratio != BVK2_SENSOR_FILTER_RATIO:
                        failures.append(
                            f"{state} BVK2 sensor cutoff ratio is "
                            f"{sensor_filter_ratio:.12g}"
                        )
                    if (
                        model_row.get("discretization_schema")
                        != BVK2_DISCRETIZATION_SCHEMA
                    ):
                        failures.append(
                            f"{state} BVK2 has the wrong discretization schema"
                        )

        ud_row = summary_by_key.get((state, "uneyama_doi"))
        if ud_row is not None and ud_row.get("status") == "accepted":
            if ud_row.get("period_local_minimum_check_pass") != "true":
                failures.append(f"{state} UD lacks the two-sided energy minimum")
            if ud_row.get("period_boundary_limited") != "false":
                failures.append(f"{state} UD is boundary limited")
            if ud_row.get("field_gate_pass") != "true":
                failures.append(f"{state} UD field solve did not converge")
            if ud_row.get("morphology_gate_pass") != "true":
                failures.append(f"{state} UD lacks a one-domain lamellar profile")
            try:
                int(float(ud_row["branch_dominant_mode"]))
            except (KeyError, ValueError):
                failures.append(f"{state} UD lacks harmonic-mode provenance")

    published_opf = [
        row for row in digitized if row["model"] == OPF_PUBLISHED_MODEL
    ]
    opf_period_errors: list[float] = []
    opf_period_errors_by_f: dict[float, list[float]] = defaultdict(list)
    for published in published_opf:
        state = case_key(published)
        calculated = summary_by_key.get((state, "liu2019_opf"))
        if calculated is None:
            failures.append(f"{state} lacks a calculated OPF marker")
            continue
        error = float(calculated["signed_period_error"]) - float(
            published["signed_period_error"]
        )
        opf_period_errors.append(error)
        opf_period_errors_by_f[state[0]].append(error)
    if len(published_opf) != OPF_PUBLISHED_COUNT:
        failures.append(
            f"digitized OPF marker count is {len(published_opf)}; "
            f"expected {OPF_PUBLISHED_COUNT}"
        )
    published_metrics = error_metrics(opf_period_errors)
    if published_metrics["rmse"] > OPF_PUBLISHED_RMSE_TOL:
        failures.append(
            f"OPF published-marker RMSE is {published_metrics['rmse']:.6g}"
        )
    if abs(published_metrics["bias"]) > OPF_PUBLISHED_BIAS_TOL:
        failures.append(
            f"OPF published-marker bias is {published_metrics['bias']:.6g}"
        )
    if (
        published_metrics["maximum_absolute_error"]
        > OPF_PUBLISHED_MAXABS_TOL
    ):
        failures.append(
            "OPF published-marker maximum absolute error is "
            f"{published_metrics['maximum_absolute_error']:.6g}"
        )

    published_ok = [
        row for row in digitized if row["model"] == OK_PUBLISHED_MODEL
    ]
    ok_period_errors: list[float] = []
    ok_period_errors_by_f: dict[float, list[float]] = defaultdict(list)
    for published in published_ok:
        state = case_key(published)
        calculated = summary_by_key.get((state, "ohta_kawasaki"))
        if calculated is None:
            failures.append(f"{state} lacks a calculated OK marker")
            continue
        error = float(calculated["signed_period_error"]) - float(
            published["signed_period_error"]
        )
        ok_period_errors.append(error)
        ok_period_errors_by_f[state[0]].append(error)
    if len(published_ok) != OK_PUBLISHED_COUNT:
        failures.append(
            f"digitized OK marker count is {len(published_ok)}; "
            f"expected {OK_PUBLISHED_COUNT}"
        )
    ok_published_metrics = error_metrics(ok_period_errors)
    if ok_published_metrics["rmse"] > OK_PUBLISHED_RMSE_TOL:
        failures.append(
            f"OK published-marker RMSE is {ok_published_metrics['rmse']:.6g}"
        )
    if abs(ok_published_metrics["bias"]) > OK_PUBLISHED_BIAS_TOL:
        failures.append(
            f"OK published-marker bias is {ok_published_metrics['bias']:.6g}"
        )
    if (
        ok_published_metrics["maximum_absolute_error"]
        > OK_PUBLISHED_MAXABS_TOL
    ):
        failures.append(
            "OK published-marker maximum absolute error is "
            f"{ok_published_metrics['maximum_absolute_error']:.6g}"
        )

    profiles_by_key: dict[tuple[tuple[float, float], str], list[dict[str, str]]] = (
        defaultdict(list)
    )
    for row in profiles:
        profiles_by_key[(case_key(row), row["model"])].append(row)
    native_by_key: dict[tuple[tuple[float, float], str], list[dict[str, str]]] = (
        defaultdict(list)
    )
    for row in native:
        native_by_key[(case_key(row), row["model"])].append(row)
    validated_native_by_key: dict[
        tuple[tuple[float, float], str], list[dict[str, str]]
    ] = {}
    native_inventory_checks = 0
    for key, summary_row in summary_by_key.items():
        if summary_row.get("status") != "accepted":
            continue
        try:
            expected_native = int(float(summary_row["grid_count"]))
        except (KeyError, TypeError, ValueError):
            failures.append(f"{key[0]} {key[1]} lacks a valid native grid count")
            continue
        native_rows = native_by_key.get(key, [])
        inventory_errors = native_profile_inventory_errors(
            native_rows, expected_native
        )
        native_inventory_checks += 1
        if inventory_errors:
            for error in inventory_errors:
                failures.append(f"{key[0]} {key[1]} native profile: {error}")
            continue
        validated_native_by_key[key] = sorted(
            native_rows, key=lambda row: int(row["native_index"])
        )

    rms_checks = 0
    polynomial_native_profile_checks = Counter()
    polynomial_spectral_tail_maxima = {
        model: {"value": 0.0, "state": None}
        for model in POLYNOMIAL_PROFILE_MODELS
    }
    for state in sorted(expected_states):
        reference_rows = sorted(
            profiles_by_key.get((state, "scft"), []), key=lambda row: float(row["s"])
        )
        if len(reference_rows) != SAMPLE_COUNT:
            failures.append(f"{state} SCFT aligned profile has {len(reference_rows)} rows")
            continue
        reference_native_rows = validated_native_by_key.get((state, "scft"), [])
        if not reference_native_rows:
            continue
        reference_native = [
            float(row["phi_a"]) for row in reference_native_rows
        ]
        for model in MODELS:
            key = (state, model)
            summary_row = summary_by_key.get(key)
            if summary_row is None or summary_row.get("status") != "accepted":
                continue
            model_rows = sorted(
                profiles_by_key.get(key, []), key=lambda row: float(row["s"])
            )
            expected_aligned_count = (
                BVK2_PROFILE_SAMPLE_COUNT if model == "bvk2" else SAMPLE_COUNT
            )
            if len(model_rows) != expected_aligned_count:
                failures.append(
                    f"{state} {model} aligned profile has {len(model_rows)} rows"
                )
                continue
            values = [float(row["phi_a"]) for row in model_rows]
            if not all(math.isfinite(value) for value in values):
                failures.append(f"{state} {model} aligned profile is nonfinite")
                continue
            target_values = periodic_resample(
                reference_native, expected_aligned_count
            )
            rms = math.sqrt(sum(
                (value - target) ** 2
                for value, target in zip(values, target_values)
            ) / expected_aligned_count)
            reported = float(summary_row["profile_rms"])
            if not math.isclose(rms, reported, rel_tol=2.0e-11, abs_tol=2.0e-11):
                failures.append(
                    f"{state} {model} RMS mismatch: {rms:.12g} != {reported:.12g}"
                )
            rms_checks += 1

            ordered_native_rows = validated_native_by_key.get(key)
            if ordered_native_rows is None:
                continue
            if model in POLYNOMIAL_PROFILE_MODELS:
                native_values = [
                    float(row["phi_a"]) for row in ordered_native_rows
                ]
                domain_count = periodic_domain_count(native_values)
                dominant_mode = dominant_nonzero_mode(native_values, state[0])
                spectral_tail = spectral_tail_rms_fraction(native_values)
                polynomial_native_profile_checks[model] += 1
                if domain_count != 1:
                    failures.append(
                        f"{state} {model} native profile has "
                        f"{domain_count} periodic domains"
                    )
                if dominant_mode != 1:
                    failures.append(
                        f"{state} {model} native dominant mode is "
                        f"{dominant_mode}"
                    )
                if spectral_tail > POLYNOMIAL_SPECTRAL_TAIL_RMS_TOL:
                    failures.append(
                        f"{state} {model} native spectral-tail RMS fraction is "
                        f"{spectral_tail:.6g}"
                    )
                if spectral_tail > polynomial_spectral_tail_maxima[model]["value"]:
                    polynomial_spectral_tail_maxima[model] = {
                        "value": spectral_tail,
                        "state": {"f": state[0], "chiN": state[1]},
                    }
            elif model in {"burp_ti", "uneyama_doi", "bvk2"}:
                native_values = [
                    float(row["phi_a"]) for row in ordered_native_rows
                ]
                if model in {"burp_ti", "bvk2"}:
                    mode = dominant_nonzero_mode(native_values, state[0])
                    if mode != 1:
                        failures.append(
                            f"{state} {model} native dominant mode is {mode}"
                        )
                else:
                    if periodic_domain_count(native_values) != 1:
                        failures.append(
                            f"{state} UD native profile is not one periodic domain"
                        )

    report = {
        "schema": "liu2019-stress-free-period-map-validation-v3",
        "status": "passed" if not failures else "failed",
        "expected_states": len(expected_states),
        "expected_models_per_state": len(MODELS),
        "summary_rows": len(summary),
        "aligned_profile_rows": len(profiles),
        "native_profile_rows": len(native),
        "native_profile_inventory_checks": native_inventory_checks,
        "status_counts": dict(Counter(row["status"] for row in summary)),
        "classified_nonaccepted_ud_rows": len(classified_nonaccepted),
        "bvk2_oversample_factor": BVK2_OVERSAMPLE,
        "bvk2_aligned_profile_sample_count": BVK2_PROFILE_SAMPLE_COUNT,
        "bvk2_sensor_filter_ratio": BVK2_SENSOR_FILTER_RATIO,
        "bvk2_discretization_schema": BVK2_DISCRETIZATION_SCHEMA,
        "polynomial_native_profile_gate": {
            "models": list(POLYNOMIAL_PROFILE_MODELS),
            "required_periodic_domains": 1,
            "required_dominant_mode": 1,
            "spectral_tail_fraction": POLYNOMIAL_SPECTRAL_TAIL_FRACTION,
            "spectral_tail_rms_tolerance": POLYNOMIAL_SPECTRAL_TAIL_RMS_TOL,
            "checks_by_model": dict(polynomial_native_profile_checks),
            "maximum_spectral_tail_rms_fraction_by_model": (
                polynomial_spectral_tail_maxima
            ),
        },
        "opf_published_marker_metrics": published_metrics,
        "opf_published_marker_metrics_by_composition": {
            f"{f_value:.2f}": error_metrics(errors)
            for f_value, errors in sorted(opf_period_errors_by_f.items())
        },
        "ok_published_marker_metrics": ok_published_metrics,
        "ok_published_marker_metrics_by_composition": {
            f"{f_value:.2f}": error_metrics(errors)
            for f_value, errors in sorted(ok_period_errors_by_f.items())
        },
        "rms_checks": rms_checks,
        "failures": failures,
    }
    output = results / "validation_report.json"
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main(parse_args().results_dir)
