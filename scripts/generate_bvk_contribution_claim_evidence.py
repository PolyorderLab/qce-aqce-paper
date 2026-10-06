#!/usr/bin/env python3
"""Freeze the numerical evidence used to frame the QCE/AQCE contribution.

This script reads only canonical, accepted CSV ledgers.  It validates cohort
membership, row uniqueness, acceptance gates, and the separation between
production and diagnostic evidence before writing one deterministic JSON
artifact.  Use ``--check`` to compare a fresh reconstruction with the checked
artifact without modifying it.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Iterable


PROJECT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = PROJECT / "results/bvk_contribution_claim_evidence/claim_evidence.json"

PRODUCTION_MAP = Path("results/liu2019_stress_free_period_map/summary.csv")
FIXED_MAP = Path("results/bvk_fixed_stiffness_period_map_nx256/summary.csv")
CELL_STRESS = Path("results/lamellar_cell_stress_mechanism/summary.csv")
COEFFICIENT_DIAGNOSTIC = Path("results/bvk2_global_c2_refit/stage2_rows.csv")
COEFFICIENT_RANKING = Path("results/bvk2_global_c2_refit/candidate_ranking.csv")
ABA_ADAPTIVE = Path("results/bvk2_aba_comprehensive_validation/benchmark_summary.csv")
ABA_FIXED = Path("results/bvk2_aba_fixed_stiffness_validation/summary.csv")
STRONG_MAP = Path("results/bvk_strong_segregation_extension/summary.csv")
STRONG_FIXED = Path("results/bvk_strong_segregation_extension/fixed_stiffness_summary.csv")
STRONG_AGGREGATE = Path("results/bvk_strong_segregation_extension/aggregate.csv")
STRONG_VALIDATION = Path("results/bvk_strong_segregation_extension/validation_report.json")
PHASE_BOUNDARIES = Path("results/bvk2_publication_phase_diagram/accepted_phase_boundaries.csv")
SCFT_PHASE_REFERENCE = Path("results/bvk2_publication_phase_diagram/scft_reference_boundaries.csv")

PRODUCTION_MODELS = (
    "ohta_kawasaki",
    "uneyama_doi",
    "liu2019_opf",
    "burp_ti",
    "bvk2",
)
STRONG_MODELS = (
    "scft",
    "ohta_kawasaki",
    "uneyama_doi",
    "liu2019_opf",
    "burp_ti",
    "bvk2",
)
COMMON_EXPECTED_STATES = frozenset(
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
FULL_EXPECTED_STATES = frozenset(
    (f_a, float(chi_n))
    for f_a, first_chi_n in (
        (0.20, 25),
        (0.25, 19),
        (0.30, 15),
        (0.35, 13),
        (0.40, 12),
        (0.45, 11),
        (0.50, 11),
    )
    for chi_n in range(first_chi_n, 36)
)
assert len(FULL_EXPECTED_STATES) == 146
ABA_EXPECTED_CHIN = (20.0, 22.0, 25.0, 30.0, 35.0, 40.0, 45.0)
STRONG_EXPECTED_STATES = frozenset(
    (f_a, float(chi_n))
    for f_a in (0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50)
    for chi_n in range(36, 51)
)
assert len(STRONG_EXPECTED_STATES) == 105
PHASE_TRANSITION_COUNTS = {
    "C/O70": 9,
    "FCC/BCC": 6,
    "FCC/DIS": 6,
    "G/C": 22,
    "G/L": 24,
    "G/O70": 6,
    "O70/L": 16,
    "S/C": 27,
    "S/DIS": 33,
}
SCFT_PHASE_TRANSITION_COUNTS = {
    "C/G": 200,
    "C/O70": 102,
    "DIS/S": 102,
    "DIS/S_cp": 200,
    "G/L": 200,
    "G/O70": 126,
    "O70/L": 102,
    "S/C": 102,
    "S_cp/S": 200,
}
ORDERED_PHASES = frozenset({"FCC", "BCC", "CYL", "GYR", "O70", "LAM"})
TRAINING_STATES_IN_PRODUCTION_MAP = frozenset(
    {"f0.5_chiN15", "f0.5_chiN25", "f0.35_chiN30"}
)


def read_rows(path: Path) -> list[dict[str, str]]:
    with (PROJECT / path).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def read_json(path: Path) -> dict[str, object]:
    payload = json.loads((PROJECT / path).read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"{path.as_posix()} is not a JSON object")
    return payload


def source_record(path: Path, evidence_class: str, role: str) -> dict[str, object]:
    payload = (PROJECT / path).read_bytes()
    report = {
        "path": path.as_posix(),
        "sha256": hashlib.sha256(payload).hexdigest(),
        "evidence_class": evidence_class,
        "role": role,
    }
    return report


def require_unique(
    rows: Iterable[dict[str, str]], fields: tuple[str, ...], label: str
) -> None:
    keys = [tuple(row[field] for field in fields) for row in rows]
    if len(keys) != len(set(keys)):
        raise ValueError(f"{label} contains duplicate {fields} rows")


def require_columns(
    rows: list[dict[str, str]], required: frozenset[str], label: str
) -> None:
    if not rows:
        raise ValueError(f"{label} is empty")
    missing = required - set(rows[0])
    if missing:
        raise ValueError(f"{label} lacks required columns: {sorted(missing)}")
    for index, row in enumerate(rows[1:], start=2):
        if set(row) != set(rows[0]):
            raise ValueError(f"{label} row {index} has inconsistent columns")


def require_value(
    rows: Iterable[dict[str, str]], field: str, expected: str, label: str
) -> None:
    actual = {row[field] for row in rows}
    if actual != {expected}:
        raise ValueError(
            f"{label} {field} contract failed: expected {expected!r}, found {sorted(actual)!r}"
        )


def require_values(
    rows: Iterable[dict[str, str]], field: str, expected: frozenset[str], label: str
) -> None:
    actual = {row[field] for row in rows}
    if actual != set(expected):
        raise ValueError(
            f"{label} {field} contract failed: expected {sorted(expected)!r}, "
            f"found {sorted(actual)!r}"
        )


def require_current_or_legacy_label(
    rows: Iterable[dict[str, str]],
    field: str,
    current: str,
    legacy: frozenset[str],
    label: str,
) -> None:
    """Accept one uniform current or frozen historical presentation label."""
    actual = {row[field] for row in rows}
    accepted = {current, *legacy}
    if len(actual) != 1 or not actual <= accepted:
        raise ValueError(
            f"{label} {field} contract failed: expected one of "
            f"{sorted(accepted)!r}, found {sorted(actual)!r}"
        )


def require_state_manifest(
    rows: Iterable[dict[str, str]],
    expected: frozenset[tuple[float, float]],
    label: str,
    *,
    f_field: str = "f",
) -> None:
    actual = {state(row, f_field) for row in rows}
    if actual != expected:
        missing = sorted(expected - actual)
        extra = sorted(actual - expected)
        raise ValueError(
            f"{label} state manifest failed: missing={missing[:3]!r}, extra={extra[:3]!r}"
        )


def require_pair_cardinality(
    rows: Iterable[dict[str, str]],
    fields: tuple[str, ...],
    expected_count: int,
    label: str,
) -> None:
    counts: dict[tuple[str, ...], int] = {}
    for row in rows:
        key = tuple(row[field] for field in fields)
        counts[key] = counts.get(key, 0) + 1
    wrong = {key: count for key, count in counts.items() if count != expected_count}
    if wrong:
        raise ValueError(
            f"{label} pair cardinality failed for {fields}: {list(wrong.items())[:3]!r}"
        )


def require_true(rows: Iterable[dict[str, str]], fields: tuple[str, ...], label: str) -> None:
    for row in rows:
        failed = [field for field in fields if row[field] != "true"]
        if failed:
            identity = row.get("case_id", "unknown")
            raise ValueError(f"{label} row {identity} fails gates: {failed}")


def require_false(rows: Iterable[dict[str, str]], fields: tuple[str, ...], label: str) -> None:
    for row in rows:
        failed = [field for field in fields if row[field] != "false"]
        if failed:
            identity = row.get("case_id", "unknown")
            raise ValueError(f"{label} row {identity} fails negative gates: {failed}")


def require_finite(rows: Iterable[dict[str, str]], fields: tuple[str, ...], label: str) -> None:
    for index, row in enumerate(rows, start=1):
        for field in fields:
            try:
                value = float(row[field])
            except (TypeError, ValueError) as exc:
                raise ValueError(f"{label} row {index} has invalid {field}") from exc
            if not math.isfinite(value):
                raise ValueError(f"{label} row {index} has nonfinite {field}")


def mean(rows: list[dict[str, str]], field: str, *, absolute: bool = False) -> float:
    values = [float(row[field]) for row in rows]
    if absolute:
        values = [abs(value) for value in values]
    if not values or not all(math.isfinite(value) for value in values):
        raise ValueError(f"cannot aggregate finite {field} values")
    return math.fsum(values) / len(values)


def state(row: dict[str, str], f_field: str = "f") -> tuple[float, float]:
    return float(row[f_field]), float(row["chiN"])


def state_list(states: Iterable[tuple[float, float]]) -> list[dict[str, float]]:
    return [{"f_A": f_a, "chiN": chi_n} for f_a, chi_n in sorted(states)]


def claim(
    claim_id: str,
    evidence_class: str,
    cohort_id: str,
    metric: str,
    value: float | int,
    unit: str,
    state_count: int,
    sources: list[str],
    *,
    model: str | None = None,
    qualifier: str | None = None,
) -> dict[str, object]:
    result: dict[str, object] = {
        "claim_id": claim_id,
        "evidence_class": evidence_class,
        "cohort_id": cohort_id,
        "metric": metric,
        "value": value,
        "unit": unit,
        "state_count": state_count,
        "sources": sources,
    }
    if model is not None:
        result["model"] = model
    if qualifier is not None:
        result["qualifier"] = qualifier
    return result


def expected_claim_ids() -> frozenset[str]:
    ids: set[str] = set()
    for model in (*PRODUCTION_MODELS, "bvk2_fixed"):
        ids.update({f"common133.{model}.period", f"common133.{model}.profile"})
    for cohort_id in ("held_out_143", "full_146_bvk2", "full_146_fixed"):
        ids.update({f"{cohort_id}.period", f"{cohort_id}.profile"})
    for c2 in (0.0, 0.04, 0.06, 0.08, 0.10, 0.12, 0.16, 0.20, 0.24):
        label = format(c2, ".2f")
        ids.update(
            {
                f"coefficient_diagnostic.full146.c2_{label}.period",
                f"coefficient_diagnostic.full146.c2_{label}.profile",
            }
        )
    for c2 in (0.0, 0.04, 0.16):
        label = format(c2, ".2f")
        ids.update(
            {
                f"coefficient_diagnostic.common133.c2_{label}.period",
                f"coefficient_diagnostic.common133.c2_{label}.profile",
            }
        )
    ids.update(
        {
            "high_chiN.state_count",
            "high_chiN.adaptive_period_wins",
            "high_chiN.fixed.period",
            "high_chiN.adaptive.period",
            "high_chiN_diagnostic.adaptive_period_wins",
            "high_chiN_diagnostic.fixed.period",
            "high_chiN_diagnostic.adaptive.period",
            "aba7.bvk2_fixed.period",
            "aba7.bvk2_fixed.profile",
            "aba7.bvk2.period",
            "aba7.bvk2.profile",
            "aba7.adaptive_period_wins",
        }
    )
    ids.update(
        f"high_chiN.fA_{f_a:.2f}.adaptive_period_wins"
        for f_a in (0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50)
    )
    for case_id in ("f0.35_chiN30", "f0.5_chiN20"):
        for model in ("bvk2", "bvk2_fixed"):
            ids.update(
                {
                    f"cell_stress.{case_id}.{model}.period_ratio",
                    f"cell_stress.{case_id}.{model}.period_error",
                }
            )
    for model in (
        "ohta_kawasaki",
        "uneyama_doi",
        "liu2019_opf",
        "burp_ti",
        "bvk2_fixed",
        "bvk2",
    ):
        ids.update({f"strong105.{model}.period", f"strong105.{model}.profile"})
    ids.update(
        {
            "strong105.state_count",
            "phase149.accepted_root_count",
            "phase149.ordered_phase_count",
            "phase149.o70_boundary_root_count",
            "phase149.o70_boundary_family_count",
            "phase149.fcc_boundary_root_count",
            "phase149.fcc_boundary_family_count",
        }
    )
    return frozenset(ids)


def require_claim_integrity(
    claims: list[dict[str, object]],
    cohorts: dict[str, dict[str, object]],
    sources: list[dict[str, object]],
) -> None:
    required_fields = frozenset(
        {
            "claim_id",
            "evidence_class",
            "cohort_id",
            "metric",
            "value",
            "unit",
            "state_count",
            "sources",
        }
    )
    allowed_fields = required_fields | {"model", "qualifier"}
    permitted_evidence_classes = frozenset(
        {
            "production",
            "production_control",
            "production_held_out",
            "oracle_diagnostic",
            "production_transfer",
            "production_mechanism",
            "production_strong_segregation",
            "production_phase_boundary",
            "topology_inference",
        }
    )
    metric_units = {
        "mean_absolute_period_error": "percent",
        "mean_aligned_profile_rms": "fraction",
        "matched_state_count": "count",
        "adaptive_lower_absolute_period_error_count": "count",
        "period_ratio_to_scft": "ratio",
        "absolute_period_error": "percent",
        "accepted_boundary_root_count": "count",
        "ordered_phase_count": "count",
        "o70_boundary_root_count": "count",
        "o70_boundary_family_count": "count",
        "fcc_boundary_root_count": "count",
        "fcc_boundary_family_count": "count",
    }
    claim_ids = [str(row["claim_id"]) for row in claims]
    if len(claim_ids) != len(set(claim_ids)):
        raise ValueError("claim table contains duplicate claim IDs")
    expected = expected_claim_ids()
    if set(claim_ids) != expected:
        missing = sorted(expected - set(claim_ids))
        extra = sorted(set(claim_ids) - expected)
        raise ValueError(f"claim-ID contract failed: missing={missing!r}, extra={extra!r}")
    source_paths = {str(row["path"]) for row in sources}
    expected_source_paths = {
        path.as_posix()
        for path in (
            PRODUCTION_MAP,
            FIXED_MAP,
            CELL_STRESS,
            COEFFICIENT_DIAGNOSTIC,
            COEFFICIENT_RANKING,
            ABA_ADAPTIVE,
            ABA_FIXED,
            STRONG_MAP,
            STRONG_FIXED,
            STRONG_AGGREGATE,
            STRONG_VALIDATION,
            PHASE_BOUNDARIES,
            SCFT_PHASE_REFERENCE,
        )
    }
    if source_paths != expected_source_paths:
        raise ValueError("source-table path contract failed")
    source_classes = {str(row["path"]): str(row["evidence_class"]) for row in sources}
    cohort_class_compatibility = {
        "common_133": {"production", "production_control", "oracle_diagnostic"},
        "held_out_143": {"production_held_out"},
        "full_146_bvk2": {"production"},
        "full_146_fixed": {"production_control"},
        "coefficient_diagnostic_full_146": {"oracle_diagnostic"},
        "high_chiN_30_35": {"production", "production_control"},
        "high_chiN_30_35_diagnostic": {"oracle_diagnostic"},
        "aba_7": {"production_transfer", "production_control"},
        "cell_stress_2_state": {"production_mechanism"},
        "strong_segregation_105": {
            "production_strong_segregation",
            "production_control",
        },
        "phase_boundary_149": {
            "production_phase_boundary",
            "topology_inference",
        },
    }
    if set(cohorts) != set(cohort_class_compatibility):
        raise ValueError("cohort evidence-class compatibility contract failed")
    for cohort_id, expected_classes in cohort_class_compatibility.items():
        actual_classes = set(cohorts[cohort_id].get("permitted_evidence_classes", []))
        if actual_classes != expected_classes:
            raise ValueError(
                f"cohort {cohort_id} evidence-class compatibility metadata failed"
            )
    source_class_compatibility = {
        "production": {"production"},
        "production_held_out": {"production"},
        "oracle_diagnostic": {"oracle_diagnostic"},
        "production_transfer": {"production_transfer"},
        "production_mechanism": {"production_mechanism"},
        "production_control": {"production", "production_control", "production_transfer"},
        "production_strong_segregation": {"production_strong_segregation"},
        "production_phase_boundary": {"production_phase_boundary"},
        "topology_inference": {"production_phase_boundary", "literature_reference"},
    }
    for row in claims:
        missing_fields = required_fields - set(row)
        extra_fields = set(row) - allowed_fields
        if missing_fields or extra_fields:
            raise ValueError(
                f"claim {row.get('claim_id', '<missing>')} field schema failed: "
                f"missing={sorted(missing_fields)!r}, extra={sorted(extra_fields)!r}"
            )
        evidence_class = str(row["evidence_class"])
        if evidence_class not in permitted_evidence_classes:
            raise ValueError(
                f"claim {row['claim_id']} has prohibited evidence class {evidence_class}"
            )
        if row["cohort_id"] not in cohorts:
            raise ValueError(
                f"claim {row['claim_id']} references unknown cohort {row['cohort_id']}"
            )
        unknown_sources = set(row["sources"]) - source_paths
        if unknown_sources:
            raise ValueError(
                f"claim {row['claim_id']} references unknown sources {sorted(unknown_sources)!r}"
            )
        cohort_id = str(row["cohort_id"])
        if evidence_class not in cohort_class_compatibility[cohort_id]:
            raise ValueError(
                f"claim {row['claim_id']} evidence class {evidence_class} is incompatible "
                f"with cohort {cohort_id}"
            )
        incompatible_sources = {
            source
            for source in row["sources"]
            if source_classes[str(source)] not in source_class_compatibility[evidence_class]
        }
        if incompatible_sources:
            raise ValueError(
                f"claim {row['claim_id']} evidence class {evidence_class} is incompatible "
                f"with sources {sorted(incompatible_sources)!r}"
            )
        metric = str(row["metric"])
        if metric not in metric_units:
            raise ValueError(f"claim {row['claim_id']} has prohibited metric {metric}")
        if row["unit"] != metric_units[metric]:
            raise ValueError(f"claim {row['claim_id']} metric/unit contract failed")
        value = row["value"]
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise ValueError(f"claim {row['claim_id']} has a nonnumeric value")
        if not math.isfinite(float(value)):
            raise ValueError(f"claim {row['claim_id']} has a nonfinite value")
        if not 0 < int(row["state_count"]):
            raise ValueError(f"claim {row['claim_id']} has invalid state_count")
        cohort = cohorts[str(row["cohort_id"])]
        cohort_limit = int(
            cohort.get("state_count", cohort.get("state_count_per_candidate", 0))
        )
        if int(row["state_count"]) > cohort_limit:
            raise ValueError(
                f"claim {row['claim_id']} exceeds cohort cardinality {cohort_limit}"
            )


def build_evidence() -> dict[str, object]:
    production = read_rows(PRODUCTION_MAP)
    fixed = read_rows(FIXED_MAP)
    cell_stress = read_rows(CELL_STRESS)
    diagnostic = read_rows(COEFFICIENT_DIAGNOSTIC)
    diagnostic_ranking = read_rows(COEFFICIENT_RANKING)
    aba_adaptive_all = read_rows(ABA_ADAPTIVE)
    aba_fixed = read_rows(ABA_FIXED)
    strong = read_rows(STRONG_MAP)
    strong_fixed = read_rows(STRONG_FIXED)
    strong_aggregate = read_rows(STRONG_AGGREGATE)
    strong_validation = read_json(STRONG_VALIDATION)
    phase_boundaries = read_rows(PHASE_BOUNDARIES)
    scft_phase_reference = read_rows(SCFT_PHASE_REFERENCE)

    require_columns(
        production,
        frozenset(
            {
                "schema",
                "claim_kind",
                "case_id",
                "f",
                "chiN",
                "model",
                "model_label",
                "status",
                "calibration_role",
                "abs_period_error",
                "profile_rms",
                "period_optimizer",
                "period_local_minimum_check_pass",
                "period_boundary_limited",
                "field_gate_pass",
                "cell_gate_pass",
                "resolution_gate_pass",
                "composition_gate_pass",
                "morphology_gate_pass",
                "reduced_field_protocol",
                "grid_count",
                "c2",
                "oversample_factor",
                "sensor_filter_ratio",
                "discretization_schema",
            }
        ),
        "production map",
    )
    require_value(
        production,
        "schema",
        "liu2019-stress-free-period-map-v2",
        "production map",
    )
    require_value(
        production,
        "claim_kind",
        "stress_free_lamellar_observable",
        "production map",
    )
    production_models = frozenset({"scft", *PRODUCTION_MODELS})
    require_values(production, "model", production_models, "production map")
    require_pair_cardinality(production, ("f", "chiN"), 6, "production map")
    production_identity = {
        "scft": {
            "model_label": "SCFT",
            "calibration_role": "reference",
            "period_optimizer": "Polyorder.cell_solve!",
            "reduced_field_protocol": "not_applicable",
        },
        "ohta_kawasaki": {
            "model_label": "Ohta-Kawasaki (recomputed)",
            "calibration_role": "liu2019_ok_asymptotic_quadratic_mapped_c3_c4",
            "period_optimizer": "StressRootBisection",
            "reduced_field_protocol": "mean-constrained Fourier nonlinear conjugate gradient",
            "grid_count": "128",
        },
        "uneyama_doi": {
            "model_label": "Uneyama-Doi",
            "calibration_role": "not_applicable",
            "accepted_period_optimizer": "BootstrapLocal",
            "reduced_field_protocol": "mean-constrained fused-spectral L-BFGS with profile continuation",
            "accepted_grid_count": "128",
        },
        "liu2019_opf": {
            "model_label": "OPF (force/stress mapped)",
            "calibration_role": "per_state_force_stress_mapping",
            "period_optimizer": "StressRootBisection",
            "reduced_field_protocol": "mean-constrained Fourier nonlinear conjugate gradient",
            "grid_count": "128",
        },
        "burp_ti": {
            "model_label": "BURP-TI",
            "calibration_role": "not_applicable",
            "period_optimizer": "BoundedKKTNewtonGMRES+AnalyticStressRoot",
            "reduced_field_protocol": "bounded-KKT Newton-GMRES",
            "grid_count": "128",
        },
    }
    for model, identity in production_identity.items():
        model_rows = [row for row in production if row["model"] == model]
        if len(model_rows) != 146:
            raise ValueError(f"production comparator {model} does not contain 146 rows")
        require_state_manifest(
            model_rows,
            FULL_EXPECTED_STATES,
            f"production comparator {model}",
        )
        for field in (
            "model_label",
            "calibration_role",
            "period_optimizer",
            "reduced_field_protocol",
            "grid_count",
        ):
            expected = identity.get(field)
            if expected is not None:
                require_value(model_rows, field, expected, f"production comparator {model}")
        if model != "uneyama_doi":
            require_value(model_rows, "status", "accepted", f"production comparator {model}")
        require_value(model_rows, "c2", "", f"production comparator {model}")
    ud_rows = [row for row in production if row["model"] == "uneyama_doi"]
    require_values(
        ud_rows,
        "status",
        frozenset({"accepted", "provisional", "rejected"}),
        "production comparator uneyama_doi",
    )
    ud_accepted = [row for row in ud_rows if row["status"] == "accepted"]
    if len(ud_accepted) != 133:
        raise ValueError("production comparator uneyama_doi accepted cohort is not 133 rows")
    require_value(
        ud_accepted,
        "period_optimizer",
        str(production_identity["uneyama_doi"]["accepted_period_optimizer"]),
        "accepted production comparator uneyama_doi",
    )
    require_value(
        ud_accepted,
        "grid_count",
        str(production_identity["uneyama_doi"]["accepted_grid_count"]),
        "accepted production comparator uneyama_doi",
    )
    bvk2_source_rows = [row for row in production if row["model"] == "bvk2"]
    if len(bvk2_source_rows) != 146:
        raise ValueError("production BVK2 source does not contain 146 rows")
    require_state_manifest(bvk2_source_rows, FULL_EXPECTED_STATES, "production BVK2 source")
    require_current_or_legacy_label(
        bvk2_source_rows,
        "model_label",
        "AQCE",
        frozenset({"APRM", "BVK2"}),
        "production AQCE source",
    )
    require_value(bvk2_source_rows, "status", "accepted", "production BVK2 source")
    require_values(
        bvk2_source_rows,
        "calibration_role",
        frozenset({"training", "holdout", "frozen_c2_evaluation"}),
        "production BVK2 source",
    )
    require_value(bvk2_source_rows, "c2", "0.16", "production BVK2 source")
    require_value(bvk2_source_rows, "grid_count", "256", "production BVK2 source")
    require_value(bvk2_source_rows, "oversample_factor", "2", "production BVK2 source")
    require_value(bvk2_source_rows, "sensor_filter_ratio", "12.0", "production BVK2 source")
    require_value(
        bvk2_source_rows,
        "discretization_schema",
        "ud-theta-spectral-adaptive-k-physical-sensor-filter-oversampled-v2",
        "production BVK2 source",
    )
    bvk2_role_cases = {
        role: {row["case_id"] for row in bvk2_source_rows if row["calibration_role"] == role}
        for role in ("training", "holdout", "frozen_c2_evaluation")
    }
    if bvk2_role_cases["training"] != TRAINING_STATES_IN_PRODUCTION_MAP:
        raise ValueError("production BVK2 training-role membership contract failed")
    if bvk2_role_cases["holdout"] != {
        "f0.5_chiN12",
        "f0.5_chiN20",
        "f0.5_chiN30",
    }:
        raise ValueError("production BVK2 holdout-role membership contract failed")
    if len(bvk2_role_cases["frozen_c2_evaluation"]) != 140:
        raise ValueError("production BVK2 evaluation-role cardinality contract failed")
    require_value(
        bvk2_source_rows,
        "period_optimizer",
        "FilteredUDThetaAnalyticCellStressRoot",
        "production BVK2 source",
    )
    require_value(
        bvk2_source_rows,
        "reduced_field_protocol",
        "oversampled filtered UD-theta L-BFGS with cosine-Newton polish",
        "production BVK2 source",
    )

    require_columns(
        fixed,
        frozenset(
            {
                "schema",
                "case_id",
                "f",
                "chiN",
                "model",
                "model_label",
                "status",
                "calibration_role",
                "adaptive",
                "c2",
                "grid_count",
                "oversample_factor",
                "sensor_filter_ratio",
                "discretization_schema",
                "abs_period_error",
                "profile_rms",
                "field_gate_pass",
                "composition_gate_pass",
                "morphology_gate_pass",
                "cell_gate_pass",
                "local_minimum_check_pass",
                "stress_orientation_pass",
                "accepted",
            }
        ),
        "fixed-stiffness map",
    )
    require_value(fixed, "schema", "bvk2-fixed-stiffness-period-map-v1", "fixed-stiffness map")
    require_value(fixed, "model", "bvk2_fixed", "fixed-stiffness map")
    require_current_or_legacy_label(
        fixed,
        "model_label",
        "QCE",
        frozenset({"PRM", "BVK, fixed stiffness"}),
        "fixed-stiffness map",
    )
    require_value(fixed, "status", "accepted", "fixed-stiffness map")
    require_value(fixed, "calibration_role", "fixed_stiffness_production_control", "fixed-stiffness map")
    require_value(fixed, "adaptive", "false", "fixed-stiffness map")
    require_value(fixed, "c2", "0.16", "fixed-stiffness map")
    require_value(fixed, "grid_count", "256", "fixed-stiffness map")
    require_value(fixed, "oversample_factor", "2", "fixed-stiffness map")
    require_value(fixed, "sensor_filter_ratio", "12.0", "fixed-stiffness map")
    require_value(
        fixed,
        "discretization_schema",
        "ud-theta-spectral-fixed-k-physical-sensor-filter-oversampled-v2",
        "fixed-stiffness map",
    )
    require_state_manifest(fixed, FULL_EXPECTED_STATES, "fixed-stiffness map")
    require_pair_cardinality(fixed, ("f", "chiN"), 1, "fixed-stiffness map")

    require_columns(
        diagnostic,
        frozenset(
            {
                "case_id",
                "split",
                "f",
                "chiN",
                "c2",
                "nx",
                "initial_grid_nx",
                "period_abs_error",
                "profile_rms",
                "solver_schema",
                "period_solver_schema",
                "evidence_tier",
                "stationarity_pass",
                "phase_identity_pass",
                "cell_stationarity_pass",
                "selection_eligible",
                "model_converged",
                "period_local_minimum_check_pass",
                "period_boundary_limited",
                "fingerprint",
            }
        ),
        "coefficient diagnostic",
    )
    require_values(
        diagnostic,
        "split",
        frozenset({"train", "holdout", "unassigned"}),
        "coefficient diagnostic",
    )
    require_value(diagnostic, "nx", "64", "coefficient diagnostic")
    require_value(diagnostic, "initial_grid_nx", "64", "coefficient diagnostic")
    require_value(
        diagnostic,
        "solver_schema",
        "mean_logit_cosine_lbfgs_residual_newton_v2",
        "coefficient diagnostic",
    )
    require_value(
        diagnostic,
        "period_solver_schema",
        "analytic_log_period_cell_stress_root_v4",
        "coefficient diagnostic",
    )
    require_value(
        diagnostic,
        "evidence_tier",
        "stationary_calibration_candidate",
        "coefficient diagnostic",
    )
    require_value(diagnostic, "model_converged", "true", "coefficient diagnostic")
    expected_c2_strings = frozenset(
        {"0.0", "0.04", "0.06", "0.08", "0.1", "0.12", "0.16", "0.2", "0.24"}
    )
    require_values(diagnostic, "c2", expected_c2_strings, "coefficient diagnostic")
    if any(
        not row["fingerprint"].startswith("model=BVK2|")
        or "|schema=bvk2-rpa-eta-v1|" not in row["fingerprint"]
        or "|law_mode=adaptive|" not in row["fingerprint"]
        or "|discretization=central_eta_forward_edge_exact_v1|" not in row["fingerprint"]
        or f"|c2={float(row['c2']):.17g}|" not in row["fingerprint"]
        for row in diagnostic
    ):
        raise ValueError("coefficient diagnostic model/adaptive/c2 fingerprint contract failed")
    require_pair_cardinality(diagnostic, ("f", "chiN"), 9, "coefficient diagnostic")

    require_columns(
        diagnostic_ranking,
        frozenset(
            {
                "c2",
                "case_count",
                "eligible",
                "stratified_balanced_loss",
                "period_mae",
                "profile_rms_mean",
                "complete",
                "rank",
            }
        ),
        "coefficient ranking",
    )

    require_columns(
        aba_adaptive_all,
        frozenset(
            {
                "schema",
                "case_id",
                "fA",
                "chiN",
                "model",
                "period_relative_error",
                "profile_rms",
                "grid",
                "field_converged",
                "accepted",
                "c2",
            }
        ),
        "ABA adaptive benchmark",
    )
    require_value(
        aba_adaptive_all,
        "schema",
        "bvk2-aba-comprehensive-validation-v2",
        "ABA adaptive benchmark",
    )
    require_values(
        aba_adaptive_all,
        "model",
        frozenset({"scft", "uneyama_doi", "burp_ti", "bvk1", "bvk2"}),
        "ABA adaptive benchmark",
    )
    require_value(aba_adaptive_all, "accepted", "true", "ABA adaptive benchmark")
    require_pair_cardinality(aba_adaptive_all, ("case_id",), 5, "ABA adaptive benchmark")
    require_pair_cardinality(aba_adaptive_all, ("model",), 7, "ABA adaptive benchmark")

    require_columns(
        aba_fixed,
        frozenset(
            {
                "schema",
                "case_id",
                "fA",
                "chiN",
                "model",
                "status",
                "adaptive",
                "c2",
                "nx",
                "period_relative_error",
                "profile_rms",
                "field_gate_pass",
                "composition_gate_pass",
                "morphology_gate_pass",
                "cell_gate_pass",
                "local_minimum_check_pass",
                "accepted",
            }
        ),
        "ABA fixed-stiffness benchmark",
    )
    require_value(
        aba_fixed,
        "schema",
        "bvk2-aba-fixed-stiffness-validation-v1",
        "ABA fixed-stiffness benchmark",
    )
    require_value(aba_fixed, "model", "bvk2_fixed", "ABA fixed-stiffness benchmark")
    require_value(aba_fixed, "status", "accepted", "ABA fixed-stiffness benchmark")
    require_value(aba_fixed, "adaptive", "false", "ABA fixed-stiffness benchmark")
    require_value(aba_fixed, "c2", "0.16", "ABA fixed-stiffness benchmark")
    require_value(aba_fixed, "nx", "128", "ABA fixed-stiffness benchmark")
    require_pair_cardinality(aba_fixed, ("case_id",), 1, "ABA fixed-stiffness benchmark")

    require_columns(
        cell_stress,
        frozenset(
            {
                "schema",
                "case_id",
                "f",
                "chiN",
                "model",
                "nx",
                "source_period_ratio",
                "signed_period_error",
                "internal_minimum_pass",
                "stress_orientation_pass",
                "root_stress_pass",
                "force_gate_pass",
                "morphology_gate_pass",
                "provenance_pass",
                "accepted",
                "source_artifact",
                "source_fingerprint",
            }
        ),
        "cell-stress mechanism",
    )
    require_value(
        cell_stress,
        "schema",
        "lamellar-cell-stress-mechanism-v2",
        "cell-stress mechanism",
    )
    require_values(
        cell_stress,
        "model",
        frozenset({"uneyama_doi", "liu2019_opf", "burp_ti", "bvk2_fixed", "bvk2"}),
        "cell-stress mechanism",
    )
    require_value(cell_stress, "accepted", "true", "cell-stress mechanism")
    require_pair_cardinality(cell_stress, ("case_id",), 5, "cell-stress mechanism")
    require_pair_cardinality(cell_stress, ("model",), 2, "cell-stress mechanism")
    expected_cell_nx = {
        "uneyama_doi": "64",
        "liu2019_opf": "128",
        "burp_ti": "64",
        "bvk2_fixed": "256",
        "bvk2": "256",
    }
    if any(row["nx"] != expected_cell_nx[row["model"]] for row in cell_stress):
        raise ValueError("cell-stress mechanism model-resolution contract failed")
    expected_bvk_tokens = {
        "bvk2_fixed": ("adaptive=false", "nx=256", "c2=0.16"),
        "bvk2": ("adaptive=true", "nx=256", "c2=0.16"),
    }
    if any(
        any(token not in row["source_fingerprint"] for token in expected_bvk_tokens[row["model"]])
        for row in cell_stress
        if row["model"] in expected_bvk_tokens
    ):
        raise ValueError("cell-stress mechanism adaptive/c2 provenance contract failed")

    require_columns(
        strong,
        frozenset(
            {
                "schema", "claim_kind", "case_id", "f", "chiN", "model",
                "model_label", "calibration_role", "abs_period_error", "profile_rms",
                "period_optimizer", "period_local_minimum_check_pass",
                "period_boundary_limited", "field_gate_pass", "cell_gate_pass",
                "resolution_gate_pass", "composition_gate_pass", "morphology_gate_pass",
                "reduced_field_protocol", "grid_count", "c2", "oversample_factor",
                "sensor_filter_ratio", "discretization_schema", "status",
            }
        ),
        "strong-segregation map",
    )
    require_value(strong, "schema", "liu2019-stress-free-period-map-v2", "strong-segregation map")
    require_value(strong, "claim_kind", "stress_free_lamellar_observable", "strong-segregation map")
    require_values(strong, "model", frozenset(STRONG_MODELS), "strong-segregation map")
    require_value(strong, "status", "accepted", "strong-segregation map")
    require_pair_cardinality(strong, ("f", "chiN"), 6, "strong-segregation map")
    strong_identity = {
        "scft": ("SCFT", "reference", "Polyorder.cell_solve!", None),
        "ohta_kawasaki": ("Ohta-Kawasaki (recomputed)", "liu2019_ok_asymptotic_quadratic_extrapolated_c3_c4", "StressRootBisection", "128"),
        "liu2019_opf": ("OPF (force/stress mapped)", "per_state_force_stress_mapping", "StressRootBisection", "128"),
        "uneyama_doi": ("Uneyama-Doi", "not_applicable", "BootstrapLocal", "128"),
        "burp_ti": ("BURP-TI", "not_applicable", "BoundedKKTNewtonGMRES+AnalyticStressRoot", "128"),
        "bvk2": ("AQCE", "frozen_c2_evaluation", "FilteredUDThetaAnalyticCellStressRoot", "256"),
    }
    for model, (label, role, optimizer, grid) in strong_identity.items():
        rows = [row for row in strong if row["model"] == model]
        if len(rows) != 105:
            raise ValueError(f"strong-segregation model {model} does not contain 105 rows")
        require_state_manifest(rows, STRONG_EXPECTED_STATES, f"strong-segregation model {model}")
        if model == "bvk2":
            require_current_or_legacy_label(
                rows,
                "model_label",
                label,
                frozenset({"APRM", "BVK2"}),
                f"strong-segregation model {model}",
            )
        else:
            require_value(rows, "model_label", label, f"strong-segregation model {model}")
        require_value(rows, "calibration_role", role, f"strong-segregation model {model}")
        require_value(rows, "period_optimizer", optimizer, f"strong-segregation model {model}")
        if grid is not None:
            require_value(rows, "grid_count", grid, f"strong-segregation model {model}")
    strong_bvk2 = [row for row in strong if row["model"] == "bvk2"]
    require_value(strong_bvk2, "c2", "0.16", "strong-segregation BVK2")
    require_value(strong_bvk2, "oversample_factor", "2", "strong-segregation BVK2")
    require_value(strong_bvk2, "sensor_filter_ratio", "12.0", "strong-segregation BVK2")
    require_value(
        strong_bvk2,
        "discretization_schema",
        "ud-theta-spectral-adaptive-k-physical-sensor-filter-oversampled-v2",
        "strong-segregation BVK2",
    )
    require_true(
        strong,
        (
            "period_local_minimum_check_pass", "field_gate_pass", "cell_gate_pass",
            "resolution_gate_pass", "composition_gate_pass", "morphology_gate_pass",
        ),
        "strong-segregation map",
    )
    require_false(strong, ("period_boundary_limited",), "strong-segregation map")
    require_finite(strong, ("f", "chiN", "abs_period_error", "profile_rms"), "strong-segregation map")

    require_columns(
        strong_fixed,
        frozenset(
            {
                "schema", "claim_kind", "case_id", "f", "chiN", "model",
                "model_label", "status", "abs_period_error", "profile_rms",
                "grid_count", "oversample_factor", "sensor_filter_ratio",
                "projected_force_norm", "projected_force_maxabs", "stress_norm",
            }
        ),
        "strong-segregation fixed-stiffness map",
    )
    require_value(strong_fixed, "schema", "bvk-strong-segregation-extension-v1", "strong-segregation fixed-stiffness map")
    require_value(strong_fixed, "claim_kind", "stress_free_lamellar_observable", "strong-segregation fixed-stiffness map")
    require_value(strong_fixed, "model", "bvk2_fixed", "strong-segregation fixed-stiffness map")
    require_current_or_legacy_label(
        strong_fixed,
        "model_label",
        "QCE",
        frozenset({"PRM"}),
        "strong-segregation fixed-stiffness map",
    )
    require_value(strong_fixed, "status", "accepted", "strong-segregation fixed-stiffness map")
    require_value(strong_fixed, "grid_count", "256", "strong-segregation fixed-stiffness map")
    require_value(strong_fixed, "oversample_factor", "2", "strong-segregation fixed-stiffness map")
    require_value(strong_fixed, "sensor_filter_ratio", "12.0", "strong-segregation fixed-stiffness map")
    require_state_manifest(strong_fixed, STRONG_EXPECTED_STATES, "strong-segregation fixed-stiffness map")
    require_finite(
        strong_fixed,
        ("f", "chiN", "abs_period_error", "profile_rms", "projected_force_norm", "stress_norm"),
        "strong-segregation fixed-stiffness map",
    )

    require_columns(
        strong_aggregate,
        frozenset(
            {"schema", "cohort", "state_count", "model", "model_label", "mean_abs_period_error_percent", "mean_profile_rms"}
        ),
        "strong-segregation aggregate",
    )
    require_value(strong_aggregate, "schema", "bvk-strong-segregation-extension-v1", "strong-segregation aggregate")
    require_value(strong_aggregate, "cohort", "integer_chiN_36_50", "strong-segregation aggregate")
    require_value(strong_aggregate, "state_count", "105", "strong-segregation aggregate")
    require_values(strong_aggregate, "model", frozenset({*PRODUCTION_MODELS, "bvk2_fixed"}), "strong-segregation aggregate")
    require_finite(strong_aggregate, ("state_count", "mean_abs_period_error_percent", "mean_profile_rms"), "strong-segregation aggregate")

    expected_validation = {
        "cohort": "integer_chiN_36_50",
        "fixed_row_count": 105,
        "main_row_count": 630,
        "models": ["ohta_kawasaki", "uneyama_doi", "liu2019_opf", "burp_ti", "bvk2_fixed", "bvk2"],
        "schema": "bvk-strong-segregation-extension-v1",
        "state_count": 105,
        "status": "passed",
    }
    if strong_validation != expected_validation:
        raise ValueError("strong-segregation validation-report contract failed")

    require_columns(
        phase_boundaries,
        frozenset({"transition", "chiN", "fA", "bracket_width", "claim_kind", "status", "acceptance_basis", "source"}),
        "phase-boundary ledger",
    )
    require_value(phase_boundaries, "claim_kind", "boundary_bracket", "phase-boundary ledger")
    require_value(phase_boundaries, "status", "accepted", "phase-boundary ledger")
    if len(phase_boundaries) != 149:
        raise ValueError("phase-boundary ledger does not contain 149 accepted roots")
    actual_transition_counts = {
        transition: sum(row["transition"] == transition for row in phase_boundaries)
        for transition in {row["transition"] for row in phase_boundaries}
    }
    if actual_transition_counts != PHASE_TRANSITION_COUNTS:
        raise ValueError("phase-boundary transition manifest failed")
    require_finite(phase_boundaries, ("chiN", "fA", "bracket_width"), "phase-boundary ledger")
    if any(float(row["bracket_width"]) <= 0.0 for row in phase_boundaries):
        raise ValueError("phase-boundary ledger contains a nonpositive bracket width")
    if any(not row["acceptance_basis"].strip() or not row["source"].strip() for row in phase_boundaries):
        raise ValueError("phase-boundary provenance contract failed")

    require_columns(
        scft_phase_reference,
        frozenset({"theory", "transition", "side", "curve_id", "point_index", "fA", "chiN", "role", "source", "source_svg"}),
        "SCFT phase reference",
    )
    require_value(scft_phase_reference, "theory", "SCFT", "SCFT phase reference")
    require_value(
        scft_phase_reference,
        "role",
        "literature reference overlay; not a BVK2 accepted root",
        "SCFT phase reference",
    )
    require_values(scft_phase_reference, "side", frozenset({"left", "right"}), "SCFT phase reference")
    actual_scft_transition_counts = {
        transition: sum(row["transition"] == transition for row in scft_phase_reference)
        for transition in {row["transition"] for row in scft_phase_reference}
    }
    if actual_scft_transition_counts != SCFT_PHASE_TRANSITION_COUNTS:
        raise ValueError("SCFT phase-reference transition manifest failed")
    require_finite(scft_phase_reference, ("curve_id", "point_index", "fA", "chiN"), "SCFT phase reference")
    if any(not row["source"].strip() or not row["source_svg"].strip() for row in scft_phase_reference):
        raise ValueError("SCFT phase-reference provenance contract failed")

    require_unique(production, ("model", "f", "chiN"), "production map")
    require_unique(fixed, ("model", "f", "chiN"), "fixed-stiffness map")
    require_unique(diagnostic, ("f", "chiN", "c2"), "coefficient diagnostic")
    require_unique(diagnostic_ranking, ("c2",), "coefficient ranking")
    require_unique(aba_adaptive_all, ("model", "fA", "chiN"), "ABA benchmark")
    require_unique(aba_fixed, ("model", "fA", "chiN"), "ABA fixed control")
    require_unique(cell_stress, ("case_id", "model"), "cell-stress mechanism")
    require_unique(strong, ("model", "f", "chiN"), "strong-segregation map")
    require_unique(strong_fixed, ("model", "f", "chiN"), "strong-segregation fixed-stiffness map")
    require_unique(strong_aggregate, ("model",), "strong-segregation aggregate")
    require_unique(phase_boundaries, ("transition", "chiN"), "phase-boundary ledger")
    require_unique(scft_phase_reference, ("side", "curve_id", "point_index"), "SCFT phase reference")

    accepted_by_model = {
        model: [
            row
            for row in production
            if row["model"] == model and row["status"] == "accepted"
        ]
        for model in PRODUCTION_MODELS
    }
    common_states = set.intersection(
        *({state(row) for row in rows} for rows in accepted_by_model.values())
    )
    if common_states != COMMON_EXPECTED_STATES:
        raise ValueError(
            "common production cohort differs from the locked 133-state manifest"
        )

    production_common = {
        model: [row for row in rows if state(row) in common_states]
        for model, rows in accepted_by_model.items()
    }
    for model, rows in production_common.items():
        if len(rows) != 133:
            raise ValueError(f"{model} common cohort has {len(rows)} rows")
        require_true(
            rows,
            (
                "field_gate_pass",
                "cell_gate_pass",
                "resolution_gate_pass",
                "composition_gate_pass",
                "morphology_gate_pass",
                "period_local_minimum_check_pass",
            ),
            f"{model} common production cohort",
        )
        require_false(rows, ("period_boundary_limited",), f"{model} common cohort")

    fixed_accepted = [row for row in fixed if row["accepted"] == "true"]
    if len(fixed_accepted) != 146 or {state(row) for row in fixed_accepted} != {
        state(row) for row in accepted_by_model["bvk2"]
    }:
        raise ValueError("fixed-stiffness production map is not the complete 146-state map")
    require_true(
        fixed_accepted,
        (
            "field_gate_pass",
            "cell_gate_pass",
            "composition_gate_pass",
            "morphology_gate_pass",
            "local_minimum_check_pass",
            "stress_orientation_pass",
        ),
        "fixed-stiffness production map",
    )
    fixed_common = [row for row in fixed_accepted if state(row) in common_states]
    if len(fixed_common) != 133:
        raise ValueError("fixed-stiffness common cohort is not 133 states")

    bvk2_full = accepted_by_model["bvk2"]
    if len(bvk2_full) != 146:
        raise ValueError("production BVK2 map is not 146 accepted states")
    require_state_manifest(bvk2_full, FULL_EXPECTED_STATES, "full BVK2 production map")
    require_true(
        bvk2_full,
        (
            "field_gate_pass",
            "cell_gate_pass",
            "resolution_gate_pass",
            "composition_gate_pass",
            "morphology_gate_pass",
            "period_local_minimum_check_pass",
        ),
        "full BVK2 production map",
    )
    require_false(bvk2_full, ("period_boundary_limited",), "full BVK2 production map")
    production_training = {row["case_id"] for row in bvk2_full if row["calibration_role"] == "training"}
    if production_training != TRAINING_STATES_IN_PRODUCTION_MAP:
        raise ValueError("production-map calibration membership has changed")
    bvk2_held_out = [row for row in bvk2_full if row["case_id"] not in production_training]
    if len(bvk2_held_out) != 143:
        raise ValueError("held-out production cohort is not 143 states")

    if len(diagnostic) != 146 * 9:
        raise ValueError("coefficient diagnostic is not a complete 146-state x 9-candidate table")
    require_true(
        diagnostic,
        (
            "stationarity_pass",
            "phase_identity_pass",
            "cell_stationarity_pass",
            "selection_eligible",
            "period_local_minimum_check_pass",
        ),
        "coefficient diagnostic",
    )
    require_false(diagnostic, ("period_boundary_limited",), "coefficient diagnostic")
    if {row["nx"] for row in diagnostic} != {"64"}:
        raise ValueError("coefficient diagnostic is not uniformly labeled nx=64")
    diagnostic_c2 = sorted({float(row["c2"]) for row in diagnostic})
    diagnostic_by_c2 = {
        c2: [row for row in diagnostic if float(row["c2"]) == c2]
        for c2 in diagnostic_c2
    }
    if any(len(rows) != 146 for rows in diagnostic_by_c2.values()):
        raise ValueError("a coefficient-diagnostic candidate lacks the full 146 states")
    for c2, rows in diagnostic_by_c2.items():
        require_state_manifest(
            rows,
            FULL_EXPECTED_STATES,
            f"coefficient diagnostic c2={c2}",
        )
    if {float(row["c2"]) for row in diagnostic_ranking} != set(diagnostic_c2):
        raise ValueError("coefficient ranking and diagnostic candidates differ")
    for row in diagnostic_ranking:
        c2 = float(row["c2"])
        rows = diagnostic_by_c2[c2]
        if int(row["case_count"]) != 146 or row["eligible"] != "True" or row["complete"] != "True":
            raise ValueError(f"coefficient-ranking c2={c2} is incomplete or ineligible")
        if not math.isclose(float(row["period_mae"]), mean(rows, "period_abs_error"), abs_tol=1.0e-15):
            raise ValueError(f"coefficient-ranking period metric drift at c2={c2}")
        if not math.isclose(float(row["profile_rms_mean"]), mean(rows, "profile_rms"), abs_tol=1.0e-15):
            raise ValueError(f"coefficient-ranking profile metric drift at c2={c2}")
    diagnostic_common_by_c2 = {
        c2: [row for row in rows if state(row) in common_states]
        for c2, rows in diagnostic_by_c2.items()
    }
    if any(len(rows) != 133 for rows in diagnostic_common_by_c2.values()):
        raise ValueError("a coefficient-diagnostic candidate lacks the common 133 states")

    adaptive_by_state = {state(row): row for row in bvk2_full}
    fixed_by_state = {state(row): row for row in fixed_accepted}
    high_chin_states = frozenset(
        state_key
        for state_key in set(adaptive_by_state) & set(fixed_by_state)
        if 30.0 <= state_key[1] <= 35.0
    )
    if len(high_chin_states) != 42:
        raise ValueError("high-chiN matched cohort is not 42 states")
    adaptive_high_wins = [
        state_key
        for state_key in high_chin_states
        if float(adaptive_by_state[state_key]["abs_period_error"])
        < float(fixed_by_state[state_key]["abs_period_error"])
    ]
    diagnostic_high_fixed = [
        row
        for row in diagnostic_by_c2[0.0]
        if 30.0 <= float(row["chiN"]) <= 35.0
    ]
    diagnostic_high_adaptive = [
        row
        for row in diagnostic_by_c2[0.16]
        if 30.0 <= float(row["chiN"]) <= 35.0
    ]
    diagnostic_high_fixed_by_state = {
        state(row): row for row in diagnostic_high_fixed
    }
    diagnostic_high_adaptive_by_state = {
        state(row): row for row in diagnostic_high_adaptive
    }
    if (
        set(diagnostic_high_fixed_by_state) != high_chin_states
        or set(diagnostic_high_adaptive_by_state) != high_chin_states
    ):
        raise ValueError("diagnostic and production high-chiN memberships differ")
    diagnostic_adaptive_high_wins = sum(
        float(diagnostic_high_adaptive_by_state[state_key]["period_abs_error"])
        < float(diagnostic_high_fixed_by_state[state_key]["period_abs_error"])
        for state_key in high_chin_states
    )

    aba_adaptive = [row for row in aba_adaptive_all if row["model"] == "bvk2"]
    aba_states = frozenset((0.5, chi_n) for chi_n in ABA_EXPECTED_CHIN)
    for model in {row["model"] for row in aba_adaptive_all}:
        require_state_manifest(
            [row for row in aba_adaptive_all if row["model"] == model],
            aba_states,
            f"ABA adaptive benchmark model={model}",
            f_field="fA",
        )
    require_state_manifest(
        aba_fixed,
        aba_states,
        "ABA fixed-stiffness benchmark",
        f_field="fA",
    )
    require_value(aba_adaptive, "c2", "0.16", "adaptive ABA BVK2 rows")
    require_value(aba_adaptive, "grid", "128", "adaptive ABA BVK2 rows")
    require_true(aba_adaptive, ("field_converged", "accepted"), "adaptive ABA")
    require_true(
        aba_fixed,
        (
            "field_gate_pass",
            "cell_gate_pass",
            "composition_gate_pass",
            "morphology_gate_pass",
            "local_minimum_check_pass",
            "accepted",
        ),
        "fixed-stiffness ABA",
    )
    aba_adaptive_by_chin = {float(row["chiN"]): row for row in aba_adaptive}
    aba_fixed_by_chin = {float(row["chiN"]): row for row in aba_fixed}
    aba_period_wins = sum(
        abs(float(aba_adaptive_by_chin[chi_n]["period_relative_error"]))
        < abs(float(aba_fixed_by_chin[chi_n]["period_relative_error"]))
        for chi_n in ABA_EXPECTED_CHIN
    )

    stress_selected = [
        row for row in cell_stress if row["model"] in {"bvk2_fixed", "bvk2"}
    ]
    if len(stress_selected) != 4 or {row["nx"] for row in stress_selected} != {"256"}:
        raise ValueError("cell-stress ablation is not two states x two nx=256 controls")
    expected_stress_states = frozenset({(0.50, 20.0), (0.35, 30.0)})
    for model in {row["model"] for row in cell_stress}:
        require_state_manifest(
            [row for row in cell_stress if row["model"] == model],
            expected_stress_states,
            f"cell-stress mechanism model={model}",
        )
    require_true(
        stress_selected,
        (
            "internal_minimum_pass",
            "stress_orientation_pass",
            "root_stress_pass",
            "force_gate_pass",
            "morphology_gate_pass",
            "provenance_pass",
            "accepted",
        ),
        "cell-stress mechanism",
    )

    strong_rows_by_model = {
        model: [row for row in strong if row["model"] == model]
        for model in PRODUCTION_MODELS
    }
    strong_rows_by_model["bvk2_fixed"] = strong_fixed
    strong_aggregate_by_model = {row["model"]: row for row in strong_aggregate}
    for model, rows in strong_rows_by_model.items():
        aggregate_row = strong_aggregate_by_model[model]
        if not math.isclose(
            float(aggregate_row["mean_abs_period_error_percent"]),
            100.0 * mean(rows, "abs_period_error"),
            abs_tol=1.0e-14,
        ):
            raise ValueError(f"strong-segregation aggregate period metric drift for {model}")
        if not math.isclose(
            float(aggregate_row["mean_profile_rms"]),
            mean(rows, "profile_rms"),
            abs_tol=1.0e-14,
        ):
            raise ValueError(f"strong-segregation aggregate profile metric drift for {model}")

    phase_name_map = {
        "FCC": "FCC", "BCC": "BCC", "S": "BCC", "C": "CYL",
        "G": "GYR", "O70": "O70", "L": "LAM", "DIS": "DIS",
    }
    phase_graph = {
        phase_name_map[token]
        for transition in PHASE_TRANSITION_COUNTS
        for token in transition.split("/")
    }
    ordered_phases = phase_graph - {"DIS"}
    if ordered_phases != ORDERED_PHASES:
        raise ValueError("phase-boundary diagram topology inference failed")
    scft_phase_name_map = {
        "S_cp": "FCC", "S": "BCC", "C": "CYL", "G": "GYR",
        "O70": "O70", "L": "LAM", "DIS": "DIS",
    }
    scft_phase_graph = {
        scft_phase_name_map[token]
        for transition in SCFT_PHASE_TRANSITION_COUNTS
        for token in transition.split("/")
    }
    if scft_phase_graph - {"DIS"} != ORDERED_PHASES:
        raise ValueError("SCFT phase diagram topology inference failed")
    o70_transitions = {transition for transition in PHASE_TRANSITION_COUNTS if "O70" in transition}
    if o70_transitions != {"C/O70", "G/O70", "O70/L"}:
        raise ValueError("O70 boundary-family manifest failed")
    o70_root_count = sum(PHASE_TRANSITION_COUNTS[transition] for transition in o70_transitions)
    fcc_transitions = {transition for transition in PHASE_TRANSITION_COUNTS if "FCC" in transition}
    if fcc_transitions != {"FCC/BCC", "FCC/DIS"}:
        raise ValueError("FCC boundary-family manifest failed")
    fcc_root_count = sum(PHASE_TRANSITION_COUNTS[transition] for transition in fcc_transitions)

    source_paths = [
        PRODUCTION_MAP,
        FIXED_MAP,
        CELL_STRESS,
        COEFFICIENT_DIAGNOSTIC,
        COEFFICIENT_RANKING,
        ABA_ADAPTIVE,
        ABA_FIXED,
        STRONG_MAP,
        STRONG_FIXED,
        STRONG_AGGREGATE,
        STRONG_VALIDATION,
        PHASE_BOUNDARIES,
        SCFT_PHASE_REFERENCE,
    ]
    source_classes = {
        PRODUCTION_MAP: ("production", "adaptive and comparator lamellar map"),
        FIXED_MAP: ("production_control", "fixed-stiffness lamellar map"),
        CELL_STRESS: ("production_mechanism", "matched cell-stress ablation"),
        COEFFICIENT_DIAGNOSTIC: ("oracle_diagnostic", "historical nx=64 coefficient sweep"),
        COEFFICIENT_RANKING: ("oracle_diagnostic", "stratified coefficient ranking"),
        ABA_ADAPTIVE: ("production_transfer", "adaptive ABA benchmark"),
        ABA_FIXED: ("production_control", "fixed-stiffness ABA benchmark"),
        STRONG_MAP: ("production_strong_segregation", "accepted 105-state strong-segregation map"),
        STRONG_FIXED: ("production_control", "matched fixed-stiffness strong-segregation map"),
        STRONG_AGGREGATE: ("production_strong_segregation", "independently assembled strong-segregation aggregates"),
        STRONG_VALIDATION: ("production_strong_segregation", "strong-segregation validation contract"),
        PHASE_BOUNDARIES: ("production_phase_boundary", "accepted AQCE phase-boundary roots"),
        SCFT_PHASE_REFERENCE: ("literature_reference", "vector SCFT phase-boundary reference"),
    }
    sources = [source_record(path, *source_classes[path]) for path in source_paths]

    claims: list[dict[str, object]] = []
    for model, rows in production_common.items():
        claims.extend(
            (
                claim(
                    f"common133.{model}.period",
                    "production",
                    "common_133",
                    "mean_absolute_period_error",
                    100.0 * mean(rows, "abs_period_error"),
                    "percent",
                    133,
                    [PRODUCTION_MAP.as_posix()],
                    model=model,
                ),
                claim(
                    f"common133.{model}.profile",
                    "production",
                    "common_133",
                    "mean_aligned_profile_rms",
                    mean(rows, "profile_rms"),
                    "fraction",
                    133,
                    [PRODUCTION_MAP.as_posix()],
                    model=model,
                ),
            )
        )
    for model, rows in (("bvk2_fixed", fixed_common),):
        claims.extend(
            (
                claim(
                    f"common133.{model}.period",
                    "production_control",
                    "common_133",
                    "mean_absolute_period_error",
                    100.0 * mean(rows, "abs_period_error"),
                    "percent",
                    133,
                    [FIXED_MAP.as_posix()],
                    model=model,
                ),
                claim(
                    f"common133.{model}.profile",
                    "production_control",
                    "common_133",
                    "mean_aligned_profile_rms",
                    mean(rows, "profile_rms"),
                    "fraction",
                    133,
                    [FIXED_MAP.as_posix()],
                    model=model,
                ),
            )
        )

    for cohort_id, rows, evidence_class in (
        ("held_out_143", bvk2_held_out, "production_held_out"),
        ("full_146_bvk2", bvk2_full, "production"),
        ("full_146_fixed", fixed_accepted, "production_control"),
    ):
        period_field = "abs_period_error"
        model = "bvk2_fixed" if rows is fixed_accepted else "bvk2"
        source_path = FIXED_MAP if rows is fixed_accepted else PRODUCTION_MAP
        claims.extend(
            (
                claim(
                    f"{cohort_id}.period",
                    evidence_class,
                    cohort_id,
                    "mean_absolute_period_error",
                    100.0 * mean(rows, period_field),
                    "percent",
                    len(rows),
                    [source_path.as_posix()],
                    model=model,
                ),
                claim(
                    f"{cohort_id}.profile",
                    evidence_class,
                    cohort_id,
                    "mean_aligned_profile_rms",
                    mean(rows, "profile_rms"),
                    "fraction",
                    len(rows),
                    [source_path.as_posix()],
                    model=model,
                ),
            )
        )

    for c2, rows in diagnostic_by_c2.items():
        label = format(c2, ".2f")
        claims.extend(
            (
                claim(
                    f"coefficient_diagnostic.full146.c2_{label}.period",
                    "oracle_diagnostic",
                    "coefficient_diagnostic_full_146",
                    "mean_absolute_period_error",
                    100.0 * mean(rows, "period_abs_error"),
                    "percent",
                    146,
                    [COEFFICIENT_DIAGNOSTIC.as_posix()],
                    model=f"bvk2_c2_{label}",
                    qualifier="nx=64 historical oracle diagnostic; not a production ranking",
                ),
                claim(
                    f"coefficient_diagnostic.full146.c2_{label}.profile",
                    "oracle_diagnostic",
                    "coefficient_diagnostic_full_146",
                    "mean_aligned_profile_rms",
                    mean(rows, "profile_rms"),
                    "fraction",
                    146,
                    [COEFFICIENT_DIAGNOSTIC.as_posix()],
                    model=f"bvk2_c2_{label}",
                    qualifier="nx=64 historical oracle diagnostic; not a production ranking",
                ),
            )
        )
    for c2 in (0.0, 0.04, 0.16):
        rows = diagnostic_common_by_c2[c2]
        label = format(c2, ".2f")
        claims.extend(
            (
                claim(
                    f"coefficient_diagnostic.common133.c2_{label}.period",
                    "oracle_diagnostic",
                    "common_133",
                    "mean_absolute_period_error",
                    100.0 * mean(rows, "period_abs_error"),
                    "percent",
                    133,
                    [COEFFICIENT_DIAGNOSTIC.as_posix()],
                    model=f"bvk2_c2_{label}",
                    qualifier="nx=64 historical oracle diagnostic; not a production ranking",
                ),
                claim(
                    f"coefficient_diagnostic.common133.c2_{label}.profile",
                    "oracle_diagnostic",
                    "common_133",
                    "mean_aligned_profile_rms",
                    mean(rows, "profile_rms"),
                    "fraction",
                    133,
                    [COEFFICIENT_DIAGNOSTIC.as_posix()],
                    model=f"bvk2_c2_{label}",
                    qualifier="nx=64 historical oracle diagnostic; not a production ranking",
                ),
            )
        )

    high_fixed = [fixed_by_state[key] for key in sorted(high_chin_states)]
    high_adaptive = [adaptive_by_state[key] for key in sorted(high_chin_states)]
    claims.extend(
        (
            claim(
                "high_chiN.state_count",
                "production_control",
                "high_chiN_30_35",
                "matched_state_count",
                42,
                "count",
                42,
                [PRODUCTION_MAP.as_posix(), FIXED_MAP.as_posix()],
            ),
            claim(
                "high_chiN.adaptive_period_wins",
                "production_control",
                "high_chiN_30_35",
                "adaptive_lower_absolute_period_error_count",
                len(adaptive_high_wins),
                "count",
                42,
                [PRODUCTION_MAP.as_posix(), FIXED_MAP.as_posix()],
            ),
            claim(
                "high_chiN.fixed.period",
                "production_control",
                "high_chiN_30_35",
                "mean_absolute_period_error",
                100.0 * mean(high_fixed, "abs_period_error"),
                "percent",
                42,
                [FIXED_MAP.as_posix()],
                model="bvk2_fixed",
            ),
            claim(
                "high_chiN.adaptive.period",
                "production",
                "high_chiN_30_35",
                "mean_absolute_period_error",
                100.0 * mean(high_adaptive, "abs_period_error"),
                "percent",
                42,
                [PRODUCTION_MAP.as_posix()],
                model="bvk2",
            ),
        )
    )
    for f_a in sorted({key[0] for key in high_chin_states}):
        composition_states = [key for key in high_chin_states if key[0] == f_a]
        wins = sum(key in adaptive_high_wins for key in composition_states)
        claims.append(
            claim(
                f"high_chiN.fA_{f_a:.2f}.adaptive_period_wins",
                "production_control",
                "high_chiN_30_35",
                "adaptive_lower_absolute_period_error_count",
                wins,
                "count",
                len(composition_states),
                [PRODUCTION_MAP.as_posix(), FIXED_MAP.as_posix()],
                qualifier=f"f_A={f_a:.2f}",
            )
        )

    claims.extend(
        (
            claim(
                "high_chiN_diagnostic.adaptive_period_wins",
                "oracle_diagnostic",
                "high_chiN_30_35_diagnostic",
                "adaptive_lower_absolute_period_error_count",
                diagnostic_adaptive_high_wins,
                "count",
                42,
                [COEFFICIENT_DIAGNOSTIC.as_posix()],
                qualifier="nx=64 c2=0 versus c2=0.16 oracle diagnostic",
            ),
            claim(
                "high_chiN_diagnostic.fixed.period",
                "oracle_diagnostic",
                "high_chiN_30_35_diagnostic",
                "mean_absolute_period_error",
                100.0 * mean(diagnostic_high_fixed, "period_abs_error"),
                "percent",
                42,
                [COEFFICIENT_DIAGNOSTIC.as_posix()],
                model="bvk2_c2_0.00",
                qualifier="nx=64 oracle diagnostic",
            ),
            claim(
                "high_chiN_diagnostic.adaptive.period",
                "oracle_diagnostic",
                "high_chiN_30_35_diagnostic",
                "mean_absolute_period_error",
                100.0 * mean(diagnostic_high_adaptive, "period_abs_error"),
                "percent",
                42,
                [COEFFICIENT_DIAGNOSTIC.as_posix()],
                model="bvk2_c2_0.16",
                qualifier="nx=64 oracle diagnostic",
            ),
        )
    )

    for model, rows, evidence_class, source_path in (
        ("bvk2_fixed", aba_fixed, "production_control", ABA_FIXED),
        ("bvk2", aba_adaptive, "production_transfer", ABA_ADAPTIVE),
    ):
        claims.extend(
            (
                claim(
                    f"aba7.{model}.period",
                    evidence_class,
                    "aba_7",
                    "mean_absolute_period_error",
                    100.0 * mean(rows, "period_relative_error", absolute=True),
                    "percent",
                    7,
                    [source_path.as_posix()],
                    model=model,
                ),
                claim(
                    f"aba7.{model}.profile",
                    evidence_class,
                    "aba_7",
                    "mean_aligned_profile_rms",
                    mean(rows, "profile_rms"),
                    "fraction",
                    7,
                    [source_path.as_posix()],
                    model=model,
                ),
            )
        )
    claims.append(
        claim(
            "aba7.adaptive_period_wins",
            "production_control",
            "aba_7",
            "adaptive_lower_absolute_period_error_count",
            aba_period_wins,
            "count",
            7,
            [ABA_ADAPTIVE.as_posix(), ABA_FIXED.as_posix()],
        )
    )

    for row in sorted(stress_selected, key=lambda item: (item["case_id"], item["model"])):
        model = row["model"]
        case_id = row["case_id"]
        claims.extend(
            (
                claim(
                    f"cell_stress.{case_id}.{model}.period_ratio",
                    "production_mechanism",
                    "cell_stress_2_state",
                    "period_ratio_to_scft",
                    float(row["source_period_ratio"]),
                    "ratio",
                    1,
                    [CELL_STRESS.as_posix()],
                    model=model,
                    qualifier=case_id,
                ),
                claim(
                    f"cell_stress.{case_id}.{model}.period_error",
                    "production_mechanism",
                    "cell_stress_2_state",
                    "absolute_period_error",
                    100.0 * abs(float(row["signed_period_error"])),
                    "percent",
                    1,
                    [CELL_STRESS.as_posix()],
                    model=model,
                    qualifier=case_id,
                ),
            )
        )

    claims.append(
        claim(
            "strong105.state_count",
            "production_strong_segregation",
            "strong_segregation_105",
            "matched_state_count",
            105,
            "count",
            105,
            [STRONG_MAP.as_posix(), STRONG_AGGREGATE.as_posix(), STRONG_VALIDATION.as_posix()],
            qualifier="seven compositions at every integer chiN from 36 through 50",
        )
    )
    for model, rows in strong_rows_by_model.items():
        evidence_class = "production_control" if model == "bvk2_fixed" else "production_strong_segregation"
        claim_sources = (
            [STRONG_FIXED.as_posix()]
            if model == "bvk2_fixed"
            else [STRONG_MAP.as_posix(), STRONG_AGGREGATE.as_posix(), STRONG_VALIDATION.as_posix()]
        )
        claims.extend(
            (
                claim(
                    f"strong105.{model}.period",
                    evidence_class,
                    "strong_segregation_105",
                    "mean_absolute_period_error",
                    100.0 * mean(rows, "abs_period_error"),
                    "percent",
                    105,
                    claim_sources,
                    model=model,
                ),
                claim(
                    f"strong105.{model}.profile",
                    evidence_class,
                    "strong_segregation_105",
                    "mean_aligned_profile_rms",
                    mean(rows, "profile_rms"),
                    "fraction",
                    105,
                    claim_sources,
                    model=model,
                ),
            )
        )

    claims.extend(
        (
            claim(
                "phase149.accepted_root_count",
                "production_phase_boundary",
                "phase_boundary_149",
                "accepted_boundary_root_count",
                149,
                "count",
                149,
                [PHASE_BOUNDARIES.as_posix()],
                model="bvk2",
                qualifier="accepted boundary roots; no adaptive-stiffness phase ablation",
            ),
            claim(
                "phase149.o70_boundary_root_count",
                "production_phase_boundary",
                "phase_boundary_149",
                "o70_boundary_root_count",
                o70_root_count,
                "count",
                149,
                [PHASE_BOUNDARIES.as_posix()],
                model="bvk2",
                qualifier="roots on C/O70, G/O70, and O70/L boundary families",
            ),
            claim(
                "phase149.o70_boundary_family_count",
                "topology_inference",
                "phase_boundary_149",
                "o70_boundary_family_count",
                len(o70_transitions),
                "count",
                149,
                [PHASE_BOUNDARIES.as_posix(), SCFT_PHASE_REFERENCE.as_posix()],
                model="bvk2",
                qualifier="inferred from the accepted transition manifest; not adaptive-stiffness causality",
            ),
            claim(
                "phase149.ordered_phase_count",
                "topology_inference",
                "phase_boundary_149",
                "ordered_phase_count",
                len(ordered_phases),
                "count",
                149,
                [PHASE_BOUNDARIES.as_posix(), SCFT_PHASE_REFERENCE.as_posix()],
                model="bvk2",
                qualifier="common six-ordered-phase set plus DIS inferred from the AQCE and SCFT boundary graphs; not adaptive-stiffness causality",
            ),
            claim(
                "phase149.fcc_boundary_root_count",
                "production_phase_boundary",
                "phase_boundary_149",
                "fcc_boundary_root_count",
                fcc_root_count,
                "count",
                149,
                [PHASE_BOUNDARIES.as_posix()],
                model="bvk2",
                qualifier="roots on FCC/BCC and FCC/DIS boundary families",
            ),
            claim(
                "phase149.fcc_boundary_family_count",
                "topology_inference",
                "phase_boundary_149",
                "fcc_boundary_family_count",
                len(fcc_transitions),
                "count",
                149,
                [PHASE_BOUNDARIES.as_posix(), SCFT_PHASE_REFERENCE.as_posix()],
                model="bvk2",
                qualifier="FCC pocket support inferred from accepted AQCE roots and the SCFT reference graph; not adaptive-stiffness causality",
            ),
        )
    )

    ranked = sorted(diagnostic_ranking, key=lambda row: int(row["rank"]))
    if [int(row["rank"]) for row in ranked] != list(range(1, 10)):
        raise ValueError("coefficient diagnostic ranks are not exactly 1 through 9")
    oracle = ranked[0]
    baseline = next(row for row in ranked if float(row["c2"]) == 0.16)

    report = {
        "schema": "bvk-contribution-claim-evidence-v1",
        "generator": "scripts/generate_bvk_contribution_claim_evidence.py",
        "sources": sources,
        "evidence_policy": {
            "production_label": "Publication metrics use accepted production-resolution ledgers.",
            "diagnostic_label": "The nx=64 coefficient sweep is an oracle diagnostic and is not a production-model ranking.",
            "topology_label": "Accepted phase-boundary roots are direct production evidence; the ordered-phase set and connectivity are inferences from their transition graph and are not attributed specifically to adaptive stiffness.",
            "simulation_decision": "No additional simulation is required for the contribution reframing; current evidence passes the factual-record gate.",
        },
        "cohorts": {
            "common_133": {
                "state_count": 133,
                "membership": state_list(common_states),
                "definition": "Exact accepted-state intersection of the five production comparator models.",
                "permitted_evidence_classes": [
                    "oracle_diagnostic",
                    "production",
                    "production_control",
                ],
            },
            "held_out_143": {
                "state_count": 143,
                "excluded_calibration_case_ids": sorted(production_training),
                "membership_case_ids": sorted(row["case_id"] for row in bvk2_held_out),
                "definition": "Complete accepted AQCE map excluding the three calibration states represented in that map.",
                "permitted_evidence_classes": ["production_held_out"],
            },
            "full_146_bvk2": {
                "state_count": 146,
                "membership": state_list(state(row) for row in bvk2_full),
                "permitted_evidence_classes": ["production"],
            },
            "full_146_fixed": {
                "state_count": 146,
                "membership": state_list(state(row) for row in fixed_accepted),
                "permitted_evidence_classes": ["production_control"],
            },
            "coefficient_diagnostic_full_146": {
                "state_count_per_candidate": 146,
                "candidate_c2": diagnostic_c2,
                "resolution_nx": 64,
                "evidence_class": "oracle_diagnostic",
                "permitted_evidence_classes": ["oracle_diagnostic"],
            },
            "high_chiN_30_35": {
                "state_count": 42,
                "membership": state_list(high_chin_states),
                "evidence_class": "production_control",
                "permitted_evidence_classes": ["production", "production_control"],
            },
            "high_chiN_30_35_diagnostic": {
                "state_count": 42,
                "membership": state_list(high_chin_states),
                "evidence_class": "oracle_diagnostic",
                "resolution_nx": 64,
                "permitted_evidence_classes": ["oracle_diagnostic"],
            },
            "aba_7": {
                "state_count": 7,
                "f_A": 0.5,
                "chiN": list(ABA_EXPECTED_CHIN),
                "permitted_evidence_classes": [
                    "production_control",
                    "production_transfer",
                ],
            },
            "cell_stress_2_state": {
                "state_count": 2,
                "membership": state_list(
                    {(float(row["f"]), float(row["chiN"])) for row in stress_selected}
                ),
                "resolution_nx": 256,
                "permitted_evidence_classes": ["production_mechanism"],
            },
            "strong_segregation_105": {
                "state_count": 105,
                "membership": state_list(STRONG_EXPECTED_STATES),
                "chiN": list(range(36, 51)),
                "permitted_evidence_classes": [
                    "production_control",
                    "production_strong_segregation",
                ],
            },
            "phase_boundary_149": {
                "state_count": 149,
                "transition_counts": PHASE_TRANSITION_COUNTS,
                "ordered_phases": sorted(ORDERED_PHASES),
                "permitted_evidence_classes": [
                    "production_phase_boundary",
                    "topology_inference",
                ],
            },
        },
        "gate_audit": {
            "common_133_production": {
                "pass": True,
                "row_count": 133 * len(PRODUCTION_MODELS),
                "gates": [
                    "field",
                    "cell_stress",
                    "resolution",
                    "composition",
                    "morphology",
                    "local_period_minimum",
                    "not_boundary_limited",
                ],
            },
            "full_146_bvk2_production": {
                "pass": True,
                "row_count": 146,
                "gates": [
                    "field",
                    "cell_stress",
                    "resolution",
                    "composition",
                    "morphology",
                    "local_period_minimum",
                    "not_boundary_limited",
                ],
            },
            "held_out_143_bvk2_production": {
                "pass": True,
                "row_count": 143,
                "gates": [
                    "field",
                    "cell_stress",
                    "resolution",
                    "composition",
                    "morphology",
                    "local_period_minimum",
                    "not_boundary_limited",
                ],
            },
            "fixed_146_production": {
                "pass": True,
                "row_count": 146,
                "gates": [
                    "field",
                    "cell_stress",
                    "composition",
                    "morphology",
                    "local_period_minimum",
                    "stress_orientation",
                ],
            },
            "coefficient_diagnostic": {
                "pass": True,
                "row_count": 1314,
                "gates": [
                    "stationarity",
                    "phase_identity",
                    "cell_stationarity",
                    "selection_eligible",
                    "local_period_minimum",
                    "not_boundary_limited",
                ],
            },
            "aba_fixed_adaptive": {
                "pass": True,
                "row_count": 14,
                "note": "The adaptive ledger exposes field-converged/accepted gates; the fixed control additionally exposes composition, morphology, cell, and local-minimum gates.",
            },
            "cell_stress_2_state": {
                "pass": True,
                "row_count": 4,
                "gates": [
                    "internal_minimum",
                    "stress_orientation",
                    "root_stress",
                    "force",
                    "morphology",
                    "provenance",
                ],
            },
            "strong_segregation_105": {
                "pass": True,
                "row_count": 735,
                "gates": [
                    "exact_state_manifest",
                    "accepted_status",
                    "field",
                    "cell_stress",
                    "resolution",
                    "composition",
                    "morphology",
                    "local_period_minimum",
                    "not_boundary_limited",
                    "aggregate_reconstruction",
                    "validation_report",
                ],
            },
            "phase_boundary_149": {
                "pass": True,
                "row_count": 149,
                "gates": [
                    "accepted_status",
                    "boundary_bracket_claim_kind",
                    "positive_bracket_width",
                    "finite_coordinates",
                    "exact_transition_manifest",
                    "nonempty_provenance",
                ],
            },
            "scft_phase_reference": {
                "pass": True,
                "row_count": 1334,
                "gates": [
                    "SCFT_identity",
                    "exact_transition_manifest",
                    "finite_coordinates",
                    "nonempty_provenance",
                    "literature_reference_not_production_root",
                ],
            },
        },
        "diagnostic_selection": {
            "full_146_oracle_c2": float(oracle["c2"]),
            "selection_metric": "stratified_balanced_loss",
            "oracle_stratified_balanced_loss": float(oracle["stratified_balanced_loss"]),
            "frozen_c2": 0.16,
            "frozen_stratified_balanced_loss": float(baseline["stratified_balanced_loss"]),
            "note": "Candidate summaries are descriptive oracle diagnostics; production remains frozen at c2=0.16.",
        },
        "claims": sorted(claims, key=lambda item: str(item["claim_id"])),
    }
    require_claim_integrity(report["claims"], report["cohorts"], sources)
    return report


def serialized_evidence() -> str:
    return json.dumps(build_evidence(), indent=2, sort_keys=True, allow_nan=False) + "\n"


def output_label(path: Path) -> str:
    try:
        return path.relative_to(PROJECT).as_posix()
    except ValueError:
        return str(path)


def write_or_check_output(output: Path, *, check: bool) -> None:
    resolved = output if output.is_absolute() else PROJECT / output
    payload = serialized_evidence()
    if check:
        if not resolved.is_file():
            raise SystemExit(f"missing checked evidence artifact: {resolved}")
        if resolved.read_text(encoding="utf-8") != payload:
            raise SystemExit(f"claim-evidence artifact is stale: {resolved}")
        print(f"PASS: {output_label(resolved)} matches canonical evidence")
        return
    resolved.parent.mkdir(parents=True, exist_ok=True)
    resolved.write_text(payload, encoding="utf-8")
    print(f"Wrote {output_label(resolved)}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--check",
        action="store_true",
        help="fail if a fresh reconstruction differs from the checked artifact",
    )
    args = parser.parse_args()
    write_or_check_output(args.output, check=args.check)


if __name__ == "__main__":
    main()
