#!/usr/bin/env python3
"""Freeze BVK2 publication-boundary resolution and uncertainty metadata.

This is an audit/aggregation script.  It never promotes or demotes canonical
phase-boundary rows and never launches a field solve.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
from collections import defaultdict
from collections.abc import Iterable
from pathlib import Path
from typing import Any

import mpltex

PROJECT_ROOT = Path(__file__).resolve().parents[1]
RESULTS = PROJECT_ROOT / "results"
PUBLICATION_LEDGER = (
    RESULTS / "bvk2_publication_phase_diagram/accepted_phase_boundaries.csv"
)
ACTIVE_PROVENANCE = RESULTS / "bvk2_publication_phase_diagram/boundary_diagnostics.csv"
CORRECTION_ROOT = RESULTS / "bvk2_richardson_correction"
DEFAULT_OUTPUT_DIR = RESULTS / "bvk2_resolution_uncertainty_table"
SCHEMA = "bvk2-resolution-uncertainty-table-v1"

STRESS_TARGET = 1.0e-3
STRESS_MAX = 2.0e-3
FILTERED_GC_SCHEMA = "bvk2-ud-theta-gc-global-sensor-filter-r6-v1"


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def read_one(path: Path) -> dict[str, str]:
    rows = read_rows(path)
    if len(rows) != 1:
        raise ValueError(f"expected one row in {path}, found {len(rows)}")
    return rows[0]


def unique_record(rows: Iterable[dict[str, str]], *, context: str) -> dict[str, str]:
    selected = list(rows)
    if len(selected) != 1:
        raise ValueError(f"expected one {context} record, found {len(selected)}")
    return selected[0]


def fvalue(row: dict[str, str], key: str) -> float:
    value = float(row[key])
    if not math.isfinite(value):
        raise ValueError(f"{key} is nonfinite in {row}")
    return value


def truth(value: object) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def resolve_repository_path(value: str | Path) -> Path:
    """Resolve repository-relative and legacy absolute artifact paths portably."""
    path = Path(value)
    if not path.is_absolute():
        path = PROJECT_ROOT / path
    resolved = path.resolve()
    try:
        resolved.relative_to(PROJECT_ROOT.resolve())
    except ValueError:
        parts = resolved.parts
        if "results" not in parts:
            raise ValueError(f"path is outside the repository: {value}") from None
        resolved = (PROJECT_ROOT / Path(*parts[parts.index("results") :])).resolve()
        try:
            resolved.relative_to(PROJECT_ROOT.resolve())
        except ValueError:
            raise ValueError(f"path is outside the repository: {value}") from None
    return resolved


def relative(path: str | Path) -> str:
    return resolve_repository_path(path).relative_to(PROJECT_ROOT.resolve()).as_posix()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def evidence(paths: Iterable[Path]) -> tuple[str, str]:
    unique: list[Path] = []
    for path in paths:
        resolved = path.resolve()
        if resolved not in unique:
            if not resolved.is_file():
                raise ValueError(f"missing evidence file: {resolved}")
            unique.append(resolved)
    return (
        ";".join(relative(path) for path in unique),
        ";".join(sha256(path) for path in unique),
    )


def keyed(rows: Iterable[dict[str, str]], key: str = "chiN") -> dict[float, dict[str, str]]:
    result: dict[float, dict[str, str]] = {}
    for row in rows:
        numeric = fvalue(row, key)
        if numeric in result:
            raise ValueError(f"duplicate {key}={numeric}")
        result[numeric] = row
    return result


def stress_tier(maximum: float | None) -> str:
    if maximum is None or not math.isfinite(maximum):
        return "not_tabulated"
    if maximum <= STRESS_TARGET:
        return "target"
    if maximum <= STRESS_MAX:
        return "fallback"
    return "above_current_maximum"


def finite_stress(value: object, *, context: str) -> float:
    if value is None or str(value).strip() == "":
        raise ValueError(f"missing cell stress for {context}")
    stress = float(value)
    if not math.isfinite(stress):
        raise ValueError(f"nonfinite cell stress for {context}")
    return abs(stress)


def validate_current_stress(values: Iterable[object], *, context: str) -> float:
    """Return the maximum stress only when a current/unified label is valid."""
    stresses = [finite_stress(value, context=context) for value in values]
    if not stresses:
        raise ValueError(f"missing cell stress for {context}")
    maximum = max(stresses)
    if maximum > STRESS_MAX:
        raise ValueError(
            f"cell stress exceeds current hard maximum for {context}: {maximum}"
        )
    return maximum


def max_stress(
    path: Path,
    chi_n: float,
    coordinates: tuple[float, float],
    column_pairs: tuple[tuple[str, str], ...],
) -> float:
    candidates = [
        row for row in read_rows(path)
        if math.isclose(fvalue(row, "chiN"), chi_n, abs_tol=1.0e-10)
    ]
    rows = exact_coordinate_rows(
        candidates,
        coordinates=coordinates,
        coordinate_key="fA",
        per_coordinate=1,
        context=f"chiN={chi_n} endpoints in {path}",
        reject_extra=True,
    )
    for first, second in column_pairs:
        if first in rows[0] and second in rows[0]:
            return max(
                finite_stress(row.get(key), context=f"{path}:{key}")
                for row in rows
                for key in (first, second)
            )
    raise ValueError(f"no recognized stress columns in {path}")


def exact_coordinate_rows(
    rows: Iterable[dict[str, str]],
    *,
    coordinates: Iterable[float],
    coordinate_key: str,
    per_coordinate: int,
    context: str,
    reject_extra: bool = False,
) -> list[dict[str, str]]:
    candidates = list(rows)
    expected = tuple(float(value) for value in coordinates)
    if len(expected) != len(set(expected)):
        raise ValueError(f"duplicate expected endpoint coordinate for {context}")
    selected: list[dict[str, str]] = []
    for coordinate in expected:
        matches = [
            row
            for row in candidates
            if math.isclose(
                fvalue(row, coordinate_key), coordinate, abs_tol=2.0e-10
            )
        ]
        if len(matches) != per_coordinate:
            raise ValueError(
                f"expected {per_coordinate} unique rows at {coordinate} for "
                f"{context}, found {len(matches)}"
            )
        selected.extend(matches)
    if len(selected) != per_coordinate * len(expected) or (
        reject_extra and len(candidates) != len(selected)
    ):
        raise ValueError(f"endpoint coordinate set is not exact for {context}")
    return selected


def _same_number(
    first: dict[str, str],
    first_key: str,
    second: dict[str, str],
    second_key: str,
    *,
    context: str,
    tolerance: float = 2.0e-10,
) -> None:
    if not math.isclose(
        fvalue(first, first_key),
        fvalue(second, second_key),
        abs_tol=tolerance,
    ):
        raise ValueError(
            f"numeric provenance mismatch for {context}: "
            f"{first_key} != {second_key}"
        )


def validate_filtered_gc_contract(
    boundary_rows: list[dict[str, str]],
    endpoint_rows: list[dict[str, str]],
    manifest_rows: list[dict[str, str]],
) -> dict[float, tuple[dict[str, str], dict[str, str]]]:
    """Validate all filtered-G/C endpoints and reconstruct canonical brackets."""
    boundaries = keyed(boundary_rows)
    accepted = [row for row in endpoint_rows if row.get("status") == "accepted"]
    identities = [(fvalue(row, "chiN"), fvalue(row, "fA")) for row in accepted]
    if len(identities) != len(set(identities)):
        raise ValueError("filtered G/C accepted endpoint identities are not unique")

    manifest_by_job: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in manifest_rows:
        manifest_by_job[row["job_id"]].append(row)
    accepted_jobs = {row["job_id"] for row in accepted}
    if set(manifest_by_job) != accepted_jobs:
        raise ValueError("filtered G/C manifest does not join accepted endpoints exactly")

    validated: dict[float, list[dict[str, str]]] = defaultdict(list)
    for endpoint in accepted:
        matches = manifest_by_job.get(endpoint["job_id"], [])
        if len(matches) != 1:
            raise ValueError(
                f"filtered G/C manifest join is not unique for {endpoint['job_id']}"
            )
        manifest = matches[0]
        if endpoint.get("transition") != "G/C" or endpoint.get("claim_kind") != "phase_endpoint_pair":
            raise ValueError("filtered G/C endpoint claim vocabulary mismatch")
        if endpoint.get("status_reason") != "same_schema_filtered_phase_pair":
            raise ValueError("filtered G/C endpoint acceptance reason mismatch")
        if manifest.get("claim_kind") != "phase_endpoint_pair":
            raise ValueError("filtered G/C manifest claim mismatch")
        _same_number(endpoint, "chiN", manifest, "chiN", context="endpoint/manifest chiN")
        _same_number(endpoint, "fA", manifest, "fA", context="endpoint/manifest fA")
        _same_number(
            endpoint,
            "sensor_filter_ratio",
            manifest,
            "sensor_filter_ratio",
            context="endpoint/manifest filter ratio",
        )
        if not math.isclose(fvalue(endpoint, "sensor_filter_ratio"), 6.0, abs_tol=1e-12):
            raise ValueError("filtered G/C endpoint does not use r=6")
        require_same_source(
            endpoint["artifact"], manifest["artifact"], context="filtered G/C artifact"
        )
        artifact = resolve_repository_path(manifest["artifact"])
        outdir = resolve_repository_path(manifest["outdir"])
        if artifact.parent != outdir:
            raise ValueError("filtered G/C pair artifact is outside its manifest outdir")
        raw = read_one(artifact)
        required_raw = {
            "schema": FILTERED_GC_SCHEMA,
            "transition": "G/C",
            "claim_kind": "phase_endpoint_pair",
            "gyr_grid": "96x96x96",
            "cyl_grid": "128x224",
            "gyr_status": "accepted_target",
            "cyl_status": "accepted_target",
            "status": "accepted",
            "status_reason": "same_schema_filtered_phase_pair",
        }
        for field, expected in required_raw.items():
            if raw.get(field) != expected:
                raise ValueError(
                    f"filtered G/C raw pair {field} mismatch for {endpoint['job_id']}"
                )
        for field in (
            "gyr_plateau_pass",
            "cyl_plateau_pass",
            "gyr_phase_identity_pass",
            "cyl_phase_identity_pass",
        ):
            if not truth(raw.get(field)):
                raise ValueError(
                    f"filtered G/C raw pair {field} failed for {endpoint['job_id']}"
                )
        for endpoint_key, raw_key in (
            ("chiN", "chiN"),
            ("fA", "fA"),
            ("sensor_filter_ratio", "sensor_filter_ratio"),
            ("gyr_energy_density", "gyr_energy_density"),
            ("cyl_energy_density", "cyl_energy_density"),
            ("delta_f", "delta_f"),
            ("gyr_cell_stress", "gyr_cell_stress"),
            ("cyl_cell_stress", "cyl_cell_stress"),
        ):
            _same_number(
                endpoint,
                endpoint_key,
                raw,
                raw_key,
                context=f"endpoint/raw {endpoint_key}",
                tolerance=2.0e-12,
            )
        if not math.isclose(fvalue(raw, "sensor_filter_ratio"), 6.0, abs_tol=1e-12):
            raise ValueError("filtered G/C raw pair does not use r=6")
        if not math.isclose(
            fvalue(raw, "gyr_energy_density") - fvalue(raw, "cyl_energy_density"),
            fvalue(raw, "delta_f"),
            abs_tol=2.0e-12,
        ):
            raise ValueError("filtered G/C raw pair delta is internally inconsistent")
        for endpoint_key, raw_key in (
            ("gyr_cell_root", "gyr_cell_root"),
            ("cyl_cell_root", "cyl_cell_root"),
        ):
            require_same_source(
                endpoint[endpoint_key], raw[raw_key], context=f"filtered G/C {raw_key}"
            )
            root = resolve_repository_path(raw[raw_key])
            if not is_within(root, outdir):
                raise ValueError(f"filtered G/C {raw_key} is outside its manifest outdir")
            if not root.is_file():
                raise ValueError(f"filtered G/C {raw_key} is missing: {root}")
        validated[fvalue(endpoint, "chiN")].append(endpoint)

    if set(validated) != set(boundaries):
        raise ValueError("filtered G/C endpoint slices do not match boundary summaries")
    selected: dict[float, tuple[dict[str, str], dict[str, str]]] = {}
    for chi_n, endpoints in validated.items():
        ordered = sorted(endpoints, key=lambda row: fvalue(row, "fA"))
        candidates: list[tuple[float, dict[str, str], dict[str, str]]] = []
        for lower, upper in zip(ordered, ordered[1:], strict=False):
            delta_lower = fvalue(lower, "delta_f")
            delta_upper = fvalue(upper, "delta_f")
            if delta_lower * delta_upper < 0.0:
                candidates.append(
                    (fvalue(upper, "fA") - fvalue(lower, "fA"), lower, upper)
                )
        if not candidates:
            raise ValueError(f"filtered G/C has no adjacent sign change at chiN={chi_n}")
        minimum = min(width for width, _lower, _upper in candidates)
        narrowest = [
            pair for pair in candidates if math.isclose(pair[0], minimum, abs_tol=2e-12)
        ]
        if len(narrowest) != 1:
            raise ValueError(
                f"filtered G/C minimum-width sign-changing pair is tied at chiN={chi_n}"
            )
        width, lower, upper = narrowest[0]
        delta_lower = fvalue(lower, "delta_f")
        delta_upper = fvalue(upper, "delta_f")
        secant = fvalue(lower, "fA") - delta_lower * width / (
            delta_upper - delta_lower
        )
        boundary = boundaries[chi_n]
        checks = (
            (fvalue(boundary, "coordinate_lower"), fvalue(lower, "fA"), "lower"),
            (fvalue(boundary, "coordinate_upper"), fvalue(upper, "fA"), "upper"),
            (fvalue(boundary, "delta_lower"), delta_lower, "delta_lower"),
            (fvalue(boundary, "delta_upper"), delta_upper, "delta_upper"),
            (2.0 * fvalue(boundary, "coordinate_half_width"), width, "width"),
            (fvalue(boundary, "coordinate"), secant, "secant"),
        )
        for reported, reconstructed, label in checks:
            if not math.isclose(reported, reconstructed, abs_tol=2e-10):
                raise ValueError(
                    f"filtered G/C boundary {label} is stale at chiN={chi_n}"
                )
        if not (
            boundary.get("transition") == "G/C"
            and boundary.get("claim_kind") == "boundary_bracket"
            and boundary.get("status") == "accepted"
            and boundary.get("energy_contract")
            == "same_schema_filtered_gyr96_minus_cyl128x224_r6"
            and math.isclose(
                fvalue(boundary, "sensor_filter_ratio"), 6.0, abs_tol=1e-12
            )
        ):
            raise ValueError(f"filtered G/C boundary contract mismatch at chiN={chi_n}")
        selected[chi_n] = (lower, upper)
    return selected


def infer_coarse_dims(fine_dims: str, dx_coarse: float, dx_fine: float) -> str:
    sizes = [int(part) for part in fine_dims.split("x")]
    ratio = dx_fine / dx_coarse
    coarse = [round(size * ratio) for size in sizes]
    if any(size <= 0 for size in coarse):
        raise ValueError(f"cannot infer coarse dimensions from {fine_dims}")
    return "x".join(str(size) for size in coarse)


def load_corrections() -> list[tuple[dict[str, str], Path]]:
    result: list[tuple[dict[str, str], Path]] = []
    for path in sorted(CORRECTION_ROOT.glob("**/corrected_roots.csv")):
        for row in read_rows(path):
            result.append((row, path))
    return result


def match_correction(
    publication: dict[str, str],
    corrections: list[tuple[dict[str, str], Path]],
    active_provenance: list[dict[str, str]],
) -> tuple[dict[str, str], Path]:
    """Match exactly one correction selected by the active provenance ledger."""
    legacy_source = RESULTS / "bvk2_gc_gl_boundary_campaign"
    if resolve_repository_path(publication["source"]) != legacy_source.resolve():
        raise ValueError(
            "Richardson correction requested for a nonlegacy canonical source: "
            f"{publication}"
        )
    provenance_matches = [
        row
        for row in active_provenance
        if row["transition"] == publication["transition"]
        and math.isclose(
            fvalue(row, "chiN"), fvalue(publication, "chiN"), abs_tol=1.0e-10
        )
        and row.get("fA_corrected", "").strip()
        and math.isclose(
            fvalue(row, "fA_corrected"),
            fvalue(publication, "fA"),
            abs_tol=2.0e-10,
        )
        and row.get("status") == "accepted"
        and row.get("gate_reason") == "accepted"
        and resolve_repository_path(row["source"]) == legacy_source.resolve()
    ]
    if len(provenance_matches) != 1:
        raise ValueError(
            "could not uniquely match active Richardson provenance by transition, "
            f"chiN, coordinate, and source for {publication}; found "
            f"{len(provenance_matches)}"
        )
    matches = [
        (row, path)
        for row, path in corrections
        if row["transition"] == publication["transition"]
        and math.isclose(
            fvalue(row, "chiN"), fvalue(publication, "chiN"), abs_tol=1.0e-10
        )
        and math.isclose(
            fvalue(row, "root_fA_corrected"),
            fvalue(publication, "fA"),
            abs_tol=2.0e-10,
        )
    ]
    if len(matches) != 1:
        raise ValueError(
            "could not uniquely match canonical correction by transition, chiN, "
            f"coordinate, and active provenance for {publication}; found "
            f"{len(matches)}"
        )
    return matches[0]


def load_sdis_points() -> list[tuple[dict[str, str], Path]]:
    directory = RESULTS / "bvk2_sdis_boundary_campaign/points"
    return [(read_one(path), path) for path in sorted(directory.glob("*.csv"))]


def match_sdis(
    publication: dict[str, str],
    points: list[tuple[dict[str, str], Path]],
) -> tuple[dict[str, str], Path]:
    chi_n = fvalue(publication, "chiN")
    f_a = fvalue(publication, "fA")
    matches = []
    for point, path in points:
        direction = point["direction"]
        if direction == "fA":
            same = (
                math.isclose(fvalue(point, "root_chiN"), chi_n, abs_tol=1.0e-10)
                and math.isclose(fvalue(point, "root_fA"), f_a, abs_tol=1.0e-10)
            )
        elif direction == "chiN":
            same = (
                math.isclose(fvalue(point, "root_chiN"), chi_n, abs_tol=1.0e-10)
                and math.isclose(fvalue(point, "target_fA"), f_a, abs_tol=1.0e-10)
            )
        else:
            raise ValueError(f"unknown S/DIS direction {direction}")
        if same:
            matches.append((point, path))
    if len(matches) != 1:
        raise ValueError(f"could not uniquely match S/DIS publication row {publication}")
    return matches[0]


def match_fixed_fa_point(
    publication: dict[str, str], directory: Path, pattern: str
) -> tuple[dict[str, str], Path]:
    matches: list[tuple[dict[str, str], Path]] = []
    for path in sorted(directory.glob(pattern)):
        point = read_one(path)
        if (
            point["transition"] == publication["transition"]
            and math.isclose(
                fvalue(point, "root_fA"),
                fvalue(publication, "fA"),
                abs_tol=2.0e-10,
            )
            and math.isclose(
                fvalue(point, "root_chiN"),
                fvalue(publication, "chiN"),
                abs_tol=2.0e-10,
            )
        ):
            matches.append((point, path))
    if len(matches) != 1:
        raise ValueError(
            f"could not uniquely match fixed-fA point {publication}; "
            f"found {len(matches)}"
        )
    source = resolve_repository_path(publication["source"])
    if source != matches[0][1].resolve() and source != directory.parent.resolve():
        raise ValueError(
            f"fixed-fA point does not match canonical source for {publication}"
        )
    return matches[0]


def is_within(path: Path, directory: Path) -> bool:
    try:
        path.resolve().relative_to(directory.resolve())
    except ValueError:
        return False
    return True


def classify_canonical_source(
    transition: str,
    source_value: str,
    *,
    direct_gl_sources: set[Path],
) -> str:
    """Return an explicit source family or fail closed."""
    source = resolve_repository_path(source_value)
    exact: dict[tuple[str, Path], str] = {
        ("S/DIS", (RESULTS / "bvk2_sdis_boundary_campaign").resolve()): "sdis",
        (
            "S/DIS",
            (RESULTS / "bvk2_sdis_fixed_fa_chi_refinement").resolve(),
        ): "sdis_fixed_fa",
        ("S/C", (RESULTS / "bvk2_gc_gl_boundary_campaign").resolve()): "richardson",
        ("G/C", (RESULTS / "bvk2_gc_gl_boundary_campaign").resolve()): "richardson",
        (
            "G/C",
            (RESULTS / "bvk2_ud_theta_gc_filtered_r6/boundary_points.csv").resolve(),
        ): "gc_filtered",
        (
            "G/C",
            (
                RESULTS
                / "bvk2_ud_theta_g_pocket_highchi_ledger/accepted_gc_boundaries.csv"
            ).resolve(),
        ): "gc_high",
        (
            "G/C",
            (RESULTS / "bvk2_ud_theta_gc_chi54_direct/boundary_point.csv").resolve(),
        ): "gc54",
        (
            "G/L",
            (RESULTS / "bvk2_fcc_o70_boundary_campaign/boundary_brackets.csv").resolve(),
        ): "gl15",
        (
            "G/L",
            (
                RESULTS
                / "bvk2_ud_theta_gl_matched_spacing_audit/matched_boundary_promotions.csv"
            ).resolve(),
        ): "gl_matched",
    }
    if (transition, source) in exact:
        return exact[(transition, source)]
    fixed_sc_root = RESULTS / "bvk2_sc_fixed_fa_chi_continuation"
    if transition == "S/C" and is_within(source, fixed_sc_root):
        return "sc_fixed_fa"
    sc_high_root = RESULTS / "bvk2_ud_theta_sc_highchi"
    if transition == "S/C" and is_within(source, sc_high_root):
        return "sc_high"
    if transition == "G/L" and source in direct_gl_sources:
        return "gl_direct"
    raise ValueError(
        f"unrecognized canonical source family for {transition}: {source_value}"
    )


def require_same_source(actual: str | Path, expected: str | Path, *, context: str) -> None:
    if resolve_repository_path(actual) != resolve_repository_path(expected):
        raise ValueError(f"canonical source mismatch for {context}")


def canonical_identity(row: dict[str, Any]) -> tuple[str, float, float, str]:
    source = row["canonical_source"] if "canonical_source" in row else row["source"]
    return (
        str(row["transition"]),
        fvalue(row, "chiN"),
        fvalue(row, "fA"),
        relative(str(source)),
    )


def validate_identity_coverage(
    publication: Iterable[dict[str, Any]], generated: Iterable[dict[str, Any]]
) -> None:
    expected = [canonical_identity(row) for row in publication]
    actual = [canonical_identity(row) for row in generated]
    if len(expected) != len(set(expected)):
        raise ValueError("canonical publication identities are not unique")
    if len(actual) != len(set(actual)):
        raise ValueError("generated canonical identities are not unique")
    if set(actual) != set(expected):
        missing = set(expected) - set(actual)
        extra = set(actual) - set(expected)
        raise ValueError(
            "generated rows do not cover canonical identities one-to-one; "
            f"missing={missing}, extra={extra}"
        )


def point_row(
    publication: dict[str, str],
    *,
    contract_id: str,
    root_direction: str,
    root_coordinate: str,
    phase_a: str,
    phase_b: str,
    phase_a_grid: str,
    phase_b_grid: str,
    energy_evaluation: str,
    resolution_treatment: str,
    richardson_used: bool,
    resolution_shift_fa: float | None,
    root_semantics: str,
    uncertainty_scope: str,
    strict_total_coordinate_bound: bool,
    maximum_cell_stress: float | None,
    stress_policy: str,
    audit_status: str,
    evidence_paths: Iterable[Path],
) -> dict[str, Any]:
    width = fvalue(publication, "bracket_width")
    if width <= 0.0:
        raise ValueError(f"nonpositive bracket width in {publication}")
    source_paths, source_hashes = evidence(evidence_paths)
    return {
        "schema": SCHEMA,
        "transition": publication["transition"],
        "chiN": fvalue(publication, "chiN"),
        "fA": fvalue(publication, "fA"),
        "contract_id": contract_id,
        "root_direction": root_direction,
        "root_coordinate": root_coordinate,
        "bracket_width": width,
        "coordinate_half_width": width / 2.0,
        "phase_a": phase_a,
        "phase_b": phase_b,
        "phase_a_resolution": phase_a_grid,
        "phase_b_resolution": phase_b_grid,
        "energy_evaluation": energy_evaluation,
        "resolution_treatment": resolution_treatment,
        "richardson_used": str(richardson_used).lower(),
        "resolution_shift_fA": (
            "" if resolution_shift_fa is None else resolution_shift_fa
        ),
        "root_semantics": root_semantics,
        "uncertainty_scope": uncertainty_scope,
        "strict_total_coordinate_bound": str(strict_total_coordinate_bound).lower(),
        "maximum_abs_cell_stress": (
            "" if maximum_cell_stress is None else maximum_cell_stress
        ),
        "stress_tier_under_current_contract": stress_tier(maximum_cell_stress),
        "stress_policy": stress_policy,
        "audit_status": audit_status,
        "claim_kind": publication["claim_kind"],
        "canonical_status": publication["status"],
        "canonical_source": relative(publication["source"]),
        "evidence_sources": source_paths,
        "evidence_sha256": source_hashes,
    }


def build_pointwise_table() -> list[dict[str, Any]]:
    audited_transitions = {"S/DIS", "S/C", "G/C", "G/L"}
    publication = [
        row for row in read_rows(PUBLICATION_LEDGER)
        if row["transition"] in audited_transitions
    ]
    expected_counts = {"S/DIS": 33, "S/C": 27, "G/C": 22, "G/L": 24}
    actual_counts = {
        transition: sum(row["transition"] == transition for row in publication)
        for transition in sorted(audited_transitions)
    }
    if len(publication) != 106 or actual_counts != expected_counts:
        raise ValueError(
            "expected the 106-root core audit "
            f"{expected_counts}, found {len(publication)} rows {actual_counts}"
        )
    if any(
        row["status"] != "accepted" or row["claim_kind"] != "boundary_bracket"
        for row in publication
    ):
        raise ValueError("publication ledger contains a non-accepted boundary row")
    sdis_points = load_sdis_points()
    corrections = load_corrections()
    active_provenance = read_rows(ACTIVE_PROVENANCE)
    legacy_root_source = (RESULTS / "bvk2_gc_gl_boundary_campaign").resolve()
    fixed_sdis_source = (RESULTS / "bvk2_sdis_fixed_fa_chi_refinement").resolve()
    fixed_sc_root = RESULTS / "bvk2_sc_fixed_fa_chi_continuation"
    filtered_gc_path = RESULTS / "bvk2_ud_theta_gc_filtered_r6/boundary_points.csv"
    filtered_gc = keyed(read_rows(filtered_gc_path))
    filtered_gc_endpoints_path = (
        RESULTS / "bvk2_ud_theta_gc_filtered_r6/endpoint_pairs.csv"
    )
    filtered_gc_endpoints = read_rows(filtered_gc_endpoints_path)
    filtered_gc_manifest_path = RESULTS / "bvk2_ud_theta_gc_filtered_r6/manifest.csv"
    filtered_gc_manifest = read_rows(filtered_gc_manifest_path)
    filtered_gc_pairs = validate_filtered_gc_contract(
        list(filtered_gc.values()), filtered_gc_endpoints, filtered_gc_manifest
    )
    matched_gl_path = (
        RESULTS
        / "bvk2_ud_theta_gl_matched_spacing_audit/matched_boundary_promotions.csv"
    )
    matched_gl = keyed(read_rows(matched_gl_path))
    sc_high_path = RESULTS / "bvk2_ud_theta_sc_highchi/accepted_boundaries.csv"
    sc_high = keyed(read_rows(sc_high_path))
    gc_high_path = (
        RESULTS
        / "bvk2_ud_theta_g_pocket_highchi_ledger/accepted_gc_boundaries.csv"
    )
    gc_high = keyed(read_rows(gc_high_path))
    gc54_path = RESULTS / "bvk2_ud_theta_gc_chi54_direct/boundary_point.csv"
    gc54 = read_one(gc54_path)
    gl_high_path = RESULTS / "bvk2_ud_theta_gl_direct_fine/accepted_boundaries.csv"
    gl_high = keyed(read_rows(gl_high_path))
    direct_gl_sources = {
        resolve_repository_path(row["source"]) for row in gl_high.values()
    }
    fcc_o70_path = RESULTS / "bvk2_fcc_o70_boundary_campaign/boundary_brackets.csv"
    gl15 = unique_record(
        (
            row
            for row in read_rows(fcc_o70_path)
            if row["boundary"] == "G/L"
            and math.isclose(fvalue(row, "chiN"), 15.0, abs_tol=1.0e-10)
            and row["status"] == "accepted"
        ),
        context="accepted G/L chiN=15",
    )

    result: list[dict[str, Any]] = []
    matched_corrections: set[tuple[str, float, float, str]] = set()
    for row in publication:
        transition = row["transition"]
        chi_n = fvalue(row, "chiN")
        root_fa = fvalue(row, "fA")
        canonical_source = resolve_repository_path(row["source"])
        source_family = classify_canonical_source(
            transition,
            row["source"],
            direct_gl_sources=direct_gl_sources,
        )

        if source_family in {"sdis", "sdis_fixed_fa"}:
            if source_family == "sdis_fixed_fa":
                point, path = match_fixed_fa_point(
                    row,
                    fixed_sdis_source / "points",
                    "S_DIS_target_fA*.csv",
                )
            else:
                point, path = match_sdis(row, sdis_points)
            if not math.isclose(
                fvalue(point, "bracket_width"),
                fvalue(row, "bracket_width"),
                abs_tol=2.0e-12,
            ):
                raise ValueError("S/DIS bracket width disagrees with publication ledger")
            direction = point["direction"]
            result.append(
                point_row(
                    row,
                    contract_id=(
                        "SDIS_BCC48_ROOT_FA"
                        if direction == "fA"
                        else "SDIS_BCC48_ROOT_CHIN_REFINED"
                    ),
                    root_direction=direction,
                    root_coordinate="fA" if direction == "fA" else "chiN",
                    phase_a="BCC",
                    phase_b="DIS",
                    phase_a_grid="48x48x48",
                    phase_b_grid="analytic homogeneous",
                    energy_evaluation="direct BCC excess energy; F_DIS=0 analytic",
                    resolution_treatment="direct grid; no fine-grid correction",
                    richardson_used=False,
                    resolution_shift_fa=None,
                    root_semantics="ordered-branch survival midpoint",
                    uncertainty_scope=(
                        "predicate-bracket coordinate half-width; excludes BCC grid bias"
                    ),
                    strict_total_coordinate_bound=False,
                    maximum_cell_stress=None,
                    stress_policy=(
                        "stress-oriented interior cell root; family observed "
                        "|stress| <= 2.80e-6"
                    ),
                    audit_status="accepted_survival_root_resolution_bias_unquantified",
                    evidence_paths=[path, PUBLICATION_LEDGER],
                )
            )
            continue

        if source_family == "sc_fixed_fa":
            canonical_path = canonical_source
            point_directory = (
                canonical_path.parent
                if canonical_path.is_file()
                else fixed_sc_root / "points"
            )
            point_pattern = (
                canonical_path.name
                if canonical_path.is_file()
                else "S_C_target_fA*.csv"
            )
            point, path = match_fixed_fa_point(
                row,
                point_directory,
                point_pattern,
            )
            if not math.isclose(
                fvalue(point, "bracket_width"),
                fvalue(row, "bracket_width"),
                abs_tol=2.0e-12,
            ):
                raise ValueError("fixed-fA S/C bracket width disagrees with ledger")
            result.append(
                point_row(
                    row,
                    contract_id="SC_FIXED_FA_BCC48_CYL48x84",
                    root_direction="chiN",
                    root_coordinate="chiN",
                    phase_a="BCC",
                    phase_b="CYL",
                    phase_a_grid="48x48x48",
                    phase_b_grid="48x84",
                    energy_evaluation="signed BCC-CYL free-energy bracket",
                    resolution_treatment="direct grid; no fine-grid correction",
                    richardson_used=False,
                    resolution_shift_fa=None,
                    root_semantics="signed endpoint root at fixed fA",
                    uncertainty_scope=(
                        "accepted signed-bracket half-width; excludes spatial-grid bias"
                    ),
                    strict_total_coordinate_bound=False,
                    maximum_cell_stress=None,
                    stress_policy="phase-specific campaign gate; no uniform recertification",
                    audit_status="legacy_phase_specific_acceptance_not_uniformly_recertified",
                    evidence_paths=[path, PUBLICATION_LEDGER],
                )
            )
            continue

        if source_family == "richardson":
            key = canonical_identity(row)
            correction, correction_path = match_correction(
                row, corrections, active_provenance
            )
            corrected = fvalue(correction, "root_fA_corrected")
            if not math.isclose(corrected, root_fa, abs_tol=2.0e-10):
                raise ValueError(
                    f"active correction does not match publication root for {key}"
                )
            matched_corrections.add(key)
            fine_a = correction["fine_dims_a"]
            fine_b = correction["fine_dims_b"]
            coarse_a = infer_coarse_dims(
                fine_a,
                fvalue(correction, "dx_a_campaign_rg"),
                fvalue(correction, "dx_a_fine_rg"),
            )
            coarse_b = infer_coarse_dims(
                fine_b,
                fvalue(correction, "dx_b_campaign_rg"),
                fvalue(correction, "dx_b_fine_rg"),
            )
            point_name = (
                f"{transition.replace('/', '_')}_target_chiN"
                f"{chi_n:.6f}".replace(".", "p")
                + ".csv"
            )
            campaign_point = (
                RESULTS / "bvk2_gc_gl_boundary_campaign/points" / point_name
            )
            phase_a = correction["phase_a"]
            phase_b = correction["phase_b"]
            if (transition, phase_a, phase_b) not in {
                ("S/C", "BCC", "CYL"),
                ("G/C", "GYR", "CYL"),
            }:
                raise ValueError(f"wrong corrected competitors for {key}")
            shift = fvalue(correction, "shift")
            half_width = fvalue(row, "bracket_width") / 2.0
            result.append(
                point_row(
                    row,
                    contract_id=(
                        f"{transition.replace('/', '')}_RICHARDSON_"
                        f"{coarse_a}_TO_{fine_a}_{coarse_b}_TO_{fine_b}"
                    ),
                    root_direction="fA",
                    root_coordinate="fA",
                    phase_a=phase_a,
                    phase_b=phase_b,
                    phase_a_grid=f"{coarse_a} -> {fine_a}",
                    phase_b_grid=f"{coarse_b} -> {fine_b}",
                    energy_evaluation="continuum estimate E=E_inf-C*dx^2",
                    resolution_treatment="second-order Richardson on both phases",
                    richardson_used=True,
                    resolution_shift_fa=shift,
                    root_semantics="corrected coordinate from legacy raw bracket",
                    uncertainty_scope=(
                        "stored raw-bracket half-width only; Richardson-model "
                        "uncertainty is not quantified"
                    ),
                    strict_total_coordinate_bound=False,
                    maximum_cell_stress=None,
                    stress_policy="legacy campaign gate; no unified recertification",
                    audit_status=(
                        "resolution_shift_exceeds_stored_half_width"
                        if abs(shift) > half_width
                        else "resolution_shift_within_stored_half_width_but_total_bound_missing"
                    ),
                    evidence_paths=[
                        campaign_point,
                        correction_path,
                        ACTIVE_PROVENANCE,
                        PUBLICATION_LEDGER,
                    ],
                )
            )
            continue

        if source_family == "gc_filtered":
            source = filtered_gc[chi_n]
            if not (
                math.isclose(fvalue(source, "coordinate"), root_fa, abs_tol=2.0e-10)
                and math.isclose(
                    2.0 * fvalue(source, "coordinate_half_width"),
                    fvalue(row, "bracket_width"),
                    abs_tol=2.0e-12,
                )
            ):
                raise ValueError(f"filtered G/C source mismatch at {chi_n}")
            endpoint_rows = list(filtered_gc_pairs[chi_n])
            maximum = validate_current_stress(
                (
                    endpoint.get(key)
                    for endpoint in endpoint_rows
                    for key in ("gyr_cell_stress", "cyl_cell_stress")
                ),
                context=f"filtered G/C chiN={chi_n}",
            )
            result.append(
                point_row(
                    row,
                    contract_id="GC_FILTERED_R6_GYR96_CYL128x224",
                    root_direction="fA",
                    root_coordinate="fA",
                    phase_a="GYR",
                    phase_b="CYL",
                    phase_a_grid="96x96x96",
                    phase_b_grid="128x224",
                    energy_evaluation="same-schema filtered factor-6 endpoint energies",
                    resolution_treatment="direct filtered-sensor bracket",
                    richardson_used=False,
                    resolution_shift_fa=None,
                    root_semantics="strict signed endpoint secant",
                    uncertainty_scope=(
                        "accepted signed-bracket half-width; excludes spatial-grid bias"
                    ),
                    strict_total_coordinate_bound=False,
                    maximum_cell_stress=maximum,
                    stress_policy="later uniform audit: target 1e-3; hard maximum 2e-3",
                    audit_status="current_unified_direct_contract",
                    evidence_paths=[
                        filtered_gc_path,
                        filtered_gc_endpoints_path,
                        filtered_gc_manifest_path,
                        *(Path(endpoint["artifact"]) for endpoint in endpoint_rows),
                        PUBLICATION_LEDGER,
                    ],
                )
            )
            continue

        if source_family == "sc_high":
            source = sc_high[chi_n]
            if not (
                math.isclose(fvalue(source, "fA_secant"), root_fa, abs_tol=2.0e-10)
                and math.isclose(
                    fvalue(source, "bracket_width"),
                    fvalue(row, "bracket_width"),
                    abs_tol=2.0e-12,
                )
            ):
                raise ValueError(f"S/C high-chi source mismatch at {chi_n}")
            boundary_path = resolve_repository_path(source["source"])
            require_same_source(
                canonical_source, boundary_path, context=f"S/C chiN={chi_n}"
            )
            endpoint_path = boundary_path.parent / "endpoint_energies.csv"
            maximum = max_stress(
                endpoint_path,
                chi_n,
                (fvalue(source, "fA_lower"), fvalue(source, "fA_upper")),
                (("bcc_cell_stress", "cyl_cell_stress"),),
            )
            tier = stress_tier(maximum)
            if tier != "above_current_maximum":
                validate_current_stress(
                    (maximum,), context=f"S/C chiN={chi_n} pass label"
                )
            result.append(
                point_row(
                    row,
                    contract_id="SC_THETA_BCC48_CYL48x84_AUDIT4",
                    root_direction="fA",
                    root_coordinate="fA",
                    phase_a="BCC",
                    phase_b="CYL",
                    phase_a_grid="48x48x48",
                    phase_b_grid="48x84",
                    energy_evaluation="factor-3 optimization; signed factor-4 audit",
                    resolution_treatment="direct; no Richardson correction",
                    richardson_used=False,
                    resolution_shift_fa=None,
                    root_semantics="strict signed endpoint secant",
                    uncertainty_scope=(
                        "accepted signed-bracket half-width; excludes spatial-grid bias"
                    ),
                    strict_total_coordinate_bound=False,
                    maximum_cell_stress=maximum,
                    stress_policy="current target 1e-3; hard maximum 2e-3",
                    audit_status=(
                        "direct_root_current_stress_gate_pass"
                        if tier != "above_current_maximum"
                        else "legacy_accepted_endpoint_exceeds_current_stress_maximum"
                    ),
                    evidence_paths=[
                        sc_high_path,
                        boundary_path,
                        endpoint_path,
                        PUBLICATION_LEDGER,
                    ],
                )
            )
            continue

        if source_family in {"gc54", "gc_high"}:
            if source_family == "gc54":
                source = gc54
                source_path = gc54_path
                endpoint_path = gc54_path.parent / "endpoint_energies.csv"
                if not math.isclose(
                    fvalue(source, "fA"), root_fa, abs_tol=2.0e-10
                ):
                    raise ValueError("G/C chi54 source mismatch")
                maximum = max_stress(
                    endpoint_path,
                    chi_n,
                    (fvalue(source, "fA_lower"), fvalue(source, "fA_upper")),
                    (("gyr_stress", "cyl_stress"),),
                )
                validate_current_stress(
                    (maximum,), context=f"G/C chiN={chi_n} unified label"
                )
                grid_a = source["gyr_resolution"]
                grid_b = source["cyl_resolution"]
                audit = source["audit_factor"]
                contract_id = "GC_THETA_GYR112_CYL96x168_AUDIT3"
                status = "current_unified_direct_contract"
                paths = [source_path, endpoint_path, PUBLICATION_LEDGER]
            else:
                source = gc_high[chi_n]
                source_path = resolve_repository_path(source["source"])
                endpoint_path = resolve_repository_path(source["contract_source"])
                if not math.isclose(
                    fvalue(source, "fA_secant"), root_fa, abs_tol=2.0e-10
                ):
                    raise ValueError(f"G/C high-chi source mismatch at {chi_n}")
                maximum = max_stress(
                    endpoint_path,
                    chi_n,
                    (fvalue(source, "fA_lower"), fvalue(source, "fA_upper")),
                    (
                        ("gyr_cell_stress", "cyl_cell_stress"),
                        ("gyr_stress", "cyl_stress"),
                    ),
                )
                grid_a = source["gyr_resolution"]
                grid_b = source["cyl_resolution"]
                audit = source["audit_factor"]
                contract_id = "GC_THETA_GYR96_CYL48x84_AUDIT3"
                status = (
                    "legacy_direct_current_stress_gate_pass"
                    if stress_tier(maximum) != "above_current_maximum"
                    else "legacy_accepted_endpoint_exceeds_current_stress_maximum"
                )
                if status == "legacy_direct_current_stress_gate_pass":
                    validate_current_stress(
                        (maximum,), context=f"G/C chiN={chi_n} pass label"
                    )
                paths = [
                    gc_high_path,
                    source_path,
                    endpoint_path,
                    PUBLICATION_LEDGER,
                ]
            result.append(
                point_row(
                    row,
                    contract_id=contract_id,
                    root_direction="fA",
                    root_coordinate="fA",
                    phase_a="GYR",
                    phase_b="CYL",
                    phase_a_grid=grid_a,
                    phase_b_grid=grid_b,
                    energy_evaluation=f"factor-2 optimization; signed factor-{audit} audit",
                    resolution_treatment="direct; no Richardson correction",
                    richardson_used=False,
                    resolution_shift_fa=None,
                    root_semantics="strict signed endpoint secant",
                    uncertainty_scope=(
                        "accepted signed-bracket half-width; excludes spatial-grid bias"
                    ),
                    strict_total_coordinate_bound=False,
                    maximum_cell_stress=maximum,
                    stress_policy="current target 1e-3; hard maximum 2e-3",
                    audit_status=status,
                    evidence_paths=paths,
                )
            )
            continue

        if source_family in {"gl_matched", "gl15", "gl_direct"}:
            if source_family == "gl_matched":
                source = matched_gl[chi_n]
                if not (
                    math.isclose(
                        fvalue(source, "fA_secant"), root_fa, abs_tol=2.0e-10
                    )
                    and math.isclose(
                        fvalue(source, "bracket_width"),
                        fvalue(row, "bracket_width"),
                        abs_tol=2.0e-12,
                    )
                ):
                    raise ValueError(f"matched-spacing G/L source mismatch at {chi_n}")
                endpoint_paths = [
                    resolve_repository_path(source["lower_endpoint"]),
                    resolve_repository_path(source["upper_endpoint"]),
                ]
                endpoint_rows = exact_coordinate_rows(
                    [read_one(path) for path in endpoint_paths],
                    coordinates=(
                        fvalue(source, "fA_lower"),
                        fvalue(source, "fA_upper"),
                    ),
                    coordinate_key="fA",
                    per_coordinate=1,
                    context=f"matched-spacing G/L chiN={chi_n}",
                )
                maximum = validate_current_stress(
                    (
                        endpoint.get(key)
                        for endpoint in endpoint_rows
                        for key in ("gyr_stress", "matched_lam_stress")
                    ),
                    context=f"matched-spacing G/L chiN={chi_n}",
                )
                result.append(
                    point_row(
                        row,
                        contract_id=(
                            f"GL_MATCHED_GYR{source['gyr_resolution']}_"
                            f"LAM{source['lam_resolution']}_AUDIT3"
                        ),
                        root_direction="fA",
                        root_coordinate="fA",
                        phase_a="GYR",
                        phase_b="LAM",
                        phase_a_grid=f"{source['gyr_resolution']}x" * 2
                        + source["gyr_resolution"],
                        phase_b_grid=source["lam_resolution"],
                        energy_evaluation=source["energy_contract"],
                        resolution_treatment="matched-spacing direct factor-3 bracket",
                        richardson_used=False,
                        resolution_shift_fa=None,
                        root_semantics="strict signed endpoint secant",
                        uncertainty_scope=(
                            "accepted signed-bracket half-width; excludes residual "
                            "matched-spacing bias"
                        ),
                        strict_total_coordinate_bound=False,
                        maximum_cell_stress=maximum,
                        stress_policy="later uniform audit: target 1e-3; hard maximum 2e-3",
                        audit_status="current_unified_matched_spacing_contract",
                        evidence_paths=[
                            matched_gl_path,
                            *endpoint_paths,
                            PUBLICATION_LEDGER,
                        ],
                    )
                )
                continue

            if source_family == "gl15":
                if not (
                    math.isclose(
                        fvalue(gl15, "fA_estimate"), root_fa, abs_tol=2.0e-10
                    )
                    and math.isclose(
                        fvalue(gl15, "width"),
                        fvalue(row, "bracket_width"),
                        abs_tol=2.0e-12,
                    )
                ):
                    raise ValueError("G/L chi15 source mismatch")
                endpoint_paths = [
                    resolve_repository_path(gl15[name])
                    for name in (
                        "lower_phaseA_endpoint",
                        "lower_phaseB_endpoint",
                        "upper_phaseA_endpoint",
                        "upper_phaseB_endpoint",
                    )
                ]
                endpoint_records = list(
                    zip(
                        [read_one(path) for path in endpoint_paths],
                        ("GYR", "LAM", "GYR", "LAM"),
                        strict=True,
                    )
                )
                endpoint_rows = exact_coordinate_rows(
                    [endpoint for endpoint, _phase in endpoint_records],
                    coordinates=(fvalue(gl15, "fA_lower"), fvalue(gl15, "fA_upper")),
                    coordinate_key="fA",
                    per_coordinate=2,
                    context="G/L chiN=15 phase endpoints",
                )
                for coordinate in (fvalue(gl15, "fA_lower"), fvalue(gl15, "fA_upper")):
                    phases = {
                        phase
                        for endpoint, phase in endpoint_records
                        if math.isclose(
                            fvalue(endpoint, "fA"), coordinate, abs_tol=2.0e-10
                        )
                    }
                    if phases != {"GYR", "LAM"}:
                        raise ValueError(
                            f"G/L chiN=15 phase endpoints are not unique at {coordinate}"
                        )
                maximum = validate_current_stress(
                    (endpoint.get("cell_stress") for endpoint in endpoint_rows),
                    context="G/L chiN=15 current label",
                )
                grid_a = "64x64x64"
                contract_id = "GL_THETA_GYR64_LAM1024_AUDIT3"
                paths = [fcc_o70_path, *endpoint_paths, PUBLICATION_LEDGER]
            else:
                source = gl_high[chi_n]
                if not (
                    math.isclose(
                        fvalue(source, "fA_secant"), root_fa, abs_tol=2.0e-10
                    )
                    and source["status"] == "accepted"
                    and math.isclose(
                        fvalue(source, "bracket_width"),
                        fvalue(row, "bracket_width"),
                        abs_tol=2.0e-12,
                    )
                ):
                    raise ValueError(f"G/L high-chi source mismatch at {chi_n}")
                require_same_source(
                    canonical_source,
                    source["source"],
                    context=f"G/L chiN={chi_n}",
                )
                maximum = validate_current_stress(
                    (
                        source.get("gyr_cell_stress_max"),
                        source.get("lam_cell_stress_max"),
                    ),
                    context=f"G/L chiN={chi_n} current label",
                )
                grid_a = source["gyr_resolution"]
                contract_id = "GL_THETA_GYR112_LAM1024_AUDIT3"
                paths = [
                    gl_high_path,
                    resolve_repository_path(source["source"]),
                    PUBLICATION_LEDGER,
                ]
            result.append(
                point_row(
                    row,
                    contract_id=contract_id,
                    root_direction="fA",
                    root_coordinate="fA",
                    phase_a="GYR",
                    phase_b="LAM",
                    phase_a_grid=grid_a,
                    phase_b_grid="1024",
                    energy_evaluation="factor-2 optimization; signed factor-3 audit",
                    resolution_treatment="direct; no Richardson correction",
                    richardson_used=False,
                    resolution_shift_fa=None,
                    root_semantics="strict signed endpoint secant",
                    uncertainty_scope=(
                        "accepted signed-bracket half-width; excludes residual "
                        "spatial-grid bias"
                    ),
                    strict_total_coordinate_bound=False,
                    maximum_cell_stress=maximum,
                    stress_policy="target 1e-3; accepted fallback through 2e-3",
                    audit_status="current_direct_contract",
                    evidence_paths=paths,
                )
            )
            continue

        raise ValueError(f"unhandled publication row {row}")

    expected_corrections = {
        canonical_identity(row)
        for row in publication
        if row["transition"] in {"S/C", "G/C"}
        and resolve_repository_path(row["source"]) == legacy_root_source
    }
    if matched_corrections != expected_corrections:
        raise ValueError("not every canonical corrected row was matched")
    validate_identity_coverage(publication, result)
    return result


def summarize_contracts(points: list[dict[str, Any]]) -> list[dict[str, Any]]:
    groups: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in points:
        groups[row["contract_id"]].append(row)

    summaries: list[dict[str, Any]] = []
    for contract_id, rows in sorted(groups.items()):
        rows.sort(key=lambda row: (row["chiN"], row["fA"]))
        invariant_fields = (
            "transition",
            "root_direction",
            "root_coordinate",
            "phase_a",
            "phase_b",
            "phase_a_resolution",
            "phase_b_resolution",
            "energy_evaluation",
            "resolution_treatment",
            "richardson_used",
            "root_semantics",
            "uncertainty_scope",
            "strict_total_coordinate_bound",
            "stress_policy",
            "claim_kind",
            "canonical_status",
        )
        heterogeneous = [
            field
            for field in invariant_fields
            if len({str(row[field]) for row in rows}) != 1
        ]
        if heterogeneous:
            raise ValueError(
                f"contract is not homogeneous for {contract_id}: {heterogeneous}"
            )
        uncertainties = [float(row["coordinate_half_width"]) for row in rows]
        stresses = [
            float(row["maximum_abs_cell_stress"])
            for row in rows
            if row["maximum_abs_cell_stress"] != ""
        ]
        statuses = sorted({row["audit_status"] for row in rows})
        sources = sorted(
            {
                source
                for row in rows
                for source in str(row["evidence_sources"]).split(";")
            }
        )
        summaries.append(
            {
                "schema": SCHEMA,
                "contract_id": contract_id,
                "transition": rows[0]["transition"],
                "point_count": len(rows),
                "chiN_min": min(row["chiN"] for row in rows),
                "chiN_max": max(row["chiN"] for row in rows),
                "fA_min": min(row["fA"] for row in rows),
                "fA_max": max(row["fA"] for row in rows),
                "root_direction": rows[0]["root_direction"],
                "phase_a": rows[0]["phase_a"],
                "phase_b": rows[0]["phase_b"],
                "phase_a_resolution": rows[0]["phase_a_resolution"],
                "phase_b_resolution": rows[0]["phase_b_resolution"],
                "energy_evaluation": rows[0]["energy_evaluation"],
                "resolution_treatment": rows[0]["resolution_treatment"],
                "coordinate_half_width_min": min(uncertainties),
                "coordinate_half_width_max": max(uncertainties),
                "maximum_abs_cell_stress": max(stresses) if stresses else "",
                "uncertainty_scope": rows[0]["uncertainty_scope"],
                "audit_statuses": ";".join(statuses),
                "evidence_sources": ";".join(sources),
            }
        )
    return summaries


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty table {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def _contract_kind(row: dict[str, Any]) -> str:
    if row["transition"] == "S/DIS":
        return "survival"
    if row["richardson_used"] == "true":
        return "richardson"
    return "direct"


def _matplotlib_pyplot():
    import matplotlib

    matplotlib.use("Agg")
    matplotlib.rcParams.update(
        {
            "svg.fonttype": "none",
            "svg.hashsalt": "bvk2-resolution-audit-v1",
        }
    )
    import matplotlib.pyplot as plt

    return plt


def write_resolution_contract_svg(
    path: Path,
    points: list[dict[str, Any]],
    contracts: list[dict[str, Any]],
) -> None:
    """Render the accepted-root resolution map with Matplotlib."""
    with mpltex.acs_decorator:
        _write_resolution_contract_svg_acs(path, points, contracts)


def _write_resolution_contract_svg_acs(
    path: Path,
    points: list[dict[str, Any]],
    contracts: list[dict[str, Any]],
) -> None:
    plt = _matplotlib_pyplot()
    from matplotlib.lines import Line2D

    colors = dict(
        zip(("S/DIS", "S/C", "G/C", "G/L"), mpltex.tableau_10[:4])
    )
    markers = {"survival": "o", "richardson": "^", "direct": "s"}
    figure, (phase_ax, contract_ax) = plt.subplots(
        1,
        2,
        figsize=(7.0, 4.5),
        gridspec_kw={"width_ratios": (1.08, 1.0)},
        constrained_layout=True,
    )
    for transition in ("S/DIS", "S/C", "G/C", "G/L"):
        for kind in ("survival", "richardson", "direct"):
            selected = [
                row
                for row in points
                if row["transition"] == transition and _contract_kind(row) == kind
            ]
            if selected:
                phase_ax.scatter(
                    [float(row["fA"]) for row in selected],
                    [float(row["chiN"]) for row in selected],
                    s=28,
                    marker=markers[kind],
                    color=colors[transition],
                    edgecolor="white",
                    linewidth=0.45,
                    zorder=3,
                )
    phase_ax.set(
        xlabel=r"A-block fraction, $f_A$",
        ylabel=r"$\chi N$",
        xlim=(0.14, 0.505),
        ylim=(9, 62),
    )
    phase_ax.text(
        -0.11,
        1.02,
        "(a)",
        transform=phase_ax.transAxes,
        fontsize=11,
        fontweight="bold",
    )

    ordered = sorted(
        contracts,
        key=lambda row: (row["transition"], row["chiN_min"], row["contract_id"]),
    )
    for index, row in enumerate(ordered):
        transition = row["transition"]
        lo = float(row["chiN_min"])
        hi = float(row["chiN_max"])
        if math.isclose(lo, hi, abs_tol=1.0e-12):
            contract_ax.scatter(
                [lo], [index], s=35, color=colors[transition], zorder=3
            )
        else:
            contract_ax.plot(
                [lo, hi],
                [index, index],
                color=colors[transition],
                linewidth=6,
                solid_capstyle="round",
            )
    contract_ax.set(
        xlabel=r"$\chi N$",
        xlim=(9, 62),
        yticks=range(len(ordered)),
        yticklabels=[
            f"{row['transition']}  ({row['point_count']} pt)" for row in ordered
        ],
    )
    contract_ax.set_ylim(len(ordered) - 0.5, -0.5)
    contract_ax.text(
        -0.02,
        1.02,
        "(b)",
        transform=contract_ax.transAxes,
        fontsize=11,
        fontweight="bold",
    )

    phase_handles = [
        Line2D([0], [0], marker="o", linestyle="", color=colors[name], label=name)
        for name in ("S/DIS", "S/C", "G/C", "G/L")
    ]
    treatment_handles = [
        Line2D(
            [0],
            [0],
            marker=markers[kind],
            linestyle="",
            markerfacecolor=mpltex.almost_black,
            markeredgecolor=mpltex.almost_black,
            label=label,
        )
        for kind, label in (
            ("survival", "survival/bifurcation"),
            ("richardson", "Richardson corrected"),
            ("direct", "direct grid"),
        )
    ]
    phase_ax.legend(
        handles=phase_handles + treatment_handles,
        ncol=2,
        loc="upper right",
        frameon=False,
    )
    figure.savefig(
        path,
        format="svg",
        dpi=160,
        metadata={"Date": None, "Creator": "Matplotlib + mpltex ACS via uv"},
    )
    plt.close(figure)
    path.write_text(
        "\n".join(line.rstrip() for line in path.read_text(encoding="utf-8").splitlines())
        + "\n",
        encoding="utf-8",
    )


def write_acceptance_audit_svg(
    path: Path,
    points: list[dict[str, Any]],
    rejected: list[dict[str, str]],
) -> None:
    """Render accepted-contract and rejected-row audits with Matplotlib."""
    with mpltex.acs_decorator:
        _write_acceptance_audit_svg_acs(path, points, rejected)


def _write_acceptance_audit_svg_acs(
    path: Path,
    points: list[dict[str, Any]],
    rejected: list[dict[str, str]],
) -> None:
    plt = _matplotlib_pyplot()
    from matplotlib.lines import Line2D

    colors = dict(
        zip(("S/DIS", "S/C", "G/C", "G/L"), mpltex.tableau_10[:4])
    )
    labels = {
        "accepted_survival_root_resolution_bias_unquantified": "S/DIS survival; grid bias unresolved",
        "current_direct_contract": "current direct G/L contract",
        "resolution_shift_exceeds_stored_half_width": "Richardson shift above bracket half-width",
        "resolution_shift_within_stored_half_width_but_total_bound_missing": r"Richardson shift $\leq$ bracket half-width",
        "legacy_accepted_endpoint_exceeds_current_stress_maximum": "legacy endpoint above uniform stress max",
        "direct_root_current_stress_gate_pass": "direct root; current stress gate pass",
        "legacy_direct_current_stress_gate_pass": "legacy direct; current stress gate pass",
        "current_unified_direct_contract": "current unified direct contract",
        "current_unified_matched_spacing_contract": "current matched-spacing contract",
        "legacy_phase_specific_acceptance_not_uniformly_recertified": (
            "phase-specific; not uniformly recertified"
        ),
    }
    statuses = sorted(
        {row["audit_status"] for row in points},
        key=lambda status: (
            -sum(row["audit_status"] == status for row in points),
            status,
        ),
    )
    figure, (accepted_ax, rejected_ax) = plt.subplots(
        1,
        2,
        figsize=(7.0, 4.8),
        gridspec_kw={"width_ratios": (1.62, 1.0)},
        constrained_layout=True,
    )
    y_positions = list(range(len(statuses)))
    left = [0] * len(statuses)
    for transition in ("S/DIS", "S/C", "G/C", "G/L"):
        counts = [
            sum(
                row["audit_status"] == status
                and row["transition"] == transition
                for row in points
            )
            for status in statuses
        ]
        accepted_ax.barh(
            y_positions,
            counts,
            left=left,
            color=colors[transition],
            height=0.68,
            label=transition,
        )
        left = [a + b for a, b in zip(left, counts)]
    accepted_ax.set(
        xlabel="number of accepted points",
        yticks=y_positions,
        yticklabels=[labels[status] for status in statuses],
    )
    accepted_ax.text(
        -0.02,
        1.02,
        "(a)",
        transform=accepted_ax.transAxes,
        fontsize=11,
        fontweight="bold",
    )
    accepted_ax.invert_yaxis()
    accepted_ax.legend(ncol=2, loc="lower right", frameon=False)

    reason_marker = {
        "status=nonfinite_iteration": "X",
        "status=no_valid_sign_bracket": "x",
    }
    for row in rejected:
        reason = row["gate_reason"]
        rejected_ax.scatter(
            float(row["chiN"]),
            float(row["bracket_width"]),
            s=80,
            marker=reason_marker[reason],
            color=colors[row["transition"]],
            linewidth=1.8,
            zorder=3,
        )
        chi_n = float(row["chiN"])
        if chi_n == 55.0:
            offset, alignment = (-5, 8), "right"
        elif chi_n == 58.0:
            offset, alignment = (5, -13), "left"
        else:
            offset, alignment = (5, 3), "left"
        rejected_ax.annotate(
            f"{row['transition']}  {float(row['chiN']):g}",
            (float(row["chiN"]), float(row["bracket_width"])),
            xytext=offset,
            textcoords="offset points",
            color=colors[row["transition"]],
            fontsize=8,
            ha=alignment,
        )
    rejected_ax.set(
        xlabel=r"$\chi N$",
        ylabel="attempted bracket width",
        xlim=(19, 62),
        ylim=(0, 0.08),
    )
    rejected_ax.text(
        -0.02,
        1.02,
        "(b)",
        transform=rejected_ax.transAxes,
        fontsize=11,
        fontweight="bold",
    )
    rejected_ax.legend(
        handles=[
            Line2D(
                [0],
                [0],
                marker=marker,
                linestyle="",
                color=mpltex.almost_black,
                label=reason.removeprefix("status=").replace("_", " "),
            )
            for reason, marker in reason_marker.items()
        ],
        loc="upper left",
        frameon=False,
    )
    figure.savefig(
        path,
        format="svg",
        dpi=160,
        metadata={"Date": None, "Creator": "Matplotlib + mpltex ACS via uv"},
    )
    plt.close(figure)
    path.write_text(
        "\n".join(line.rstrip() for line in path.read_text(encoding="utf-8").splitlines())
        + "\n",
        encoding="utf-8",
    )


def write_readme(
    path: Path,
    points: list[dict[str, Any]],
    contracts: list[dict[str, Any]],
) -> None:
    counts = defaultdict(int)
    for row in points:
        counts[row["transition"]] += 1
    corrected = [row for row in points if row["richardson_used"] == "true"]
    corrected_shift_exceeds = [
        row
        for row in corrected
        if row["audit_status"] == "resolution_shift_exceeds_stored_half_width"
    ]
    stress_failures = [
        row
        for row in points
        if row["stress_tier_under_current_contract"] == "above_current_maximum"
    ]
    fixed_fa_sdis = [
        row
        for row in points
        if row["transition"] == "S/DIS" and row["root_direction"] == "chiN"
    ]
    corrected_counts = {
        transition: sum(row["transition"] == transition for row in corrected)
        for transition in ("S/C", "G/C")
    }
    shifted_counts = {
        transition: sum(
            row["transition"] == transition for row in corrected_shift_exceeds
        )
        for transition in ("S/C", "G/C")
    }

    table = [
        (
            "| Boundary / contract | Points | Resolutions | Energy treatment | "
            "Coordinate half-width | Audit status |"
        ),
        "|---|---:|---|---|---:|---|",
    ]
    for row in contracts:
        uncertainty = (
            f"{row['coordinate_half_width_min']:.3g}"
            if math.isclose(
                float(row["coordinate_half_width_min"]),
                float(row["coordinate_half_width_max"]),
                abs_tol=1.0e-15,
            )
            else (
                f"{row['coordinate_half_width_min']:.3g}–"
                f"{row['coordinate_half_width_max']:.3g}"
            )
        )
        table.append(
            f"| `{row['contract_id']}` | {row['point_count']} | "
            f"{row['phase_a']} {row['phase_a_resolution']}; "
            f"{row['phase_b']} {row['phase_b_resolution']} | "
            f"{row['resolution_treatment']} | {uncertainty} | "
            f"{row['audit_statuses']} |"
        )

    text = f"""# BVK2 phase-boundary resolution and uncertainty audit

This packet freezes the numerical metadata for all **{len(points)} canonical accepted
plot rows**: {counts['S/DIS']} S/DIS, {counts['S/C']} S/C,
{counts['G/C']} G/C, and {counts['G/L']} G/L. It is generated from
`results/bvk2_publication_phase_diagram/accepted_phase_boundaries.csv` and the
exact campaign/correction ledgers named in each row.

{chr(10).join(table)}

## What the uncertainty column means

`coordinate_half_width = bracket_width/2`. For fixed-`chiN` roots it is in
`fA`; for the {len(fixed_fa_sdis)} near-critical fixed-`fA` S/DIS roots it is in `chiN`.
It is **not** automatically a total continuum uncertainty:

- S/DIS is a survival/bifurcation bracket. Its half-width excludes BCC `48^3`
  spatial bias.
- Direct S/C, G/C, and G/L rows use signed endpoint brackets, but their
  half-widths exclude remaining spatial-grid bias.
- The {len(corrected)} legacy Richardson-corrected S/C and G/C rows retain the raw numerical
  bracket width while their plotted coordinate is shifted by a continuum
  extrapolation. No uncertainty of that extrapolation is tabulated.

Of the {len(corrected)} corrected rows ({corrected_counts['S/C']} S/C and
{corrected_counts['G/C']} G/C), **{len(corrected_shift_exceeds)}** have a resolution
shift larger than the stored half-width: {shifted_counts['G/C']} G/C and
{shifted_counts['S/C']} S/C rows. These rows therefore do not possess a strict
total coordinate bound under the current evidence.

## Unified-gate audit

The canonical plot remains unchanged by this read-only audit. The original
campaigns used phase-specific acceptance criteria. Applying the later uniform
stress audit retrospectively identifies **{len(stress_failures)}** direct
S/C/G/C rows with at least one endpoint above `2e-3`. This comparison does not
by itself invalidate their original phase-specific acceptance, but it prevents
those rows from being described as uniformly recertified under the later audit.
They are identified point-by-point in
`phase_boundary_resolution_uncertainty.csv` with
`audit_status=legacy_accepted_endpoint_exceeds_current_stress_maximum`.

All {counts['G/L']} G/L rows use direct factor-3 energies. The matched-spacing
subset uses the GYR and LAM grids recorded point by point; the remaining
direct rows use LAM `1024`. Their audited endpoints remain below `2e-3`.
Filtered G/C uses GYR `96^3` and CYL `128x224`; G/C at `chiN=54` uses the
direct GYR `112^3` / CYL `96x168` contract.

Consequently, this table completes the metadata-freezing work package but does
**not** close the manuscript's numerical-certification gate. The unresolved
items are total resolution uncertainty for legacy corrected S/C/G/C and
unified stress recertification for the flagged direct S/C/G/C slices.

## Files

- `phase_boundary_resolution_uncertainty.csv`: all {len(points)} canonical points.
- `resolution_contracts.csv`: grouped resolution/solver contracts.
- `resolution_contract_map.svg`: pointwise and contract-range visualization.
- `acceptance_rejection_audit.svg`: accepted-contract and rejected-row audit.
- `README.md`: this interpretation and reproduction record.

All `canonical_source` and `evidence_sources` entries are repository-relative;
the generator resolves both these portable paths and legacy absolute paths.

Reproduce from the `MonteCarlo` directory:

```bash
uv run python scripts/write_bvk2_resolution_uncertainty_table.py
uv run python test/test_write_bvk2_resolution_uncertainty_table.py
```

Schema: `{SCHEMA}`.
"""
    path.write_text(text, encoding="utf-8")


def write_outputs(output_dir: Path) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    points = build_pointwise_table()
    contracts = summarize_contracts(points)
    rejected = read_rows(
        RESULTS / "bvk2_publication_phase_diagram/rejected_points.csv"
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(output_dir / "phase_boundary_resolution_uncertainty.csv", points)
    write_csv(output_dir / "resolution_contracts.csv", contracts)
    write_resolution_contract_svg(
        output_dir / "resolution_contract_map.svg", points, contracts
    )
    write_acceptance_audit_svg(
        output_dir / "acceptance_rejection_audit.svg", points, rejected
    )
    write_readme(output_dir / "README.md", points, contracts)
    return points, contracts


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="destination directory",
    )
    args = parser.parse_args()
    points, contracts = write_outputs(args.output_dir.resolve())
    print(
        f"wrote {len(points)} canonical points across {len(contracts)} "
        f"resolution contracts to {args.output_dir.resolve()}"
    )


if __name__ == "__main__":
    main()
