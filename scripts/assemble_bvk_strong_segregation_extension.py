#!/usr/bin/env python3
"""Validate and assemble the integer chiN=36,...,50 lamellar extension."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import tempfile
from collections import defaultdict
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
STATE_SET = {
    (f, float(chi))
    for f in (0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50)
    for chi in range(36, 51)
}
MAIN_MODELS = (
    "scft",
    "ohta_kawasaki",
    "uneyama_doi",
    "liu2019_opf",
    "burp_ti",
    "bvk2",
)
AGGREGATE_MODELS = (
    "ohta_kawasaki",
    "uneyama_doi",
    "liu2019_opf",
    "burp_ti",
    "bvk2_fixed",
    "bvk2",
)
FIXED_SCHEMA = "bvk2-fixed-stiffness-figure6-root-v1"
OUTPUT_SCHEMA = "bvk-strong-segregation-extension-v1"
SAMPLE_COUNT = 1024
LABELS = {
    "liu2019_opf": "OPF",
    "ohta_kawasaki": "OK",
    "uneyama_doi": "UD",
    "burp_ti": "BURP",
    "bvk2_fixed": "QCE",
    "bvk2": "AQCE",
}


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def truth(value: str) -> bool:
    return value.strip().lower() == "true"


def atomic_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w", newline="", encoding="utf-8", dir=path.parent, delete=False
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(rows[0]),
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
        temporary = Path(handle.name)
    os.replace(temporary, path)


def atomic_text(path: Path, value: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", dir=path.parent, delete=False
    ) as handle:
        handle.write(value)
        temporary = Path(handle.name)
    os.replace(temporary, path)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def periodic_resample(values: list[float], count: int = SAMPLE_COUNT) -> np.ndarray:
    source = np.asarray(values, dtype=float)
    coordinate = np.arange(count, dtype=float) * len(source) / count
    left = np.floor(coordinate).astype(int)
    fraction = coordinate - left
    return ((1.0 - fraction) * source[left % len(source)]
            + fraction * source[(left + 1) % len(source)])


def align(reference: list[float], model: list[float], allow_swap: bool) -> tuple[np.ndarray, float]:
    target = periodic_resample(reference)
    base = periodic_resample(model)
    best: tuple[float, np.ndarray] | None = None
    candidates = (base, 1.0 - base) if allow_swap else (base,)
    for candidate in candidates:
        for shift in range(SAMPLE_COUNT):
            shifted = np.roll(candidate, -shift)
            rms = float(np.sqrt(np.mean((shifted - target) ** 2)))
            if best is None or rms < best[0]:
                best = (rms, shifted)
    assert best is not None
    return best[1], best[0]


def ok_case_dir(ok_root: Path, f: float, chi: float) -> Path:
    branch = f"f{round(100 * f):03d}"
    return ok_root / branch / "cases" / f"f{f:g}_chiN{chi:g}"


def assemble(root: Path, ok_root: Path) -> None:
    main: list[dict[str, str]] = []
    profiles: list[dict[str, str]] = []
    native: list[dict[str, str]] = []
    required_main_gates = (
        "field_gate_pass",
        "cell_gate_pass",
        "resolution_gate_pass",
        "composition_gate_pass",
        "morphology_gate_pass",
        "period_local_minimum_check_pass",
    )
    for f, chi in sorted(STATE_SET):
        case_id = f"f{f:g}_chiN{chi:g}"
        case_dir = root / "cases" / case_id
        case_summary = read_rows(case_dir / "summary.csv")
        case_profiles = read_rows(case_dir / "profiles.csv")
        case_native = read_rows(case_dir / "native_profiles.csv")
        ok_dir = ok_case_dir(ok_root, f, chi)
        ok_summary = [
            row for row in read_rows(ok_dir / "summary.csv")
            if row["model"] == "ohta_kawasaki"
        ]
        ok_profiles = [
            row for row in read_rows(ok_dir / "profiles.csv")
            if row["model"] == "ohta_kawasaki"
        ]
        ok_native = [
            row for row in read_rows(ok_dir / "native_profiles.csv")
            if row["model"] == "ohta_kawasaki"
        ]
        if (
            len(ok_summary) != 1
            or len(ok_profiles) != 256
            or len(ok_native) != 128
            or ok_summary[0]["calibration_role"]
            != "liu2019_ok_asymptotic_quadratic_extrapolated_c3_c4"
        ):
            raise ValueError(f"{case_id}: extrapolated OK inventory is incomplete")
        case_summary.extend(ok_summary)
        case_profiles.extend(ok_profiles)
        case_native.extend(ok_native)
        case_models = {row["model"] for row in case_summary}
        if len(case_summary) != len(MAIN_MODELS) or case_models != set(MAIN_MODELS):
            raise ValueError(f"{case_id}: main-model inventory is incomplete")
        if len(case_profiles) != len(MAIN_MODELS) * 256:
            raise ValueError(f"{case_id}: aligned-profile inventory is incomplete")
        if any(
            row["case_id"] != case_id
            or float(row["f"]) != f
            or float(row["chiN"]) != chi
            or row["status"] != "accepted"
            or row["period_boundary_limited"].lower() == "true"
            or not all(truth(row[name]) for name in required_main_gates)
            for row in case_summary
        ):
            raise ValueError(f"{case_id}: main-model gates failed")
        expected_native = {
            row["model"]: int(float(row["grid_count"])) for row in case_summary
        }
        observed_native: dict[str, int] = defaultdict(int)
        for row in case_native:
            observed_native[row["model"]] += 1
        if dict(observed_native) != expected_native:
            raise ValueError(f"{case_id}: native-profile inventory is incomplete")
        main.extend(case_summary)
        profiles.extend(case_profiles)
        native.extend(case_native)

    keys = {(float(row["f"]), float(row["chiN"]), row["model"]) for row in main}
    expected = {(f, chi, model) for f, chi in STATE_SET for model in MAIN_MODELS}
    if len(main) != len(expected) or keys != expected:
        raise ValueError(f"main extension inventory is incomplete: {len(main)}/{len(expected)}")
    scft_native: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in native:
        if row["model"] == "scft":
            scft_native[row["case_id"]].append(row)

    fixed_summary: list[dict[str, object]] = []
    fixed_profiles: list[dict[str, object]] = []
    for f, chi in sorted(STATE_SET):
        case_id = f"f{f:g}_chiN{chi:g}"
        case_dir = root / "fixed_cases" / case_id
        roots = read_rows(case_dir / "root.csv")
        fields = read_rows(case_dir / "profile.csv")
        if len(roots) != 1 or len(fields) != 256:
            raise ValueError(f"{case_id}: incomplete fixed-stiffness artifacts")
        row = roots[0]
        required_gates = ("field_gate_pass", "composition_gate_pass",
                          "morphology_gate_pass", "cell_gate_pass",
                          "local_minimum_check_pass", "stress_orientation_pass",
                          "accepted")
        valid = (row["schema"] == FIXED_SCHEMA and row["status"] == "accepted"
                 and row["model"] == "bvk2_fixed"
                 and row["adaptive"].lower() == "false"
                 and row["case_id"] == case_id
                 and float(row["f"]) == f
                 and float(row["chiN"]) == chi
                 and int(row["nx"]) == 256
                 and int(row["oversample_factor"]) == 2
                 and float(row["sensor_filter_ratio"]) == 12.0
                 and all(truth(row[name]) for name in required_gates)
                 and abs(float(row["root_cell_stress"])) <= 1.0e-6
                 and float(row["projected_force_norm"]) <= 1.0e-5
                 and float(row["projected_force_maxabs"]) <= 2.0e-5)
        if not valid:
            raise ValueError(f"{case_id}: fixed-stiffness gates failed")
        source_summary = ROOT / row["source_summary"]
        source_profile = ROOT / row["source_profile"].split("#", 1)[0]
        if (
            not source_summary.is_file()
            or not source_profile.is_file()
            or sha256(source_summary) != row["source_summary_sha256"]
            or sha256(source_profile) != row["source_profile_sha256"]
        ):
            raise ValueError(f"{case_id}: fixed-stiffness source provenance failed")
        reference = next(item for item in main if item["case_id"] == case_id
                         and item["model"] == "scft")
        adaptive = next(item for item in main if item["case_id"] == case_id
                        and item["model"] == "bvk2")
        if (
            row["source_model"] != "bvk2"
            or not math.isclose(
                float(row["source_period_rg"]),
                float(adaptive["period_rg"]),
                rel_tol=0.0,
                abs_tol=1.0e-12,
            )
        ):
            raise ValueError(f"{case_id}: fixed/adaptive source mismatch")
        reference_rows = sorted(scft_native[case_id], key=lambda item: int(item["native_index"]))
        if len(reference_rows) != int(float(reference["grid_count"])):
            raise ValueError(f"{case_id}: SCFT native profile is incomplete")
        ordered_fields = sorted(fields, key=lambda item: int(item["index"]))
        if [int(item["index"]) for item in ordered_fields] != list(range(1, 257)):
            raise ValueError(f"{case_id}: fixed-stiffness profile indices are incomplete")
        reference_profile = [float(item["phi_a"]) for item in reference_rows]
        model_profile = [float(item["phi_a"]) for item in ordered_fields]
        aligned, rms = align(reference_profile, model_profile, math.isclose(f, 0.5))
        period = float(row["period_rg"])
        scft_period = float(reference["period_rg"])
        error = period / scft_period - 1.0
        fixed_summary.append({
            "schema": OUTPUT_SCHEMA, "claim_kind": "stress_free_lamellar_observable",
            "case_id": case_id, "f": f, "chiN": chi, "model": "bvk2_fixed",
            "model_label": LABELS["bvk2_fixed"], "status": "accepted",
            "period_rg": period, "scft_period_rg": scft_period,
            "signed_period_error": error, "abs_period_error": abs(error),
            "profile_rms": rms, "grid_count": 256, "oversample_factor": 2,
            "sensor_filter_ratio": 12.0,
            "projected_force_norm": float(row["projected_force_norm"]),
            "projected_force_maxabs": float(row["projected_force_maxabs"]),
            "stress_norm": float(row["root_cell_stress"]),
        })
        for index, value in enumerate(aligned):
            fixed_profiles.append({
                "schema": OUTPUT_SCHEMA, "case_id": case_id, "f": f,
                "chiN": chi, "model": "bvk2_fixed", "period_rg": period,
                "s": index / SAMPLE_COUNT, "phi_a": float(value),
            })

    aggregate: list[dict[str, object]] = []
    for model in AGGREGATE_MODELS:
        selected = fixed_summary if model == "bvk2_fixed" else [
            row for row in main if row["model"] == model
        ]
        if len(selected) != len(STATE_SET):
            raise ValueError(f"{model}: extension cohort is incomplete")
        aggregate.append({
            "schema": OUTPUT_SCHEMA, "cohort": "integer_chiN_36_50",
            "state_count": len(STATE_SET), "model": model,
            "model_label": LABELS[model],
            "mean_abs_period_error_percent": 100.0 * float(np.mean([
                float(row["abs_period_error"]) for row in selected
            ])),
            "mean_profile_rms": float(np.mean([
                float(row["profile_rms"]) for row in selected
            ])),
        })

    atomic_csv(root / "summary.csv", main)
    atomic_csv(root / "profiles.csv", profiles)
    atomic_csv(root / "native_profiles.csv", native)
    atomic_csv(root / "fixed_stiffness_summary.csv", fixed_summary)
    atomic_csv(root / "fixed_stiffness_profiles.csv", fixed_profiles)
    atomic_csv(root / "aggregate.csv", aggregate)
    report = {
        "schema": OUTPUT_SCHEMA, "status": "passed",
        "cohort": "integer_chiN_36_50", "state_count": len(STATE_SET),
        "main_row_count": len(main), "fixed_row_count": len(fixed_summary),
        "models": [row["model"] for row in aggregate],
    }
    atomic_text(
        root / "validation_report.json",
        json.dumps(report, indent=2, sort_keys=True) + "\n",
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path,
                        default=ROOT / "results/bvk_strong_segregation_extension")
    parser.add_argument(
        "--ok-results-dir",
        type=Path,
        default=ROOT / "results/ok_strong_segregation_extension",
    )
    args = parser.parse_args()
    assemble(args.results_dir.resolve(), args.ok_results_dir.resolve())
    print(f"assembled {OUTPUT_SCHEMA} in {args.results_dir.resolve()}")


if __name__ == "__main__":
    main()
