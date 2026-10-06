#!/usr/bin/env python3
"""Validate refined low-composition SCFT lamellae against canonical references."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import numpy as np


PROJECT = Path(__file__).resolve().parents[1]
CANONICAL = PROJECT / "results" / "liu2019_stress_free_period_map"
AUDIT = PROJECT / "results" / "scft_low_fa_low_chi_resolution_audit"
CASES = ((0.20, 25), (0.20, 26), (0.20, 27),
         (0.25, 19), (0.25, 20), (0.25, 21),
         (0.30, 15), (0.30, 16), (0.30, 17))
PERIOD_TOLERANCE = 1.5e-3
PROFILE_RMS_TOLERANCE = 1.0e-3
START_SPREAD_TOLERANCE = 5.0e-4
SAMPLE_COUNT = 4096


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def truth(value: str) -> bool:
    return value.strip().lower() == "true"


def resample_periodic(values: np.ndarray, count: int = SAMPLE_COUNT) -> np.ndarray:
    source = np.arange(values.size + 1, dtype=float) / values.size
    closed = np.concatenate((values, values[:1]))
    target = np.arange(count, dtype=float) / count
    return np.interp(target, source, closed)


def aligned_rms(reference: np.ndarray, candidate: np.ndarray) -> tuple[float, int]:
    ref = resample_periodic(reference)
    model = resample_periodic(candidate)
    correlation = np.fft.ifft(np.conj(np.fft.fft(ref)) * np.fft.fft(model)).real
    shift = int(np.argmax(correlation))
    aligned = np.roll(model, -shift)
    return float(np.sqrt(np.mean(np.square(aligned - ref)))), shift


def audit_case(directory: Path) -> tuple[dict[str, str], np.ndarray]:
    roots = read_rows(directory / "root.csv")
    profiles = read_rows(directory / "profile.csv")
    if len(roots) != 1:
        raise ValueError(f"expected one root in {directory}")
    profiles.sort(key=lambda row: int(row["index"]))
    if [int(row["index"]) for row in profiles] != list(range(1, len(profiles) + 1)):
        raise ValueError(f"invalid profile inventory in {directory}")
    return roots[0], np.asarray([float(row["phi_a"]) for row in profiles])


def canonical_data() -> tuple[dict[tuple[float, int], dict[str, str]],
                              dict[tuple[float, int], np.ndarray]]:
    summaries: dict[tuple[float, int], dict[str, str]] = {}
    for row in read_rows(CANONICAL / "summary.csv"):
        if row["model"] != "scft":
            continue
        key = (round(float(row["f"]), 8), int(round(float(row["chiN"]))))
        if key in CASES:
            if key in summaries:
                raise ValueError(f"duplicate canonical SCFT summary row for {key}")
            summaries[key] = row

    profile_rows: dict[tuple[float, int], list[dict[str, str]]] = {
        key: [] for key in CASES
    }
    for row in read_rows(CANONICAL / "native_profiles.csv"):
        if row["model"] != "scft":
            continue
        key = (round(float(row["f"]), 8), int(round(float(row["chiN"]))))
        if key in profile_rows:
            profile_rows[key].append(row)

    profiles: dict[tuple[float, int], np.ndarray] = {}
    for key, rows in profile_rows.items():
        rows.sort(key=lambda row: int(row["native_index"]))
        expected = list(range(1, len(rows) + 1))
        actual = [int(row["native_index"]) for row in rows]
        if not rows or actual != expected:
            raise ValueError(f"invalid canonical profile inventory for {key}")
        profiles[key] = np.asarray([float(row["phi_a"]) for row in rows])
    if set(summaries) != set(CASES) or set(profiles) != set(CASES):
        raise ValueError("canonical SCFT inventory does not match the nine-state contract")
    return summaries, profiles


def main() -> None:
    summaries, canonical_profiles = canonical_data()
    comparison: list[dict[str, object]] = []
    failures: list[str] = []

    for f, chi_n in CASES:
        key = (f, chi_n)
        case_dir = AUDIT / "cases" / f"f{f:.2f}_chiN{chi_n}"
        root_rows = read_rows(case_dir / "root.csv")
        attempts = read_rows(case_dir / "attempts.csv")
        profile_rows = read_rows(case_dir / "profile.csv")
        if len(root_rows) != 1 or len(attempts) != 2:
            raise ValueError(f"invalid refined inventory for {key}")
        root = root_rows[0]
        profile_rows.sort(key=lambda row: int(row["index"]))
        indices = [int(row["index"]) for row in profile_rows]
        if indices != list(range(1, len(profile_rows) + 1)):
            raise ValueError(f"invalid refined profile inventory for {key}")
        refined_profile = np.asarray([float(row["phi_a"]) for row in profile_rows])

        canonical = summaries[key]
        canonical_period = float(canonical["period_rg"])
        refined_period = float(root["period_rg"])
        relative_period_change = abs(refined_period / canonical_period - 1.0)
        profile_rms, alignment_shift = aligned_rms(
            canonical_profiles[key], refined_profile
        )
        attempt_periods = [float(row["period_rg"]) for row in attempts]
        independently_computed_spread = max(attempt_periods) / min(attempt_periods) - 1.0

        solver_gate = (
            truth(root["accepted"])
            and truth(root["all_starts_accepted"])
            and root["convergence"] == "Polyorder.Successful()"
            and float(root["residual_norm"]) < 1.0e-6
            and float(root["stress_norm"]) <= 1.0e-5
            and float(root["incompressibility_rms"]) <= 1.0e-6
            and abs(float(root["mean_phi_a"]) - f) <= 1.0e-6
            and int(root["domain_count"]) == 1
            and int(root["dominant_mode"]) == 1
            and len(profile_rows) == int(root["grid_count"])
            and math.isclose(float(root["requested_spacing_rg"]), 0.05,
                             rel_tol=0.0, abs_tol=1.0e-12)
            and math.isclose(float(root["contour_step"]), 0.005,
                             rel_tol=0.0, abs_tol=1.0e-12)
        )
        start_gate = independently_computed_spread <= START_SPREAD_TOLERANCE
        period_gate = relative_period_change <= PERIOD_TOLERANCE
        profile_gate = profile_rms <= PROFILE_RMS_TOLERANCE
        accepted = solver_gate and start_gate and period_gate and profile_gate
        if not accepted:
            failures.append(
                f"f={f:.2f}, chiN={chi_n}: solver={solver_gate}, "
                f"start={start_gate}, period={period_gate}, profile={profile_gate}"
            )

        comparison.append({
            "schema": "scft-low-fa-resolution-audit-v1",
            "f": f,
            "chiN": chi_n,
            "canonical_grid_count": int(canonical["grid_count"]),
            "refined_grid_count": int(root["grid_count"]),
            "canonical_spacing_rg": float(canonical["scft_actual_spacing_rg"]),
            "refined_spacing_rg": float(root["actual_spacing_rg"]),
            "canonical_contour_step": float(canonical["scft_contour_step"]),
            "refined_contour_step": float(root["contour_step"]),
            "canonical_period_rg": canonical_period,
            "refined_period_rg": refined_period,
            "relative_period_change": relative_period_change,
            "profile_rms": profile_rms,
            "alignment_shift_fraction": alignment_shift / SAMPLE_COUNT,
            "two_start_period_spread": independently_computed_spread,
            "refined_residual_norm": float(root["residual_norm"]),
            "refined_stress_norm": float(root["stress_norm"]),
            "solver_gate_pass": solver_gate,
            "start_gate_pass": start_gate,
            "period_gate_pass": period_gate,
            "profile_gate_pass": profile_gate,
            "accepted": accepted,
        })

    AUDIT.mkdir(parents=True, exist_ok=True)
    with (AUDIT / "comparison.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(comparison[0]))
        writer.writeheader()
        writer.writerows(comparison)

    report = {
        "schema": "scft-low-fa-resolution-audit-v1",
        "status": "passed" if not failures else "failed",
        "case_count": len(comparison),
        "accepted_case_count": sum(bool(row["accepted"]) for row in comparison),
        "maximum_relative_period_change": max(
            float(row["relative_period_change"]) for row in comparison
        ),
        "maximum_profile_rms": max(float(row["profile_rms"]) for row in comparison),
        "maximum_two_start_period_spread": max(
            float(row["two_start_period_spread"]) for row in comparison
        ),
        "period_tolerance": PERIOD_TOLERANCE,
        "profile_rms_tolerance": PROFILE_RMS_TOLERANCE,
        "two_start_period_spread_tolerance": START_SPREAD_TOLERANCE,
        "failures": failures,
    }

    decomposition: list[dict[str, object]] = []
    for f, chi_n in CASES[:3]:
        case = f"f{f:.2f}_chiN{chi_n}"
        combined, combined_profile = audit_case(AUDIT / "cases" / case)
        spatial, _ = audit_case(AUDIT / "diagnostics" / "spatial" / case)
        contour, contour_profile = audit_case(
            AUDIT / "diagnostics" / "contour" / case
        )
        finer, finer_profile = audit_case(AUDIT / "diagnostics" / "ds0025" / case)
        full_finer, full_finer_profile = audit_case(
            AUDIT / "diagnostics" / "dx005_ds0025" / case
        )
        canonical_period = float(summaries[(f, chi_n)]["period_rg"])
        contour_period = float(contour["period_rg"])
        finer_period = float(finer["period_rg"])
        contour_refinement_rms, _ = aligned_rms(contour_profile, finer_profile)
        combined_contour_rms, _ = aligned_rms(combined_profile, finer_profile)
        full_refinement_rms, _ = aligned_rms(combined_profile, full_finer_profile)
        decomposition.append({
            "schema": "scft-low-fa-resolution-audit-v1",
            "f": f,
            "chiN": chi_n,
            "canonical_period_rg_dx0p1_ds0p01": canonical_period,
            "spatial_only_period_rg_dx0p05_ds0p01": float(spatial["period_rg"]),
            "contour_only_period_rg_dx0p1_ds0p005": contour_period,
            "combined_period_rg_dx0p05_ds0p005": float(combined["period_rg"]),
            "finer_contour_period_rg_dx0p1_ds0p0025": finer_period,
            "full_finer_period_rg_dx0p05_ds0p0025": float(
                full_finer["period_rg"]
            ),
            "spatial_only_relative_period_change": abs(
                float(spatial["period_rg"]) / canonical_period - 1.0
            ),
            "contour_only_relative_period_change": abs(
                contour_period / canonical_period - 1.0
            ),
            "combined_vs_contour_relative_period_change": abs(
                float(combined["period_rg"]) / contour_period - 1.0
            ),
            "ds0p005_vs_ds0p0025_relative_period_change": abs(
                finer_period / contour_period - 1.0
            ),
            "ds0p005_vs_ds0p0025_profile_rms": contour_refinement_rms,
            "combined_vs_finer_profile_rms": combined_contour_rms,
            "combined_vs_full_finer_relative_period_change": abs(
                float(full_finer["period_rg"]) /
                float(combined["period_rg"]) - 1.0
            ),
            "combined_vs_full_finer_profile_rms": full_refinement_rms,
            "finer_contour_residual_norm": float(finer["residual_norm"]),
            "finer_contour_stress_norm": float(finer["stress_norm"]),
            "finer_contour_accepted": truth(finer["accepted"]),
            "full_finer_accepted": truth(full_finer["accepted"]),
        })
    with (AUDIT / "resolution_decomposition.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(decomposition[0]))
        writer.writeheader()
        writer.writerows(decomposition)

    report["fA_0p20_diagnosis"] = {
        "dominant_error_source": "contour_discretization",
        "maximum_spatial_only_relative_period_change": max(
            float(row["spatial_only_relative_period_change"])
            for row in decomposition
        ),
        "minimum_contour_only_relative_period_change": min(
            float(row["contour_only_relative_period_change"])
            for row in decomposition
        ),
        "maximum_ds0p005_vs_ds0p0025_relative_period_change": max(
            float(row["ds0p005_vs_ds0p0025_relative_period_change"])
            for row in decomposition
        ),
        "maximum_ds0p005_vs_ds0p0025_profile_rms": max(
            float(row["ds0p005_vs_ds0p0025_profile_rms"])
            for row in decomposition
        ),
        "maximum_combined_vs_full_finer_relative_period_change": max(
            float(row["combined_vs_full_finer_relative_period_change"])
            for row in decomposition
        ),
        "maximum_combined_vs_full_finer_profile_rms": max(
            float(row["combined_vs_full_finer_profile_rms"])
            for row in decomposition
        ),
        "recommended_reference_contour_step": 0.005,
    }
    (AUDIT / "validation_report.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    if failures:
        raise SystemExit("SCFT refinement audit failed:\n" + "\n".join(failures))
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
