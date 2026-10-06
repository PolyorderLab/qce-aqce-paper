#!/usr/bin/env python3
"""Validate the OPF grid-doubling check for promoted low-fA SCFT states."""

from __future__ import annotations

import csv
import json
import math
import os
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CANONICAL = ROOT / "results" / "liu2019_stress_free_period_map"
AUDIT = ROOT / "results" / "scft_low_fa_low_chi_promotion" / "opf_nx256"
OUTPUT = ROOT / "results" / "scft_low_fa_low_chi_promotion"
CASES = (25, 26, 27)
MODEL = "liu2019_opf"
INTERPOLATION_COUNT = 4096
PERIOD_TOLERANCE = 1.0e-6
PROFILE_METRIC_TOLERANCE = 1.0e-4
NATIVE_PROFILE_TOLERANCE = 1.0e-4


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def truth(value: str) -> bool:
    return value.strip().lower() == "true"


def model_row(path: Path) -> dict[str, str]:
    selected = [row for row in rows(path) if row["model"] == MODEL]
    if len(selected) != 1:
        raise ValueError(f"{path}: expected exactly one {MODEL} row")
    return selected[0]


def native_profile(path: Path) -> tuple[list[float], list[float], list[float]]:
    selected = [row for row in rows(path) if row["model"] == MODEL]
    selected.sort(key=lambda row: int(row["native_index"]))
    if len(selected) < 64:
        raise ValueError(f"{path}: incomplete {MODEL} native profile")
    expected = list(range(1, len(selected) + 1))
    if [int(row["native_index"]) for row in selected] != expected:
        raise ValueError(f"{path}: noncontiguous native profile inventory")
    return (
        [float(row["s"]) for row in selected],
        [float(row["phi_a"]) for row in selected],
        [float(row["phi_b"]) for row in selected],
    )


def periodic_interpolate(
    coordinates: list[float], values: list[float], count: int
) -> list[float]:
    if len(coordinates) != len(values) or not coordinates:
        raise ValueError("invalid periodic profile")
    result: list[float] = []
    right = 1
    for index in range(count):
        target = index / count
        while right < len(coordinates) and coordinates[right] <= target:
            right += 1
        left = right - 1
        if right == len(coordinates):
            x0, y0 = coordinates[left], values[left]
            x1, y1 = 1.0, values[0]
        else:
            x0, y0 = coordinates[left], values[left]
            x1, y1 = coordinates[right], values[right]
        fraction = 0.0 if x1 == x0 else (target - x0) / (x1 - x0)
        result.append((1.0 - fraction) * y0 + fraction * y1)
    return result


def rms(left: list[float], right: list[float], shift: int = 0) -> float:
    count = len(left)
    return math.sqrt(
        sum((left[index] - right[(index - shift) % count]) ** 2 for index in range(count))
        / count
    )


def aligned_native_rms(
    coarse: tuple[list[float], list[float], list[float]],
    fine: tuple[list[float], list[float], list[float]],
) -> tuple[float, float, float]:
    coarse_a = periodic_interpolate(coarse[0], coarse[1], INTERPOLATION_COUNT)
    coarse_b = periodic_interpolate(coarse[0], coarse[2], INTERPOLATION_COUNT)
    fine_a = periodic_interpolate(fine[0], fine[1], INTERPOLATION_COUNT)
    fine_b = periodic_interpolate(fine[0], fine[2], INTERPOLATION_COUNT)
    best = min(
        (
            0.5 * (rms(coarse_a, fine_a, shift) + rms(coarse_b, fine_b, shift)),
            shift,
        )
        for shift in range(INTERPOLATION_COUNT)
    )
    shift = best[1]
    return rms(coarse_a, fine_a, shift), rms(coarse_b, fine_b, shift), shift / INTERPOLATION_COUNT


def atomic_csv(path: Path, data: list[dict[str, object]]) -> None:
    fields = list(data[0])
    with tempfile.NamedTemporaryFile(
        "w", newline="", encoding="utf-8", dir=path.parent, delete=False
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(data)
        temporary = Path(handle.name)
    os.replace(temporary, path)


def atomic_json(path: Path, payload: dict[str, object]) -> None:
    with tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", dir=path.parent, delete=False
    ) as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
        temporary = Path(handle.name)
    os.replace(temporary, path)


def main() -> None:
    comparison: list[dict[str, object]] = []
    failures: list[str] = []
    required_gates = (
        "field_gate_pass",
        "cell_gate_pass",
        "resolution_gate_pass",
        "composition_gate_pass",
        "morphology_gate_pass",
        "period_local_minimum_check_pass",
    )
    for chi_n in CASES:
        case_id = f"f0.2_chiN{chi_n}"
        coarse_dir = CANONICAL / "cases" / case_id
        fine_dir = AUDIT / "cases" / case_id
        coarse = model_row(coarse_dir / "summary.csv")
        fine = model_row(fine_dir / "summary.csv")
        gates_pass = all(
            row["status"] == "accepted" and all(truth(row[gate]) for gate in required_gates)
            for row in (coarse, fine)
        )
        period_change = abs(float(fine["period_rg"]) / float(coarse["period_rg"]) - 1.0)
        metric_change = abs(float(fine["profile_rms"]) - float(coarse["profile_rms"]))
        phi_a_rms, phi_b_rms, shift = aligned_native_rms(
            native_profile(coarse_dir / "native_profiles.csv"),
            native_profile(fine_dir / "native_profiles.csv"),
        )
        accepted = (
            gates_pass
            and int(coarse["grid_count"]) == 128
            and int(fine["grid_count"]) == 256
            and period_change < PERIOD_TOLERANCE
            and metric_change < PROFILE_METRIC_TOLERANCE
            and max(phi_a_rms, phi_b_rms) < NATIVE_PROFILE_TOLERANCE
        )
        if not accepted:
            failures.append(case_id)
        comparison.append(
            {
                "case_id": case_id,
                "f": 0.20,
                "chiN": chi_n,
                "coarse_grid_count": int(coarse["grid_count"]),
                "fine_grid_count": int(fine["grid_count"]),
                "coarse_period_rg": float(coarse["period_rg"]),
                "fine_period_rg": float(fine["period_rg"]),
                "relative_period_change": period_change,
                "coarse_profile_rms": float(coarse["profile_rms"]),
                "fine_profile_rms": float(fine["profile_rms"]),
                "profile_rms_metric_change": metric_change,
                "aligned_phi_a_rms": phi_a_rms,
                "aligned_phi_b_rms": phi_b_rms,
                "alignment_shift_fraction": shift,
                "all_solution_gates_pass": str(gates_pass).lower(),
                "accepted": str(accepted).lower(),
            }
        )

    OUTPUT.mkdir(parents=True, exist_ok=True)
    atomic_csv(OUTPUT / "opf_nx256_comparison.csv", comparison)
    atomic_json(
        OUTPUT / "opf_nx256_validation.json",
        {
            "schema": "promoted-opf-resolution-audit-v1",
            "status": "passed" if not failures else "failed",
            "case_count": len(comparison),
            "accepted_case_count": len(comparison) - len(failures),
            "failures": failures,
            "tolerances": {
                "relative_period_change": PERIOD_TOLERANCE,
                "profile_rms_metric_change": PROFILE_METRIC_TOLERANCE,
                "aligned_native_profile_rms": NATIVE_PROFILE_TOLERANCE,
            },
            "max_relative_period_change": max(
                float(row["relative_period_change"]) for row in comparison
            ),
            "max_profile_rms_metric_change": max(
                float(row["profile_rms_metric_change"]) for row in comparison
            ),
            "max_aligned_native_profile_rms": max(
                max(float(row["aligned_phi_a_rms"]), float(row["aligned_phi_b_rms"]))
                for row in comparison
            ),
        },
    )
    if failures:
        raise SystemExit(f"OPF resolution audit failed: {', '.join(failures)}")


if __name__ == "__main__":
    main()
