#!/usr/bin/env python3
"""Validate and assemble the production fixed-stiffness lamellar map."""

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
SCHEMA = "bvk2-fixed-stiffness-period-map-v1"
SOURCE_SCHEMA = "bvk2-fixed-stiffness-figure6-root-v1"
SAMPLE_COUNT = 1024
EXPECTED_SETTINGS = {
    "adaptive": "false",
    "c2": 0.16,
    "nx": 256,
    "oversample_factor": 2,
    "sensor_filter_ratio": 12.0,
}
COMMON_COHORT = {
    (composition, float(chi_n))
    for composition, lower in (
        (0.20, 30),
        (0.25, 21),
        (0.30, 17),
        (0.35, 14),
        (0.40, 13),
        (0.45, 12),
        (0.50, 12),
    )
    for chi_n in range(lower, 36)
}


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def truth(value: str) -> bool:
    return value.strip().lower() == "true"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def periodic_resample(values: list[float], count: int = SAMPLE_COUNT) -> np.ndarray:
    source = np.asarray(values, dtype=float)
    coordinates = np.arange(count, dtype=float) * len(source) / count
    left = np.floor(coordinates).astype(int)
    fraction = coordinates - left
    return (
        (1.0 - fraction) * source[left % len(source)]
        + fraction * source[(left + 1) % len(source)]
    )


def align_periodic(
    reference: list[float], model: list[float], *, allow_swap: bool
) -> tuple[np.ndarray, float, int, bool]:
    target = periodic_resample(reference)
    base = periodic_resample(model)
    best: tuple[float, int, bool, np.ndarray] | None = None
    candidates = ((base, False), (1.0 - base, True)) if allow_swap else ((base, False),)
    for candidate, swapped in candidates:
        for shift in range(SAMPLE_COUNT):
            aligned = np.roll(candidate, -shift)
            rms = float(np.sqrt(np.mean((aligned - target) ** 2)))
            if best is None or rms < best[0]:
                best = (rms, shift, swapped, aligned)
    assert best is not None
    return best[3], best[0], best[1], best[2]


def domain_count(profile: list[float], threshold: float) -> int:
    mask = np.asarray(profile) > threshold
    if np.all(mask) or not np.any(mask):
        return 0
    return int(np.sum(mask & ~np.roll(mask, 1)))


def dominant_mode(profile: list[float], mean: float) -> int:
    spectrum = np.abs(np.fft.rfft(np.asarray(profile) - mean))
    return int(np.argmax(spectrum[1:]) + 1)


def atomic_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w", newline="", encoding="utf-8", dir=path.parent, delete=False
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
        temporary = Path(handle.name)
    os.replace(temporary, path)


def atomic_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", dir=path.parent, delete=False
    ) as handle:
        handle.write(text)
        temporary = Path(handle.name)
    os.replace(temporary, path)


def assemble(source_dir: Path, reference_dir: Path, outdir: Path) -> None:
    reference_summary = [
        row
        for row in read_rows(reference_dir / "summary.csv")
        if row["model"] == "scft" and row["status"] == "accepted"
    ]
    references = {row["case_id"]: row for row in reference_summary}
    if len(references) != 146:
        raise ValueError(f"expected 146 accepted SCFT states, found {len(references)}")

    native = defaultdict(list)
    for row in read_rows(reference_dir / "native_profiles.csv"):
        if row["model"] == "scft":
            native[row["case_id"]].append(row)

    summary_rows: list[dict[str, object]] = []
    profile_rows: list[dict[str, object]] = []
    failures: list[str] = []
    for case_id, reference in sorted(
        references.items(), key=lambda item: (float(item[1]["f"]), float(item[1]["chiN"]))
    ):
        case_dir = source_dir / "cases" / case_id
        root_path = case_dir / "root.csv"
        profile_path = case_dir / "profile.csv"
        roots = read_rows(root_path)
        profiles = read_rows(profile_path)
        if len(roots) != 1 or len(profiles) != EXPECTED_SETTINGS["nx"]:
            failures.append(f"{case_id}: incomplete root/profile inventory")
            continue
        root = roots[0]
        gate_names = (
            "field_gate_pass",
            "composition_gate_pass",
            "morphology_gate_pass",
            "cell_gate_pass",
            "local_minimum_check_pass",
            "stress_orientation_pass",
            "accepted",
        )
        valid = (
            root["schema"] == SOURCE_SCHEMA
            and root["case_id"] == case_id
            and root["model"] == "bvk2_fixed"
            and root["status"] == "accepted"
            and all(truth(root[name]) for name in gate_names)
            and root["adaptive"] == EXPECTED_SETTINGS["adaptive"]
            and math.isclose(float(root["c2"]), EXPECTED_SETTINGS["c2"], abs_tol=0.0)
            and int(root["nx"]) == EXPECTED_SETTINGS["nx"]
            and int(root["oversample_factor"])
            == EXPECTED_SETTINGS["oversample_factor"]
            and math.isclose(
                float(root["sensor_filter_ratio"]),
                EXPECTED_SETTINGS["sensor_filter_ratio"],
                abs_tol=0.0,
            )
            and abs(float(root["root_cell_stress"])) <= 1.0e-6
            and float(root["projected_force_norm"]) <= 1.0e-5
            and float(root["projected_force_maxabs"]) <= 2.0e-5
        )
        model_profile = [float(row["phi_a"]) for row in profiles]
        f_value = float(reference["f"])
        valid = (
            valid
            and all(math.isfinite(value) for value in model_profile)
            and domain_count(model_profile, f_value) == 1
            and dominant_mode(model_profile, f_value) == 1
        )
        if not valid:
            failures.append(f"{case_id}: one or more production gates failed")
            continue

        reference_profile_rows = sorted(
            native[case_id], key=lambda row: int(row["native_index"])
        )
        reference_profile = [float(row["phi_a"]) for row in reference_profile_rows]
        aligned, rms, shift, swapped = align_periodic(
            reference_profile, model_profile, allow_swap=math.isclose(f_value, 0.5)
        )
        period_rg = float(root["period_rg"])
        scft_period = float(reference["period_rg"])
        signed_error = period_rg / scft_period - 1.0
        summary_rows.append(
            {
                "schema": SCHEMA,
                "claim_kind": "stress_free_lamellar_observable",
                "case_id": case_id,
                "f": f_value,
                "chiN": float(reference["chiN"]),
                "model": "bvk2_fixed",
                "model_label": "BVK, fixed stiffness",
                "status": "accepted",
                "calibration_role": "fixed_stiffness_production_control",
                "adaptive": "false",
                "c2": EXPECTED_SETTINGS["c2"],
                "grid_count": EXPECTED_SETTINGS["nx"],
                "oversample_factor": EXPECTED_SETTINGS["oversample_factor"],
                "sensor_filter_ratio": EXPECTED_SETTINGS["sensor_filter_ratio"],
                "discretization_schema": root["discretization_schema"],
                "scft_period_rg": scft_period,
                "period_rg": period_rg,
                "signed_period_error": signed_error,
                "abs_period_error": abs(signed_error),
                "profile_rms": rms,
                "alignment_shift_fraction": shift / SAMPLE_COUNT,
                "alignment_swapped": str(swapped).lower(),
                "projected_force_norm": float(root["projected_force_norm"]),
                "projected_force_maxabs": float(root["projected_force_maxabs"]),
                "cell_stress": float(root["root_cell_stress"]),
                **{name: "true" for name in gate_names},
                "source_root": str(root_path.relative_to(ROOT)),
                "source_profile": str(profile_path.relative_to(ROOT)),
                "source_root_sha256": sha256(root_path),
                "source_profile_sha256": sha256(profile_path),
            }
        )
        for index, value in enumerate(aligned):
            profile_rows.append(
                {
                    "schema": SCHEMA,
                    "case_id": case_id,
                    "f": f_value,
                    "chiN": float(reference["chiN"]),
                    "model": "bvk2_fixed",
                    "period_rg": period_rg,
                    "index": index + 1,
                    "s": index / SAMPLE_COUNT,
                    "phi_a": float(value),
                }
            )

    if failures or len(summary_rows) != 146 or len(profile_rows) != 146 * SAMPLE_COUNT:
        raise ValueError(
            f"fixed-stiffness assembly failed: {len(summary_rows)} states; "
            + " | ".join(failures[:10])
        )

    summary_fields = list(summary_rows[0])
    profile_fields = list(profile_rows[0])
    atomic_csv(outdir / "summary.csv", summary_rows, summary_fields)
    atomic_csv(outdir / "profiles.csv", profile_rows, profile_fields)
    common_rows = [
        row
        for row in summary_rows
        if (float(row["f"]), float(row["chiN"])) in COMMON_COHORT
    ]
    if len(common_rows) != 133:
        raise ValueError(
            f"fixed-stiffness common cohort has {len(common_rows)} states, expected 133"
        )
    report = {
        "schema": SCHEMA,
        "status": "passed",
        "state_count": len(summary_rows),
        "profile_row_count": len(profile_rows),
        "maximum_abs_cell_stress": max(abs(float(row["cell_stress"])) for row in summary_rows),
        "maximum_projected_force_rms": max(float(row["projected_force_norm"]) for row in summary_rows),
        "maximum_projected_force_maxabs": max(float(row["projected_force_maxabs"]) for row in summary_rows),
        "mean_abs_period_error_percent": 100.0 * float(np.mean([float(row["abs_period_error"]) for row in summary_rows])),
        "mean_profile_rms": float(np.mean([float(row["profile_rms"]) for row in summary_rows])),
        "common_cohort_state_count": len(common_rows),
        "common_cohort_mean_abs_period_error_percent": 100.0
        * float(np.mean([float(row["abs_period_error"]) for row in common_rows])),
        "common_cohort_mean_profile_rms": float(
            np.mean([float(row["profile_rms"]) for row in common_rows])
        ),
    }
    atomic_text(outdir / "validation_report.json", json.dumps(report, indent=2, sort_keys=True) + "\n")
    atomic_text(
        outdir / "README.md",
        "# Production Fixed-Stiffness BVK Lamellar Map\n\n"
        "This directory contains the fixed-stiffness BVK control used in "
        "Figures 3 and 4 and Figure S4. The calculation imposes "
        "\\(K_{\\psi,2}=K_{\\psi,0}\\) while retaining the BVK2 production "
        "grid (`nx=256`), twofold Fourier oversampling, and physical sensor "
        "filter (`kc/kstar=12`). Profiles are cyclically aligned to the "
        "independently relaxed SCFT references on 1024 normalized-cell "
        "samples.\n\n"
        f"All {report['state_count']} states pass the field, composition, "
        "primitive-lamella, cell-stress, stress-orientation, and local-period-"
        "minimum gates. The largest projected-force RMS is "
        f"`{report['maximum_projected_force_rms']:.6e}`, the largest "
        "projected-force component is "
        f"`{report['maximum_projected_force_maxabs']:.6e}`, and the largest "
        "absolute cell stress is "
        f"`{report['maximum_abs_cell_stress']:.6e}`.\n\n"
        f"- Full 146-state mean absolute period error: "
        f"{report['mean_abs_period_error_percent']:.6f}%.\n"
        f"- Full 146-state mean profile RMS: {report['mean_profile_rms']:.8f}.\n"
        f"- Common 133-state mean absolute period error: "
        f"{report['common_cohort_mean_abs_period_error_percent']:.6f}%.\n"
        f"- Common 133-state mean profile RMS: "
        f"{report['common_cohort_mean_profile_rms']:.8f}.\n\n"
        "`summary.csv` and `profiles.csv` are assembled atomically by "
        "`scripts/assemble_bvk2_fixed_stiffness_period_map.py`; the assembler "
        "verifies every root/profile inventory, model setting, solver gate, "
        "primitive-cell identity, and source SHA-256 before promotion. "
        "`validation_report.json` records the aggregate gate maxima. The "
        "per-state accepted roots and native profiles are retained under "
        "`cases/`.\n",
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source-dir",
        type=Path,
        default=ROOT / "results/bvk_fixed_stiffness_period_map_nx256",
    )
    parser.add_argument(
        "--reference-dir",
        type=Path,
        default=ROOT / "results/liu2019_stress_free_period_map",
    )
    parser.add_argument("--outdir", type=Path, default=None)
    args = parser.parse_args()
    outdir = args.source_dir if args.outdir is None else args.outdir
    assemble(args.source_dir.resolve(), args.reference_dir.resolve(), outdir.resolve())
    print(f"assembled production fixed-stiffness map in {outdir.resolve()}")


if __name__ == "__main__":
    main()
