#!/usr/bin/env python3
"""Assemble the accepted seven-state fixed-stiffness ABA control."""

from __future__ import annotations

import csv
import hashlib
import math
import os
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTDIR = ROOT / "results" / "bvk2_aba_fixed_stiffness_validation"
SOURCE_DIR = ROOT / "results" / "bvk2_aba_comprehensive_validation"
SCHEMA = "bvk2-aba-fixed-stiffness-validation-v1"
CASE_IDS = (
    "f0.50_chiN20",
    "f0.50_chiN22",
    "f0.50_chiN25",
    "f0.50_chiN30",
    "f0.50_chiN35",
    "f0.50_chiN40",
    "f0.50_chiN45",
)
PROFILE_SAMPLES = 256


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def truth(value: str) -> bool:
    return value.lower() == "true"


def validate_case(
    case_id: str, root: list[dict[str, str]], profiles: list[dict[str, str]]
) -> None:
    if len(root) != 1:
        raise ValueError(f"{case_id} requires exactly one root row")
    row = root[0]
    if row.get("schema") != SCHEMA or row.get("case_id") != case_id:
        raise ValueError(f"{case_id} root metadata is inconsistent")
    if row.get("model") != "bvk2_fixed" or truth(row.get("adaptive", "true")):
        raise ValueError(f"{case_id} is not the fixed-stiffness control")
    if float(row.get("c2", "nan")) != 0.16 or int(row.get("nx", "0")) != 128:
        raise ValueError(f"{case_id} numerical settings are inconsistent")
    gates = (
        "field_gate_pass",
        "composition_gate_pass",
        "morphology_gate_pass",
        "cell_gate_pass",
        "local_minimum_check_pass",
        "accepted",
    )
    if any(not truth(row.get(gate, "false")) for gate in gates):
        raise ValueError(f"{case_id} did not pass every acceptance gate")
    if row.get("status") != "accepted":
        raise ValueError(f"{case_id} is not accepted")
    numeric = (
        "fA",
        "chiN",
        "period_rg",
        "scft_period_rg",
        "period_relative_error",
        "profile_rms",
        "profile_correlation",
        "energy_density",
        "projected_force_norm",
        "projected_force_maxabs",
    )
    if not all(math.isfinite(float(row[field])) for field in numeric):
        raise ValueError(f"{case_id} root contains nonfinite data")
    source_summary = ROOT / row["source_summary"]
    source_profiles = ROOT / row["source_profiles"]
    if sha256(source_summary) != row["source_summary_sha256"]:
        raise ValueError(f"{case_id} summary provenance mismatch")
    if sha256(source_profiles) != row["source_profiles_sha256"]:
        raise ValueError(f"{case_id} profile provenance mismatch")

    if len(profiles) != PROFILE_SAMPLES:
        raise ValueError(f"{case_id} requires {PROFILE_SAMPLES} profile samples")
    indices: set[int] = set()
    for profile in profiles:
        if (
            profile.get("schema") != SCHEMA
            or profile.get("case_id") != case_id
            or profile.get("model") != "bvk2_fixed"
            or truth(profile.get("adaptive", "true"))
        ):
            raise ValueError(f"{case_id} profile metadata is inconsistent")
        index = int(profile["index"])
        indices.add(index)
        values = (float(profile["s"]), float(profile["phi_a"]))
        if not all(math.isfinite(value) for value in values):
            raise ValueError(f"{case_id} profile contains nonfinite data")
    if indices != set(range(1, PROFILE_SAMPLES + 1)):
        raise ValueError(f"{case_id} profile indices are incomplete")


def atomic_csv(path: Path, rows: list[dict[str, str]]) -> None:
    if not rows:
        raise ValueError(f"cannot write empty table {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w", newline="", encoding="utf-8", dir=path.parent, delete=False
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
        temporary = Path(handle.name)
    os.replace(temporary, path)


def assemble(outdir: Path = OUTDIR) -> tuple[Path, Path]:
    roots: list[dict[str, str]] = []
    profiles: list[dict[str, str]] = []
    for case_id in CASE_IDS:
        case_dir = outdir / "cases" / case_id
        root_rows = read_rows(case_dir / "root.csv")
        profile_rows = read_rows(case_dir / "profile.csv")
        validate_case(case_id, root_rows, profile_rows)
        roots.extend(root_rows)
        profiles.extend(profile_rows)
    roots.sort(key=lambda row: float(row["chiN"]))
    profiles.sort(key=lambda row: (float(row["chiN"]), int(row["index"])))
    summary_path = outdir / "summary.csv"
    profile_path = outdir / "profiles.csv"
    atomic_csv(summary_path, roots)
    atomic_csv(profile_path, profiles)
    return summary_path, profile_path


if __name__ == "__main__":
    summary, profiles = assemble()
    print(summary.relative_to(ROOT))
    print(profiles.relative_to(ROOT))
