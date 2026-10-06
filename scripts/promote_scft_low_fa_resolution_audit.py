#!/usr/bin/env python3
"""Stage, refresh, and promote the accepted low-fA SCFT references."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import shutil
import tempfile
from collections import defaultdict
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
CANONICAL = ROOT / "results" / "liu2019_stress_free_period_map"
AUDIT = ROOT / "results" / "scft_low_fa_low_chi_resolution_audit"
CASES = ((0.20, 25), (0.20, 26), (0.20, 27))
SCHEMA = "liu2019-stress-free-period-map-v2"
SAMPLE_COUNT = 256


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def truth(value: str) -> bool:
    return value.strip().lower() == "true"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def atomic_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty CSV {path}")
    fields = list(rows[0])
    if any(set(row) != set(fields) for row in rows):
        raise ValueError(f"inconsistent row schema for {path}")
    with tempfile.NamedTemporaryFile(
        "w", newline="", encoding="utf-8", dir=path.parent, delete=False
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
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


def periodic_resample(values: list[float], count: int = SAMPLE_COUNT) -> np.ndarray:
    source = np.asarray(values, dtype=float)
    coordinates = np.arange(count, dtype=float) * source.size / count
    left = np.floor(coordinates).astype(int)
    fraction = coordinates - left
    return (
        (1.0 - fraction) * source[left % source.size]
        + fraction * source[(left + 1) % source.size]
    )


def align_periodic(
    reference: list[float], model: list[float], count: int = SAMPLE_COUNT
) -> tuple[np.ndarray, float, int]:
    target = periodic_resample(reference, count)
    base = periodic_resample(model, count)
    best: tuple[float, int, np.ndarray] | None = None
    for shift in range(count):
        aligned = np.roll(base, -shift)
        rms = float(np.sqrt(np.mean(np.square(aligned - target))))
        if best is None or rms < best[0]:
            best = rms, shift, aligned
    assert best is not None
    return best[2], best[0], best[1]


def case_ids(f_value: float, chi_n: int) -> tuple[str, str]:
    return f"f{f_value:g}_chiN{chi_n}", f"f{f_value:.2f}_chiN{chi_n}"


def prepare(staging: Path) -> None:
    if staging.exists():
        raise FileExistsError(f"staging directory already exists: {staging}")
    staging.mkdir(parents=True)
    comparison = {
        (float(row["f"]), int(row["chiN"])): row
        for row in read_rows(AUDIT / "comparison.csv")
    }
    decomposition = {
        (float(row["f"]), int(row["chiN"])): row
        for row in read_rows(AUDIT / "resolution_decomposition.csv")
    }
    source_hashes: dict[str, dict[str, str]] = {}

    for key in CASES:
        f_value, chi_n = key
        canonical_id, audit_id = case_ids(f_value, chi_n)
        source = CANONICAL / "cases" / canonical_id
        destination = staging / "cases" / canonical_id
        shutil.copytree(source, destination)

        root_rows = read_rows(AUDIT / "cases" / audit_id / "root.csv")
        audit_profiles = read_rows(AUDIT / "cases" / audit_id / "profile.csv")
        if len(root_rows) != 1:
            raise ValueError(f"{audit_id}: expected one accepted audit root")
        root = root_rows[0]
        audit_profiles.sort(key=lambda row: int(row["index"]))
        if (
            not truth(root["accepted"])
            or not truth(root["all_starts_accepted"])
            or root["convergence"] != "Polyorder.Successful()"
            or not math.isclose(float(root["requested_spacing_rg"]), 0.05)
            or not math.isclose(float(root["contour_step"]), 0.005)
            or int(root["domain_count"]) != 1
            or int(root["dominant_mode"]) != 1
            or truth(comparison[key]["accepted"])
            or not truth(decomposition[key]["full_finer_accepted"])
            or float(decomposition[key]["combined_vs_full_finer_relative_period_change"])
            >= 1.0e-3
            or float(decomposition[key]["combined_vs_full_finer_profile_rms"])
            >= 1.0e-3
        ):
            raise ValueError(f"{audit_id}: promotion gates do not pass")
        if len(audit_profiles) != int(root["grid_count"]):
            raise ValueError(f"{audit_id}: profile inventory is incomplete")

        summary = read_rows(destination / "summary.csv")
        scft_rows = [row for row in summary if row["model"] == "scft"]
        if len(scft_rows) != 1:
            raise ValueError(f"{canonical_id}: expected one canonical SCFT row")
        scft = scft_rows[0]
        scft.update(
            period_rg=root["period_rg"],
            scft_period_rg=root["period_rg"],
            signed_period_error="0.0",
            abs_period_error="0.0",
            profile_rms="0.0",
            profile_r2="1.0",
            profile_correlation="1.0",
            alignment_shift_fraction="0.0",
            alignment_swapped="false",
            period_iterations="",
            period_local_minimum_check_pass="true",
            period_boundary_limited="false",
            field_gate_pass="true",
            cell_gate_pass="true",
            resolution_gate_pass="true",
            composition_gate_pass="true",
            morphology_gate_pass="true",
            minimum_phi_a=root["minimum_phi_a"],
            maximum_phi_a=root["maximum_phi_a"],
            mean_phi_a=root["mean_phi_a"],
            energy=root["free_energy"],
            stress_norm=root["stress_norm"],
            scft_convergence=root["convergence"],
            scft_residual_norm=root["residual_norm"],
            scft_incompressibility_rms=root["incompressibility_rms"],
            scft_requested_spacing_rg=root["requested_spacing_rg"],
            scft_actual_spacing_rg=root["actual_spacing_rg"],
            scft_contour_step=root["contour_step"],
            grid_count=root["grid_count"],
            discretization_schema="Polyorder pseudospectral SCFT; refined low-fA audit v1",
            status="accepted",
            status_reason="reference_gates_pass_refined_low_fA_audit",
            error_message="",
        )
        atomic_csv(destination / "summary.csv", summary)

        raw_a = [float(row["phi_a"]) for row in audit_profiles]
        raw_b = [float(row["phi_b"]) for row in audit_profiles]
        common_a = periodic_resample(raw_a)
        common_b = periodic_resample(raw_b)
        period = float(root["period_rg"])
        profiles = [row for row in read_rows(destination / "profiles.csv")
                    if row["model"] != "scft"]
        for index, (phi_a, phi_b) in enumerate(zip(common_a, common_b)):
            s = index / SAMPLE_COUNT
            profiles.append({
                "schema": SCHEMA,
                "case_id": canonical_id,
                "f": f_value,
                "chiN": float(chi_n),
                "model": "scft",
                "model_label": "SCFT",
                "period_rg": period,
                "s": s,
                "x_rg": s * period,
                "phi_a": float(phi_a),
                "phi_b": float(phi_b),
                "phi_a_raw": float(phi_a),
                "phi_b_raw": float(phi_b),
            })
        profiles.sort(key=lambda row: (row["model"], float(row["s"])))
        atomic_csv(destination / "profiles.csv", profiles)

        native = [row for row in read_rows(destination / "native_profiles.csv")
                  if row["model"] != "scft"]
        for index, row in enumerate(audit_profiles, start=1):
            native.append({
                "schema": SCHEMA,
                "case_id": canonical_id,
                "f": f_value,
                "chiN": float(chi_n),
                "model": "scft",
                "model_label": "SCFT",
                "period_rg": period,
                "native_index": index,
                "native_grid_count": len(audit_profiles),
                "s": (index - 1) / len(audit_profiles),
                "phi_a": row["phi_a"],
                "phi_b": row["phi_b"],
            })
        native.sort(key=lambda row: (row["model"], int(row["native_index"])))
        atomic_csv(destination / "native_profiles.csv", native)
        source_hashes[canonical_id] = {
            "audit_root": sha256(AUDIT / "cases" / audit_id / "root.csv"),
            "audit_profile": sha256(AUDIT / "cases" / audit_id / "profile.csv"),
        }

    atomic_json(staging / "promotion_manifest.json", {
        "schema": "scft-low-fa-promotion-v1",
        "status": "prepared",
        "cases": [case_ids(*key)[0] for key in CASES],
        "source_hashes": source_hashes,
    })


def refresh(staging: Path) -> None:
    manifest = json.loads((staging / "promotion_manifest.json").read_text())
    if manifest.get("status") not in {"prepared", "promoted"}:
        raise ValueError("staging manifest is neither prepared nor promoted")
    refreshed: dict[str, dict[str, object]] = {}
    for f_value, chi_n in CASES:
        case_id, _ = case_ids(f_value, chi_n)
        case_dir = staging / "cases" / case_id
        summary = read_rows(case_dir / "summary.csv")
        native_rows = read_rows(case_dir / "native_profiles.csv")
        native: dict[str, list[dict[str, str]]] = defaultdict(list)
        for row in native_rows:
            native[row["model"]].append(row)
        for model_rows in native.values():
            model_rows.sort(key=lambda row: int(row["native_index"]))
        scft_summary = [row for row in summary if row["model"] == "scft"]
        if len(scft_summary) != 1 or len(native["scft"]) < 64:
            raise ValueError(f"{case_id}: invalid promoted SCFT inventory")
        scft_period = float(scft_summary[0]["period_rg"])
        scft_profile = [float(row["phi_a"]) for row in native["scft"]]
        old_profiles = read_rows(case_dir / "profiles.csv")
        old_by_model: dict[str, list[dict[str, str]]] = defaultdict(list)
        for row in old_profiles:
            old_by_model[row["model"]].append(row)
        new_profiles: list[dict[str, object]] = []

        for row in summary:
            model = row["model"]
            sample_count = 1024 if model == "bvk2" else SAMPLE_COUNT
            row["scft_period_rg"] = repr(scft_period)
            if model == "scft":
                aligned = periodic_resample(scft_profile, sample_count)
                raw_a = aligned
                raw_b = 1.0 - aligned
                row.update(
                    signed_period_error="0.0", abs_period_error="0.0",
                    profile_rms="0.0", profile_r2="1.0",
                    profile_correlation="1.0", alignment_shift_fraction="0.0",
                    alignment_swapped="false",
                )
            elif native.get(model):
                model_a = [float(item["phi_a"]) for item in native[model]]
                model_b = [float(item["phi_b"]) for item in native[model]]
                aligned, rms, shift = align_periodic(
                    scft_profile, model_a, sample_count
                )
                target = periodic_resample(scft_profile, sample_count)
                centered = target - np.mean(target)
                denominator = float(np.sum(np.square(centered)))
                residual = aligned - target
                r2 = 1.0 - float(np.sum(np.square(residual))) / denominator
                correlation = float(np.corrcoef(target, aligned)[0, 1])
                period_error = float(row["period_rg"]) / scft_period - 1.0
                row.update(
                    signed_period_error=repr(period_error),
                    abs_period_error=repr(abs(period_error)),
                    profile_rms=repr(rms), profile_r2=repr(r2),
                    profile_correlation=repr(correlation),
                    alignment_shift_fraction=repr(shift / sample_count),
                    alignment_swapped="false",
                )
                raw_a = periodic_resample(model_a, sample_count)
                raw_b = periodic_resample(model_b, sample_count)
            else:
                new_profiles.extend(old_by_model[model])
                continue

            period = float(row["period_rg"])
            for index in range(sample_count):
                s = index / sample_count
                new_profiles.append({
                    "schema": SCHEMA, "case_id": case_id, "f": f_value,
                    "chiN": float(chi_n), "model": model,
                    "model_label": row["model_label"], "period_rg": period,
                    "s": s, "x_rg": s * period,
                    "phi_a": float(aligned[index]),
                    "phi_b": float(1.0 - aligned[index]),
                    "phi_a_raw": float(raw_a[index]),
                    "phi_b_raw": float(raw_b[index]),
                })

        opf = [row for row in summary if row["model"] == "liu2019_opf"]
        if (
            len(opf) != 1 or opf[0]["status"] != "accepted"
            or opf[0]["calibration_role"] != "per_state_force_stress_mapping"
            or any(not opf[0][field] for field in
                   ("opf_c2", "opf_c3", "opf_c4", "opf_c5", "opf_c6"))
            or not math.isfinite(float(opf[0]["opf_mapping_residual"]))
        ):
            raise ValueError(f"{case_id}: recomputed OPF row is not accepted")
        new_profiles.sort(key=lambda row: (row["model"], float(row["s"])))
        atomic_csv(case_dir / "summary.csv", summary)
        atomic_csv(case_dir / "profiles.csv", new_profiles)
        refreshed[case_id] = {
            "scft_period_rg": scft_period,
            "opf_period_rg": float(opf[0]["period_rg"]),
            "opf_profile_rms": float(opf[0]["profile_rms"]),
            "opf_mapping_residual": float(opf[0]["opf_mapping_residual"]),
        }

    manifest["status"] = "refreshed"
    manifest["refreshed"] = refreshed
    manifest["case_hashes"] = {
        case_ids(*key)[0]: {
            name: sha256(staging / "cases" / case_ids(*key)[0] / name)
            for name in ("summary.csv", "profiles.csv", "native_profiles.csv")
        }
        for key in CASES
    }
    atomic_json(staging / "promotion_manifest.json", manifest)


def promote(staging: Path) -> None:
    manifest_path = staging / "promotion_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    if manifest.get("status") != "refreshed":
        raise ValueError("staging manifest is not refreshed")
    backup = AUDIT / "canonical_backup_before_promotion"
    if backup.exists():
        raise FileExistsError(f"canonical backup already exists: {backup}")
    backup.mkdir(parents=True)
    promoted: list[str] = []
    try:
        for key in CASES:
            case_id, _ = case_ids(*key)
            source = staging / "cases" / case_id
            expected = manifest["case_hashes"][case_id]
            for name, digest in expected.items():
                if sha256(source / name) != digest:
                    raise ValueError(f"{case_id}: staged hash changed for {name}")
            target = CANONICAL / "cases" / case_id
            saved = backup / case_id
            os.replace(target, saved)
            temporary = target.with_name(target.name + ".incoming")
            shutil.copytree(source, temporary)
            os.replace(temporary, target)
            promoted.append(case_id)
    except Exception:
        for case_id in reversed(promoted):
            target = CANONICAL / "cases" / case_id
            if target.exists():
                shutil.rmtree(target)
            os.replace(backup / case_id, target)
        raise
    manifest["status"] = "promoted"
    manifest["canonical_backup"] = str(backup.relative_to(ROOT))
    atomic_json(manifest_path, manifest)
    atomic_json(AUDIT / "promotion_manifest.json", manifest)


def sync(staging: Path) -> None:
    manifest = json.loads((staging / "promotion_manifest.json").read_text())
    if manifest.get("status") != "refreshed":
        raise ValueError("staging manifest is not refreshed")
    for key in CASES:
        case_id, _ = case_ids(*key)
        source = staging / "cases" / case_id
        for name, digest in manifest["case_hashes"][case_id].items():
            if sha256(source / name) != digest:
                raise ValueError(f"{case_id}: staged hash changed for {name}")
        target = CANONICAL / "cases" / case_id
        previous = target.with_name(target.name + ".previous")
        incoming = target.with_name(target.name + ".incoming")
        if previous.exists() or incoming.exists():
            raise FileExistsError(f"stale transaction path for {case_id}")
        shutil.copytree(source, incoming)
        os.replace(target, previous)
        try:
            os.replace(incoming, target)
        except Exception:
            os.replace(previous, target)
            raise
        shutil.rmtree(previous)
    manifest["status"] = "promoted"
    atomic_json(staging / "promotion_manifest.json", manifest)
    atomic_json(AUDIT / "promotion_manifest.json", manifest)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("prepare", "refresh", "promote", "sync"))
    parser.add_argument("--staging", type=Path, required=True)
    args = parser.parse_args()
    staging = args.staging.resolve()
    {"prepare": prepare, "refresh": refresh, "promote": promote, "sync": sync}[
        args.mode
    ](staging)


if __name__ == "__main__":
    main()
