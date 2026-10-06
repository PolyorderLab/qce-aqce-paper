from __future__ import annotations

import csv
import importlib.util
import json
import math
import subprocess
import sys
from pathlib import Path

import pytest


PROJECT = Path(__file__).resolve().parents[1]


def load_script(name: str):
    path = PROJECT / "scripts" / name
    spec = importlib.util.spec_from_file_location(path.stem, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    sys.path.insert(0, str(path.parent))
    try:
        spec.loader.exec_module(module)
    finally:
        sys.path.remove(str(path.parent))
    return module


def test_liu2019_digitization_is_traceable_and_reproducible() -> None:
    module = load_script("digitize_liu2019_figure2.py")
    rows = module.digitize()
    assert len(rows) == 203
    assert {row["source_panel"] for row in rows} == {"b", "c"}
    assert {row["visibility_status"] for row in rows} == {"visible_marker"}
    assert {row["source_sha256"] for row in rows} == {module.EXPECTED_SHA256}
    symmetric_opf = [
        row
        for row in rows
        if row["model"] == "liu2019_opf_mapped_published"
        and row["f"] == "0.50"
    ]
    assert len(symmetric_opf) == 25
    assert [row["chiN"] for row in symmetric_opf] == list(range(11, 36))
    assert abs(float(symmetric_opf[-1]["signed_period_error"]) - 0.01075) < 8e-4


def test_allmodel_pilot_profiles_reproduce_reported_rms() -> None:
    outdir = PROJECT / "results" / "liu2019_stress_free_period_map_allmodel_pilot"
    summary_path = outdir / "summary.csv"
    profiles_path = outdir / "profiles.csv"
    if not summary_path.exists() or not profiles_path.exists():
        return
    with summary_path.open(newline="", encoding="utf-8") as handle:
        summary = {row["model"]: row for row in csv.DictReader(handle)}
    with profiles_path.open(newline="", encoding="utf-8") as handle:
        profiles: dict[str, list[float]] = {}
        for row in csv.DictReader(handle):
            profiles.setdefault(row["model"], []).append(float(row["phi_a"]))
    reference = profiles["scft"]
    assert len(reference) == 256
    for model, values in profiles.items():
        rms = (sum((value - target) ** 2 for value, target in zip(values, reference)) / 256) ** 0.5
        assert abs(rms - float(summary[model]["profile_rms"])) < 2e-12


def test_polyorder_protocol_pilot_records_requested_resolution_and_solvers() -> None:
    outdir = (
        PROJECT
        / "results"
        / "liu2019_stress_free_period_map_polyorder_protocol_pilot"
    )
    summary_path = outdir / "summary.csv"
    native_path = outdir / "native_profiles.csv"
    assert summary_path.exists()
    assert native_path.exists()
    with summary_path.open(newline="", encoding="utf-8") as handle:
        row = next(csv.DictReader(handle))
    assert row["schema"] == "liu2019-stress-free-period-map-v2"
    assert row["model"] == "scft"
    assert row["status"] == "accepted"
    assert float(row["scft_residual_norm"]) < 1e-6
    assert 0.08 <= float(row["scft_actual_spacing_rg"]) <= 0.12
    assert float(row["scft_requested_spacing_rg"]) == 0.1
    assert float(row["scft_contour_step"]) == 0.01
    assert row["scft_field_updater"] == (
        "Anderson(SD(0.1); m=10, warmup=100, αw=0.2)"
    )
    assert row["scft_cell_updater"] == (
        "VariableCell(BB(1.0), field_updater)"
    )
    with native_path.open(newline="", encoding="utf-8") as handle:
        native_rows = list(csv.DictReader(handle))
    assert len(native_rows) == int(row["grid_count"]) == 40


def test_reproduction_workflow_is_discoverable() -> None:
    runbook = (PROJECT / "docs" / "liu2019_stress_free_period_map.md").read_text(
        encoding="utf-8"
    )
    entry_points = (
        "scripts/write_liu2019_stress_free_period_map.jl",
        "scripts/digitize_liu2019_figure2.py",
        "scripts/validate_liu2019_stress_free_period_map.py",
        "scripts/render_liu2019_stress_free_period_map.py",
    )
    for entry_point in entry_points:
        assert entry_point in runbook
        assert (PROJECT / entry_point).exists()

    campaign_source = (
        PROJECT / "scripts" / "write_liu2019_stress_free_period_map.jl"
    ).read_text(encoding="utf-8")
    assert "const LPC_USAGE" in campaign_source
    assert "--reuse-existing-scft" in campaign_source
    assert "--assemble-only" in campaign_source
    assert "--max-chiN=35" in campaign_source

    for script, expected_option in (
        ("digitize_liu2019_figure2.py", "--outdir"),
        ("validate_liu2019_stress_free_period_map.py", "--results-dir"),
        ("render_liu2019_stress_free_period_map.py", "--results-dir"),
    ):
        completed = subprocess.run(
            [sys.executable, str(PROJECT / "scripts" / script), "--help"],
            check=True,
            capture_output=True,
            text=True,
        )
        assert expected_option in completed.stdout


def test_fixed_stiffness_period_diagnostic_is_complete_and_gated() -> None:
    module = load_script("render_macromolecules_figures.py")
    fixed = module.fixed_stiffness_period_rows()
    assert len(fixed) == 146
    assert {row["model"] for row in fixed} == {"bvk2_fixed"}
    assert {row["model_label"] for row in fixed} == {"BVK, fixed stiffness"}
    assert {row["grid_count"] for row in fixed} == {"256"}
    assert len({(float(row["f"]), float(row["chiN"])) for row in fixed}) == 146
    assert abs(
        100.0 * sum(float(row["abs_period_error"]) for row in fixed) / 146
        - 3.643018
    ) < 5.0e-4
    assert abs(
        sum(float(row["profile_rms"]) for row in fixed) / 146 - 0.01282003
    ) < 5.0e-6


def test_profile_rms_figure_includes_fixed_stiffness_panel() -> None:
    svg = (
        PROJECT
        / "results/liu2019_stress_free_period_map/"
        "liu2019_stress_free_profile_comparison.svg"
    ).read_text(encoding="utf-8")
    labels = ("UD", "OPF", "BURP", "QCE", "AQCE")
    for label in labels:
        assert f"<!-- {label} -->" in svg
    assert "<!-- PRM -->" not in svg
    assert "<!-- APRM -->" not in svg
    for panel_label in ("(a)", "(b)", "(c)", "(d)", "(e)"):
        assert f"<!-- {panel_label} -->" in svg


def test_fixed_stiffness_low_fa_resolution_audit_passes() -> None:
    path = (
        PROJECT
        / "results/bvk_fixed_low_fa_profile_rms_audit/resolution_summary.csv"
    )
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    expected = {
        (0.20, 25.0), (0.20, 26.0), (0.20, 27.0),
        (0.25, 26.0), (0.25, 27.0), (0.25, 28.0),
        (0.30, 28.0), (0.30, 29.0), (0.30, 30.0),
    }
    assert len(rows) == 9
    assert {(float(row["f"]), float(row["chiN"])) for row in rows} == expected
    assert {row["status"] for row in rows} == {"accepted"}
    assert {row["resolution_status"] for row in rows} == {"accepted"}
    assert max(float(row["relative_period_change_percent"]) for row in rows) < 0.156
    assert max(float(row["absolute_profile_rms_change"]) for row in rows) < 0.00125


def test_fixed_stiffness_production_map_passes() -> None:
    result_dir = PROJECT / "results/bvk_fixed_stiffness_period_map_nx256"
    with (result_dir / "summary.csv").open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    report = json.loads((result_dir / "validation_report.json").read_text())

    assert len(rows) == 146
    assert len({(float(row["f"]), float(row["chiN"])) for row in rows}) == 146
    assert {row["schema"] for row in rows} == {
        "bvk2-fixed-stiffness-period-map-v1"
    }
    assert {row["status"] for row in rows} == {"accepted"}
    assert {row["adaptive"] for row in rows} == {"false"}
    assert {int(row["grid_count"]) for row in rows} == {256}
    assert {int(row["oversample_factor"]) for row in rows} == {2}
    assert {float(row["sensor_filter_ratio"]) for row in rows} == {12.0}
    for row in rows:
        for gate in (
            "field_gate_pass",
            "composition_gate_pass",
            "morphology_gate_pass",
            "cell_gate_pass",
            "local_minimum_check_pass",
            "stress_orientation_pass",
            "accepted",
        ):
            assert row[gate] == "true"
        assert (PROJECT / row["source_root"]).is_file()
        assert (PROJECT / row["source_profile"]).is_file()

    assert report["status"] == "passed"
    assert report["state_count"] == 146
    assert report["profile_row_count"] == 146 * 1024
    assert report["mean_abs_period_error_percent"] == pytest.approx(
        3.6430175717, abs=1.0e-9
    )
    assert report["mean_profile_rms"] == pytest.approx(0.0128200277, abs=1.0e-10)
    assert report["common_cohort_state_count"] == 133
    assert report["common_cohort_mean_abs_period_error_percent"] == pytest.approx(
        3.6164403888, abs=1.0e-9
    )
    assert report["common_cohort_mean_profile_rms"] == pytest.approx(
        0.0131030183, abs=1.0e-10
    )
    assert max(float(row["projected_force_norm"]) for row in rows) < 1.0e-5
    assert max(float(row["projected_force_maxabs"]) for row in rows) < 2.0e-5
    assert max(abs(float(row["cell_stress"])) for row in rows) < 1.0e-6


def test_fixed_stiffness_period_diagnostic_fails_closed(tmp_path: Path) -> None:
    module = load_script("render_macromolecules_figures.py")
    source = PROJECT / "results/bvk_fixed_stiffness_period_map_nx256/summary.csv"
    with source.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
        fields = list(rows[0])
    first_fixed = rows[0]
    first_fixed["field_gate_pass"] = "false"
    malformed = tmp_path / "candidate_rows.csv"
    with malformed.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    with pytest.raises(ValueError, match="fails a required gate"):
        module.fixed_stiffness_period_rows(malformed)


def test_validator_accepts_only_resolved_nonaccepted_ud_classifications() -> None:
    module = load_script("validate_liu2019_stress_free_period_map.py")
    provisional = {
        "model": "uneyama_doi",
        "f": "0.20",
        "chiN": "27",
        "status": "provisional",
        "status_reason": "one_or_more_model_gates_failed",
        "reduced_field_protocol": module.REDUCED_FIELD_PROTOCOLS["uneyama_doi"],
        "field_gate_pass": "true",
        "cell_gate_pass": "false",
        "resolution_gate_pass": "true",
        "composition_gate_pass": "true",
        "morphology_gate_pass": "true",
        "period_local_minimum_check_pass": "false",
        "period_boundary_limited": "false",
        "branch_dominant_mode": "1",
    }
    assert module.is_classified_nonaccepted_ud(provisional)

    rejected = {
        "model": "uneyama_doi",
        "f": "0.50",
        "chiN": "11",
        "status": "rejected",
        "status_reason": "solver_failure",
        "error_message": (
            "all stress-free Uneyama-Doi period optimization evaluations failed"
        ),
    }
    assert module.is_classified_nonaccepted_ud(rejected)

    unresolved = dict(provisional, chiN="28")
    assert not module.is_classified_nonaccepted_ud(unresolved)

    no_minimum = dict(
        unresolved,
        status_reason="no_two_sided_primitive_period_minimum",
    )
    assert module.is_classified_nonaccepted_ud(no_minimum)


def test_polynomial_native_profile_gates_detect_topology_and_resolution() -> None:
    module = load_script("validate_liu2019_stress_free_period_map.py")
    count = 128
    mean = 0.2
    primitive = [
        mean + 0.15 * math.cos(2.0 * math.pi * index / count)
        for index in range(count)
    ]
    repeated = [
        mean + 0.15 * math.cos(4.0 * math.pi * index / count)
        for index in range(count)
    ]
    unresolved = [
        value + 0.01 * math.cos(2.0 * math.pi * (count // 2) * index / count)
        for index, value in enumerate(primitive)
    ]

    assert module.periodic_domain_count(primitive) == 1
    assert module.dominant_nonzero_mode(primitive, mean) == 1
    assert (
        module.spectral_tail_rms_fraction(primitive)
        < module.POLYNOMIAL_SPECTRAL_TAIL_RMS_TOL
    )
    assert module.periodic_domain_count(repeated) == 2
    assert module.dominant_nonzero_mode(repeated, mean) == 2
    assert (
        module.spectral_tail_rms_fraction(unresolved)
        > module.POLYNOMIAL_SPECTRAL_TAIL_RMS_TOL
    )


def test_native_profile_inventory_rejects_duplicate_and_missing_indices() -> None:
    module = load_script("validate_liu2019_stress_free_period_map.py")

    def row(index: int, count: int = 4) -> dict[str, str]:
        return {
            "native_index": str(index),
            "native_grid_count": str(count),
            "s": str((index - 1) / count),
            "phi_a": "0.4",
            "phi_b": "0.6",
        }

    assert module.native_profile_inventory_errors(
        [row(3), row(1), row(4), row(2)], 4
    ) == []

    errors = module.native_profile_inventory_errors(
        [row(1), row(2), row(2), row(4)], 4
    )
    assert "duplicate native_index values: [2]" in errors
    assert "native_index inventory is not the contiguous range 1:4" in errors

    errors = module.native_profile_inventory_errors([row(1), row(2), row(3)], 4)
    assert "has 3 rows; expected 4" in errors
    assert "native_index inventory is not the contiguous range 1:4" in errors

    malformed = [row(1), row(2), row(3), row(4)]
    malformed[2]["native_index"] = "3.5"
    errors = module.native_profile_inventory_errors(malformed, 4)
    assert "row 3 has an invalid native_index" in errors
    assert "native_index inventory is not the contiguous range 1:4" in errors


def test_native_profile_inventory_rejects_bad_grid_coordinates_and_values() -> None:
    module = load_script("validate_liu2019_stress_free_period_map.py")
    rows = [
        {
            "native_index": str(index),
            "native_grid_count": "4",
            "s": str((index - 1) / 4),
            "phi_a": "0.4",
            "phi_b": "0.6",
        }
        for index in range(1, 5)
    ]
    rows[0]["native_grid_count"] = "5"
    rows[1]["s"] = "nan"
    rows[2]["phi_a"] = "inf"
    rows[3]["s"] = "0.6"

    errors = module.native_profile_inventory_errors(rows, 4)
    assert "native_index 1 records native_grid_count 5; expected 4" in errors
    assert "native_index 2 has a nonfinite coordinate" in errors
    assert "native_index 3 has a nonfinite phi_a" in errors
    assert "native_index 4 coordinate is 0.6; expected 0.75" in errors


def test_canonical_opf_ok_native_profiles_pass_publication_gates() -> None:
    module = load_script("validate_liu2019_stress_free_period_map.py")
    rows = module.read_csv(
        PROJECT / "results" / "liu2019_stress_free_period_map",
        "native_profiles.csv",
    )
    profiles: dict[tuple[tuple[float, float], str], list[dict[str, str]]] = {}
    for row in rows:
        if row["model"] not in module.POLYNOMIAL_PROFILE_MODELS:
            continue
        profiles.setdefault((module.case_key(row), row["model"]), []).append(row)

    assert len(profiles) == 2 * 146
    maximum_tail = {model: 0.0 for model in module.POLYNOMIAL_PROFILE_MODELS}
    for (state, model), profile_rows in profiles.items():
        ordered = sorted(profile_rows, key=lambda row: int(row["native_index"]))
        values = [float(row["phi_a"]) for row in ordered]
        assert len(values) == 128
        assert module.native_profile_inventory_errors(profile_rows, 128) == []
        assert module.periodic_domain_count(values) == 1
        assert module.dominant_nonzero_mode(values, state[0]) == 1
        tail = module.spectral_tail_rms_fraction(values)
        maximum_tail[model] = max(maximum_tail[model], tail)
        assert tail <= module.POLYNOMIAL_SPECTRAL_TAIL_RMS_TOL

    assert maximum_tail["liu2019_opf"] < 5.0e-9
    assert maximum_tail["ohta_kawasaki"] < 3.0e-9


def test_promoted_low_fa_scft_references_and_opf_resolution_audit_pass() -> None:
    canonical = PROJECT / "results" / "liu2019_stress_free_period_map"
    expected_periods = {
        25: 3.3772815752408936,
        26: 3.6724051029251457,
        27: 3.8282917974818926,
    }
    for chi_n, expected_period in expected_periods.items():
        case = canonical / "cases" / f"f0.2_chiN{chi_n}"
        with (case / "summary.csv").open(newline="", encoding="utf-8") as handle:
            summary = list(csv.DictReader(handle))
        scft = next(row for row in summary if row["model"] == "scft")
        opf = next(row for row in summary if row["model"] == "liu2019_opf")
        assert scft["status"] == "accepted"
        assert float(scft["period_rg"]) == pytest.approx(expected_period, abs=1e-12)
        assert float(scft["scft_requested_spacing_rg"]) == 0.05
        assert float(scft["scft_contour_step"]) == 0.005
        assert 0.04 <= float(scft["scft_actual_spacing_rg"]) <= 0.06
        assert int(scft["grid_count"]) >= 64
        assert opf["status"] == "accepted"
        assert float(opf["scft_period_rg"]) == pytest.approx(expected_period, abs=1e-12)

    audit = PROJECT / "results" / "scft_low_fa_low_chi_promotion"
    report = json.loads((audit / "opf_nx256_validation.json").read_text())
    assert report["status"] == "passed"
    assert report["case_count"] == report["accepted_case_count"] == 3
    assert report["failures"] == []
    assert report["max_relative_period_change"] < 1e-6
    assert report["max_profile_rms_metric_change"] < 1e-4
    assert report["max_aligned_native_profile_rms"] < 1e-4

    with (audit / "opf_nx256_comparison.csv").open(
        newline="", encoding="utf-8"
    ) as handle:
        comparison = list(csv.DictReader(handle))
    assert len(comparison) == 3
    assert {int(row["chiN"]) for row in comparison} == {25, 26, 27}
    assert {row["accepted"] for row in comparison} == {"true"}
