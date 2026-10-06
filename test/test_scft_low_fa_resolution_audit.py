from __future__ import annotations

import csv
import json
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
AUDIT = PROJECT / "results" / "scft_low_fa_low_chi_resolution_audit"


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def test_low_composition_scft_resolution_verdict() -> None:
    comparison = rows(AUDIT / "comparison.csv")
    assert len(comparison) == 9
    assert len({(float(row["f"]), int(row["chiN"])) for row in comparison}) == 9

    accepted = {
        (float(row["f"]), int(row["chiN"]))
        for row in comparison
        if row["accepted"].lower() == "true"
    }
    assert accepted == {
        (0.25, 19), (0.25, 20), (0.25, 21),
        (0.30, 15), (0.30, 16), (0.30, 17),
    }
    rejected = [row for row in comparison if row["accepted"].lower() == "false"]
    assert {(float(row["f"]), int(row["chiN"])) for row in rejected} == {
        (0.20, 25), (0.20, 26), (0.20, 27)
    }
    assert min(float(row["relative_period_change"]) for row in rejected) > 0.0025

    report = json.loads((AUDIT / "validation_report.json").read_text())
    assert report["status"] == "failed"
    assert report["case_count"] == 9
    assert report["accepted_case_count"] == 6
    assert len(report["failures"]) == 3


def test_fA_0p20_error_is_contour_dominated_and_refined_results_converge() -> None:
    decomposition = rows(AUDIT / "resolution_decomposition.csv")
    assert len(decomposition) == 3
    assert {int(row["chiN"]) for row in decomposition} == {25, 26, 27}
    assert max(
        float(row["spatial_only_relative_period_change"])
        for row in decomposition
    ) < 1.5e-4
    assert min(
        float(row["contour_only_relative_period_change"])
        for row in decomposition
    ) > 2.5e-3
    assert max(
        float(row["combined_vs_full_finer_relative_period_change"])
        for row in decomposition
    ) < 1.0e-3
    assert max(
        float(row["combined_vs_full_finer_profile_rms"])
        for row in decomposition
    ) < 1.0e-3
    assert all(row["full_finer_accepted"].lower() == "true" for row in decomposition)
