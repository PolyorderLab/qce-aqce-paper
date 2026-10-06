from __future__ import annotations

import csv
import importlib.util
import inspect
from pathlib import Path
import re
import sys

import pytest


PROJECT = Path(__file__).resolve().parents[1]
SCRIPT = PROJECT / "scripts" / "render_macromolecules_figures.py"
sys.path.insert(0, str(SCRIPT.parent))
SPEC = importlib.util.spec_from_file_location("render_macromolecules_figures", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def accepted_rows() -> list[dict[str, str]]:
    source = PROJECT / "results/liu2019_stress_free_period_map/summary.csv"
    with source.open(newline="", encoding="utf-8") as handle:
        return [row for row in csv.DictReader(handle) if row["status"] == "accepted"]


def test_ok_strong_extension_requires_explicit_extrapolation() -> None:
    source = (
        PROJECT / "scripts/write_liu2019_stress_free_period_map.jl"
    ).read_text(encoding="utf-8")
    assert '--allow-ok-extrapolation", "false"' in source
    assert "OK states above chiN=35 require --allow-ok-extrapolation=true" in source
    assert "check_domain=!allow_ok_extrapolation" in source
    assert 'model == "ohta_kawasaki" ? (0.1, 0.2, 0.5)' in source
    assert '--ok-gradient-tol", "8e-6"' in source
    assert "0.0 < ok_gradient_tolerance <= 8.0e-6" in source
    assert "liu2019_ok_asymptotic_quadratic_extrapolated_c3_c4" in source


def test_aggregate_cohort_adds_fixed_stiffness_to_locked_intersection() -> None:
    subsets = MODULE.matched_lamellar_aggregate_subsets_with_fixed(accepted_rows())

    assert tuple(subsets) == MODULE.LAMELLAR_AGGREGATE_PLOT_MODELS
    assert len(MODULE.LAMELLAR_AGGREGATE_EXPECTED_STATES) == 133
    assert {len(rows) for rows in subsets.values()} == {133}
    for rows in subsets.values():
        states = {(float(row["f"]), float(row["chiN"])) for row in rows}
        assert states == MODULE.LAMELLAR_AGGREGATE_EXPECTED_STATES

    expected = {
        "ohta_kawasaki": (15.978250, 10.491931),
        "uneyama_doi": (11.183954, 4.536149),
        "liu2019_opf": (4.309027, 5.789067),
        "burp_ti": (16.529954, 2.236156),
        "bvk2_fixed": (3.616440, 1.310302),
        "bvk2": (1.463605, 1.723532),
    }
    for model, model_rows in subsets.items():
        period_error = 100 * sum(
            float(row["abs_period_error"]) for row in model_rows
        ) / len(model_rows)
        profile_error = 100 * sum(
            float(row["profile_rms"]) for row in model_rows
        ) / len(model_rows)
        assert (period_error, profile_error) == pytest.approx(
            expected[model], abs=5e-7
        )


def test_strong_segregation_extension_is_complete_and_separate() -> None:
    extension = MODULE.strong_extension_plot_rows()
    aggregate = MODULE.strong_extension_aggregate_rows()

    assert len(extension) == 735
    assert len(MODULE.STRONG_EXTENSION_STATES) == 105
    assert {
        (float(row["f"]), float(row["chiN"])) for row in extension
    } == MODULE.STRONG_EXTENSION_STATES
    assert tuple(row["model"] for row in aggregate) == (
        "ohta_kawasaki",
        "uneyama_doi",
        "liu2019_opf",
        "burp_ti",
        "bvk2_fixed",
        "bvk2",
    )
    assert tuple(row["model"] for row in aggregate) == (
        *MODULE.LAMELLAR_AGGREGATE_PLOT_MODELS,
    )
    expected = {
        "ohta_kawasaki": (22.057585, 0.129004),
        "liu2019_opf": (6.381167, 0.099438),
        "uneyama_doi": (14.966467, 0.018547),
        "burp_ti": (18.833016, 0.018490),
        "bvk2_fixed": (5.392299, 0.014049),
        "bvk2": (1.187399, 0.016363),
    }
    for row in aggregate:
        assert (
            float(row["mean_abs_period_error_percent"]),
            float(row["mean_profile_rms"]),
        ) == pytest.approx(expected[row["model"]], abs=5e-7)


def test_figure5_uses_one_axis_with_the_complete_model_order() -> None:
    svg = (
        PROJECT
        / "results/liu2019_stress_free_period_map/aggregate_error_bars.svg"
    ).read_text(encoding="utf-8")
    labels = ("OK", "UD", "OPF", "BURP", "QCE", "AQCE")
    occurrences = [svg.index(f"<!-- {label} -->") for label in labels]
    assert occurrences == sorted(occurrences)
    assert "<!-- PRM -->" not in svg
    assert "<!-- APRM -->" not in svg
    assert all(svg.count(f"<!-- {label} -->") == 1 for label in labels)
    assert "<!-- (a) -->" not in svg
    assert "<!-- (b) -->" not in svg
    for legend_label in (
        "period",
        "profile",
        "weak to intermediate segregation",
        "strong segregation",
    ):
        assert f"<!-- {legend_label} -->" in svg
    for repeated_label in (
        "period, weak to intermediate segregation",
        "profile, weak to intermediate segregation",
        "period, strong segregation",
        "profile, strong segregation",
    ):
        assert repeated_label not in svg
    renderer = inspect.getsource(MODULE.render_lamellar_aggregate_errors)
    assert renderer.index("common[2]") < renderer.index("strong_cohort[2]")
    assert renderer.index("strong_cohort[2]") < renderer.index("common[3]")
    assert renderer.index("common[3]") < renderer.index("strong_cohort[3]")
    assert "ax.legend(" in renderer
    assert 'loc="upper right"' in renderer
    assert 'loc="outside upper center"' not in renderer
    assert 'spines[["top", "right"]].set_visible(False)' not in renderer
    assert "tick_params(top=False, right=False)" not in renderer


def test_removed_six_state_composite_is_not_rendered() -> None:
    assert "macromolecules_figure2.svg" not in inspect.getsource(
        MODULE.render_lamellar
    )


def test_period_error_uses_epsilon_in_renderer_and_figure() -> None:
    label = r"$\epsilon=(D_0-D_{0,\mathrm{SCFT}})/D_{0,\mathrm{SCFT}}$"
    svg = (
        PROJECT / "results/liu2019_stress_free_period_map"
        / "liu2019_stress_free_period_comparison.svg"
    ).read_text(encoding="utf-8")
    for text in (inspect.getsource(MODULE.render_liu2019_stress_free_period_map), svg):
        assert label in text
        assert r"\eta=" not in text


def test_aggregate_cohort_rejects_duplicate_model_state() -> None:
    accepted = accepted_rows()
    duplicate = next(row for row in accepted if row["model"] == "bvk2")

    with pytest.raises(ValueError, match="duplicate accepted row"):
        MODULE.matched_lamellar_aggregate_subsets([*accepted, duplicate.copy()])


def test_aggregate_cohort_rejects_a_shrunken_intersection() -> None:
    accepted = accepted_rows()
    removed = next(
        row
        for row in accepted
        if row["model"] == "bvk2"
        and float(row["f"]) == 0.50
        and float(row["chiN"]) == 35.0
    )

    with pytest.raises(ValueError, match="locked 133-state manifest"):
        MODULE.matched_lamellar_aggregate_subsets(
            [row for row in accepted if row is not removed]
        )


@pytest.mark.parametrize(
    ("metric", "malformed_value"),
    (("abs_period_error", "nan"), ("profile_rms", "inf")),
)
def test_aggregate_cohort_rejects_nonfinite_metrics(
    metric: str, malformed_value: str
) -> None:
    accepted = [row.copy() for row in accepted_rows()]
    row = next(row for row in accepted if row["model"] == "bvk2")
    row[metric] = malformed_value

    with pytest.raises(ValueError, match=rf"{metric} is not finite"):
        MODULE.matched_lamellar_aggregate_subsets(accepted)


@pytest.mark.parametrize(
    ("field", "malformed_value", "message"),
    (
        ("schema", "legacy-map-v1", "schema/claim fingerprint mismatch"),
        ("calibration_role", "unknown", "source fingerprint mismatch"),
        ("period_optimizer", "UnverifiedSolver", "source fingerprint mismatch"),
        ("reduced_field_protocol", "legacy field solve", "source fingerprint mismatch"),
    ),
)
def test_aggregate_cohort_rejects_source_provenance_drift(
    field: str, malformed_value: str, message: str
) -> None:
    accepted = [row.copy() for row in accepted_rows()]
    row = next(row for row in accepted if row["model"] == "bvk2")
    row[field] = malformed_value

    with pytest.raises(ValueError, match=message):
        MODULE.matched_lamellar_aggregate_subsets(accepted)


def test_aggregate_cohort_rejects_missing_required_metric() -> None:
    accepted = [row.copy() for row in accepted_rows()]
    row = next(row for row in accepted if row["model"] == "bvk2")
    del row["profile_rms"]

    with pytest.raises(ValueError, match="lacks required fields"):
        MODULE.matched_lamellar_aggregate_subsets(accepted)
