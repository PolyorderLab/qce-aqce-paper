from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys

import pytest


PROJECT = Path(__file__).resolve().parents[1]
SCRIPT = PROJECT / "scripts/generate_bvk_contribution_claim_evidence.py"
ARTIFACT = PROJECT / "results/bvk_contribution_claim_evidence/claim_evidence.json"
SPEC = importlib.util.spec_from_file_location("bvk_claim_evidence", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def claims_by_id(report: dict[str, object]) -> dict[str, dict[str, object]]:
    return {row["claim_id"]: row for row in report["claims"]}


def mutate_source(
    monkeypatch: pytest.MonkeyPatch,
    target: Path,
    mutation,
) -> None:
    original = MODULE.read_rows

    def altered(path: Path) -> list[dict[str, str]]:
        rows = original(path)
        if path == target:
            mutation(rows)
        return rows

    monkeypatch.setattr(MODULE, "read_rows", altered)


def test_checked_claim_evidence_is_a_fresh_deterministic_build() -> None:
    assert ARTIFACT.read_text(encoding="utf-8") == MODULE.serialized_evidence()


@pytest.mark.parametrize("label", ("QCE", "PRM"))
def test_frozen_fixed_stiffness_labels_remain_compatible(label: str) -> None:
    MODULE.require_current_or_legacy_label(
        [{"model_label": label}],
        "model_label",
        "QCE",
        frozenset({"PRM"}),
        "fixed-stiffness map",
    )


def test_mixed_current_and_legacy_labels_are_rejected() -> None:
    with pytest.raises(ValueError, match="expected one of"):
        MODULE.require_current_or_legacy_label(
            [{"model_label": "QCE"}, {"model_label": "PRM"}],
            "model_label",
            "QCE",
            frozenset({"PRM"}),
            "fixed-stiffness map",
        )


def test_cohorts_and_evidence_classes_are_locked() -> None:
    report = MODULE.build_evidence()
    cohorts = report["cohorts"]

    assert cohorts["common_133"]["state_count"] == 133
    assert len(cohorts["common_133"]["membership"]) == 133
    assert cohorts["held_out_143"]["state_count"] == 143
    assert cohorts["held_out_143"]["excluded_calibration_case_ids"] == [
        "f0.35_chiN30",
        "f0.5_chiN15",
        "f0.5_chiN25",
    ]
    assert cohorts["full_146_bvk2"]["state_count"] == 146
    assert cohorts["full_146_fixed"]["state_count"] == 146
    assert cohorts["coefficient_diagnostic_full_146"]["evidence_class"] == (
        "oracle_diagnostic"
    )
    assert cohorts["coefficient_diagnostic_full_146"]["resolution_nx"] == 64
    assert cohorts["high_chiN_30_35"]["state_count"] == 42
    assert cohorts["high_chiN_30_35"]["evidence_class"] == "production_control"
    assert cohorts["high_chiN_30_35_diagnostic"]["evidence_class"] == (
        "oracle_diagnostic"
    )
    assert cohorts["aba_7"]["chiN"] == [20.0, 22.0, 25.0, 30.0, 35.0, 40.0, 45.0]
    assert cohorts["cell_stress_2_state"]["resolution_nx"] == 256
    assert cohorts["strong_segregation_105"]["state_count"] == 105
    assert len(cohorts["strong_segregation_105"]["membership"]) == 105
    assert cohorts["strong_segregation_105"]["chiN"] == list(range(36, 51))
    assert cohorts["phase_boundary_149"]["state_count"] == 149
    assert cohorts["phase_boundary_149"]["ordered_phases"] == [
        "BCC", "CYL", "FCC", "GYR", "LAM", "O70"
    ]
    assert all(gate["pass"] for gate in report["gate_audit"].values())
    assert report["diagnostic_selection"]["full_146_oracle_c2"] == 0.06
    assert report["diagnostic_selection"]["frozen_c2"] == 0.16
    assert {row["claim_id"] for row in report["claims"]} == MODULE.expected_claim_ids()


def test_headline_claims_reproduce_canonical_numbers() -> None:
    claims = claims_by_id(MODULE.build_evidence())
    expected = {
        "common133.uneyama_doi.period": 11.183954025075062,
        "common133.uneyama_doi.profile": 0.04536149356528059,
        "common133.bvk2_fixed.period": 3.6164403888156227,
        "common133.bvk2_fixed.profile": 0.013103018273959008,
        "common133.bvk2.period": 1.4636050378236758,
        "common133.bvk2.profile": 0.017235322541175477,
        "held_out_143.period": 1.793064029815808,
        "held_out_143.profile": 0.018946176692923367,
        "full_146_bvk2.period": 1.7719077076960121,
        "full_146_bvk2.profile": 0.018784872703090053,
        "full_146_fixed.period": 3.6430175717205873,
        "full_146_fixed.profile": 0.012820027744768044,
        "coefficient_diagnostic.common133.c2_0.04.period": 2.6662944619803355,
        "coefficient_diagnostic.common133.c2_0.04.profile": 0.010307228174783729,
        "coefficient_diagnostic.common133.c2_0.16.period": 1.4157789043848228,
        "coefficient_diagnostic.common133.c2_0.16.profile": 0.016711512237123256,
        "coefficient_diagnostic.full146.c2_0.00.period": 3.5451191653738556,
        "coefficient_diagnostic.full146.c2_0.00.profile": 0.013744997660387916,
        "coefficient_diagnostic.full146.c2_0.06.period": 2.5500113157037734,
        "coefficient_diagnostic.full146.c2_0.06.profile": 0.01166767113370193,
        "coefficient_diagnostic.full146.c2_0.16.period": 1.7299920899117365,
        "coefficient_diagnostic.full146.c2_0.16.profile": 0.018251913271985592,
        "high_chiN.adaptive_period_wins": 37,
        "high_chiN.fixed.period": 4.952013096365367,
        "high_chiN.adaptive.period": 1.9940519304072608,
        "high_chiN_diagnostic.adaptive_period_wins": 36,
        "high_chiN_diagnostic.fixed.period": 4.750452564410283,
        "high_chiN_diagnostic.adaptive.period": 1.8960290919283738,
        "aba7.bvk2_fixed.period": 3.8311825556241703,
        "aba7.bvk2_fixed.profile": 0.022672651655082292,
        "aba7.bvk2.period": 1.570168150483837,
        "aba7.bvk2.profile": 0.0062082230125418076,
        "aba7.adaptive_period_wins": 7,
        "cell_stress.f0.5_chiN20.bvk2_fixed.period_ratio": 0.9734495119665264,
        "cell_stress.f0.5_chiN20.bvk2.period_ratio": 1.0019366309622748,
        "cell_stress.f0.35_chiN30.bvk2_fixed.period_ratio": 0.9573363979336021,
        "cell_stress.f0.35_chiN30.bvk2.period_ratio": 0.9928910329387016,
        "strong105.state_count": 105,
        "strong105.ohta_kawasaki.period": 22.057584987295538,
        "strong105.ohta_kawasaki.profile": 0.12900406204412226,
        "strong105.liu2019_opf.period": 6.381166799830562,
        "strong105.liu2019_opf.profile": 0.09943754418791283,
        "strong105.uneyama_doi.period": 14.966467453988624,
        "strong105.uneyama_doi.profile": 0.018547289807699655,
        "strong105.burp_ti.period": 18.833015932233092,
        "strong105.burp_ti.profile": 0.018489734352508544,
        "strong105.bvk2_fixed.period": 5.392298722801199,
        "strong105.bvk2_fixed.profile": 0.01404880441425547,
        "strong105.bvk2.period": 1.187399011804396,
        "strong105.bvk2.profile": 0.01636250981874579,
        "phase149.accepted_root_count": 149,
        "phase149.ordered_phase_count": 6,
        "phase149.o70_boundary_root_count": 31,
        "phase149.o70_boundary_family_count": 3,
        "phase149.fcc_boundary_root_count": 12,
        "phase149.fcc_boundary_family_count": 2,
    }
    for claim_id, value in expected.items():
        assert claims[claim_id]["value"] == pytest.approx(value, abs=1.0e-14)

    assert claims["high_chiN.fA_0.20.adaptive_period_wins"]["value"] == 1
    for f_a in (0.25, 0.30, 0.35, 0.40, 0.45, 0.50):
        assert claims[f"high_chiN.fA_{f_a:.2f}.adaptive_period_wins"]["value"] == 6


def test_artifact_is_machine_readable_and_has_no_nonfinite_number() -> None:
    report = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    assert report["schema"] == "bvk-contribution-claim-evidence-v1"
    assert len(report["claims"]) == 88
    assert "NaN" not in ARTIFACT.read_text(encoding="utf-8")


def test_new_source_hashes_are_exact() -> None:
    report = MODULE.build_evidence()
    hashes = {row["path"]: row["sha256"] for row in report["sources"]}
    assert hashes[MODULE.STRONG_MAP.as_posix()] == "ca9b605df17673344ee93b48b2f2e04aba25feecb6d5cd0883b11dcff8d2fc66"
    assert hashes[MODULE.STRONG_FIXED.as_posix()] == "c65d57a2cb39583507adbb39d5a110d3a75eaaf342f61f5bca870804f5e58fc6"
    assert hashes[MODULE.STRONG_AGGREGATE.as_posix()] == "f54c7cb930acca1ad35e8673ed79eb6a51d59faf0cd17b2dcf1819a093385398"
    assert hashes[MODULE.STRONG_VALIDATION.as_posix()] == "3bac9fb4b59497f9ea3df35a91b7295cbecc5ff1a09a3e3d02d7fc39ef5df8f4"
    assert hashes[MODULE.PHASE_BOUNDARIES.as_posix()] == "e1bd7408c4c65a3f6d001f7515557b75904676198c86dcb7c81553f1e700c454"
    assert hashes[MODULE.SCFT_PHASE_REFERENCE.as_posix()] == "c2130bd398e7342c73536b6c5069b1ad371ce6bf3d81d367b8abc7a4f33d5620"


@pytest.mark.parametrize(
    ("target", "field"),
    (
        (MODULE.STRONG_MAP, "resolution_gate_pass"),
        (MODULE.STRONG_FIXED, "projected_force_norm"),
        (MODULE.STRONG_AGGREGATE, "mean_profile_rms"),
        (MODULE.PHASE_BOUNDARIES, "acceptance_basis"),
        (MODULE.SCFT_PHASE_REFERENCE, "role"),
    ),
)
def test_new_csv_sources_reject_missing_columns(
    monkeypatch: pytest.MonkeyPatch, target: Path, field: str
) -> None:
    mutate_source(monkeypatch, target, lambda rows: [row.pop(field) for row in rows])
    with pytest.raises(ValueError, match="lacks required columns"):
        MODULE.build_evidence()


def test_strong_state_manifest_and_gate_mutations_fail_closed(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    mutate_source(monkeypatch, MODULE.STRONG_MAP, lambda rows: rows[0].__setitem__("chiN", "99.0"))
    with pytest.raises(ValueError, match="state manifest failed|pair cardinality failed"):
        MODULE.build_evidence()


@pytest.mark.parametrize(
    ("field", "value", "message"),
    (
        ("schema", "legacy", "schema contract failed"),
        ("model_label", "wrong", "model_label contract failed"),
        ("field_gate_pass", "false", "fails gates"),
        ("abs_period_error", "nan", "nonfinite abs_period_error"),
    ),
)
def test_strong_source_contract_mutations_fail_closed(
    monkeypatch: pytest.MonkeyPatch, field: str, value: str, message: str
) -> None:
    mutate_source(
        monkeypatch,
        MODULE.STRONG_MAP,
        lambda rows: rows[0].__setitem__(field, value),
    )
    with pytest.raises(ValueError, match=message):
        MODULE.build_evidence()


def test_strong_aggregate_mutation_fails_closed(monkeypatch: pytest.MonkeyPatch) -> None:
    mutate_source(
        monkeypatch,
        MODULE.STRONG_AGGREGATE,
        lambda rows: rows[0].__setitem__("mean_profile_rms", "0.5"),
    )
    with pytest.raises(ValueError, match="aggregate profile metric drift"):
        MODULE.build_evidence()


def test_strong_validation_mutation_fails_closed(monkeypatch: pytest.MonkeyPatch) -> None:
    original = MODULE.read_json

    def altered(path: Path) -> dict[str, object]:
        payload = original(path)
        if path == MODULE.STRONG_VALIDATION:
            payload["status"] = "failed"
        return payload

    monkeypatch.setattr(MODULE, "read_json", altered)
    with pytest.raises(ValueError, match="validation-report contract failed"):
        MODULE.build_evidence()


def test_phase_manifest_and_finite_mutations_fail_closed(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    mutate_source(
        monkeypatch,
        MODULE.PHASE_BOUNDARIES,
        lambda rows: rows[0].__setitem__("transition", "rogue"),
    )
    with pytest.raises(ValueError, match="transition manifest failed"):
        MODULE.build_evidence()


def test_phase_nonfinite_coordinate_fails_closed(monkeypatch: pytest.MonkeyPatch) -> None:
    mutate_source(
        monkeypatch,
        MODULE.PHASE_BOUNDARIES,
        lambda rows: rows[0].__setitem__("fA", "nan"),
    )
    with pytest.raises(ValueError, match="nonfinite fA"):
        MODULE.build_evidence()


def test_phase_acceptance_gate_fails_closed(monkeypatch: pytest.MonkeyPatch) -> None:
    mutate_source(
        monkeypatch,
        MODULE.PHASE_BOUNDARIES,
        lambda rows: rows[0].__setitem__("status", "provisional"),
    )
    with pytest.raises(ValueError, match="status contract failed"):
        MODULE.build_evidence()


def test_scft_phase_reference_manifest_fails_closed(monkeypatch: pytest.MonkeyPatch) -> None:
    mutate_source(
        monkeypatch,
        MODULE.SCFT_PHASE_REFERENCE,
        lambda rows: rows[0].__setitem__("transition", "rogue"),
    )
    with pytest.raises(ValueError, match="transition manifest failed"):
        MODULE.build_evidence()


def test_generator_rejects_duplicate_source_rows(monkeypatch: pytest.MonkeyPatch) -> None:
    original = MODULE.read_rows

    def duplicated(path: Path) -> list[dict[str, str]]:
        rows = original(path)
        if path == MODULE.PRODUCTION_MAP:
            rows.append(rows[0].copy())
        return rows

    monkeypatch.setattr(MODULE, "read_rows", duplicated)
    with pytest.raises(ValueError, match="duplicate|pair cardinality"):
        MODULE.build_evidence()


def test_generator_rejects_failed_production_gate(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    original = MODULE.read_rows

    def failed_gate(path: Path) -> list[dict[str, str]]:
        rows = original(path)
        if path == MODULE.FIXED_MAP:
            rows[0]["cell_gate_pass"] = "false"
        return rows

    monkeypatch.setattr(MODULE, "read_rows", failed_gate)
    with pytest.raises(ValueError, match="fails gates"):
        MODULE.build_evidence()


def test_generator_rejects_common_cohort_membership_drift(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    original = MODULE.read_rows

    def missing_common_state(path: Path) -> list[dict[str, str]]:
        rows = original(path)
        if path == MODULE.PRODUCTION_MAP:
            return [
                row
                for row in rows
                if not (
                    row["model"] == "uneyama_doi"
                    and row["f"] == "0.5"
                    and row["chiN"] == "35.0"
                )
            ]
        return rows

    monkeypatch.setattr(MODULE, "read_rows", missing_common_state)
    with pytest.raises(ValueError, match="locked 133-state manifest|pair cardinality"):
        MODULE.build_evidence()


@pytest.mark.parametrize(
    ("mutation", "message"),
    (
        (
            lambda rows: rows[0].__setitem__("schema", "legacy-production-v1"),
            "schema contract failed",
        ),
        (
            lambda rows: next(row for row in rows if row["model"] == "bvk2").__setitem__(
                "c2", "0.12"
            ),
            "c2 contract failed",
        ),
        (
            lambda rows: next(row for row in rows if row["model"] == "bvk2").__setitem__(
                "grid_count", "128"
            ),
            "grid_count contract failed",
        ),
        (
            lambda rows: next(
                row for row in rows if row["model"] == "uneyama_doi"
            ).__setitem__("model_label", "UD drift"),
            "model_label contract failed",
        ),
    ),
)
def test_production_source_contract_mutations_fail_closed(
    monkeypatch: pytest.MonkeyPatch, mutation, message: str
) -> None:
    mutate_source(monkeypatch, MODULE.PRODUCTION_MAP, mutation)
    with pytest.raises(ValueError, match=message):
        MODULE.build_evidence()


def test_production_source_missing_c2_column_fails_closed(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    mutate_source(
        monkeypatch,
        MODULE.PRODUCTION_MAP,
        lambda rows: [row.pop("c2") for row in rows],
    )
    with pytest.raises(ValueError, match="lacks required columns"):
        MODULE.build_evidence()


@pytest.mark.parametrize(
    ("field", "value", "message"),
    (
        ("schema", "legacy-fixed-v0", "schema contract failed"),
        ("model", "bvk2", "model contract failed"),
        ("status", "provisional", "status contract failed"),
        ("adaptive", "true", "adaptive contract failed"),
        ("c2", "0.12", "c2 contract failed"),
        ("grid_count", "128", "grid_count contract failed"),
    ),
)
def test_fixed_source_contract_mutations_fail_closed(
    monkeypatch: pytest.MonkeyPatch, field: str, value: str, message: str
) -> None:
    mutate_source(monkeypatch, MODULE.FIXED_MAP, lambda rows: rows[0].__setitem__(field, value))
    with pytest.raises(ValueError, match=message):
        MODULE.build_evidence()


def test_fixed_source_missing_required_column_fails_closed(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    mutate_source(
        monkeypatch,
        MODULE.FIXED_MAP,
        lambda rows: [row.pop("adaptive") for row in rows],
    )
    with pytest.raises(ValueError, match="lacks required columns"):
        MODULE.build_evidence()


@pytest.mark.parametrize(
    ("target", "field"),
    (
        (MODULE.COEFFICIENT_DIAGNOSTIC, "fingerprint"),
        (MODULE.ABA_ADAPTIVE, "accepted"),
        (MODULE.ABA_FIXED, "adaptive"),
        (MODULE.CELL_STRESS, "nx"),
    ),
)
def test_each_source_rejects_a_missing_required_column(
    monkeypatch: pytest.MonkeyPatch, target: Path, field: str
) -> None:
    mutate_source(
        monkeypatch,
        target,
        lambda rows: [row.pop(field) for row in rows],
    )
    with pytest.raises(ValueError, match="lacks required columns"):
        MODULE.build_evidence()


def test_fixed_source_state_manifest_mutation_fails_closed(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    mutate_source(
        monkeypatch,
        MODULE.FIXED_MAP,
        lambda rows: rows[0].__setitem__("chiN", "999.0"),
    )
    with pytest.raises(ValueError, match="state manifest failed"):
        MODULE.build_evidence()


@pytest.mark.parametrize(
    ("mutation", "message"),
    (
        (lambda rows: rows[0].__setitem__("solver_schema", "legacy"), "solver_schema contract failed"),
        (lambda rows: rows[0].__setitem__("nx", "128"), "nx contract failed"),
        (lambda rows: rows[0].__setitem__("c2", "0.99"), "c2 contract failed"),
        (
            lambda rows: rows[0].__setitem__("selection_eligible", "false"),
            "fails gates",
        ),
        (
            lambda rows: rows[0].__setitem__(
                "fingerprint", rows[0]["fingerprint"].replace("model=BVK2", "model=UD")
            ),
            "fingerprint contract failed",
        ),
        (lambda rows: rows.pop(), "pair cardinality failed"),
    ),
)
def test_diagnostic_contract_mutations_fail_closed(
    monkeypatch: pytest.MonkeyPatch, mutation, message: str
) -> None:
    mutate_source(monkeypatch, MODULE.COEFFICIENT_DIAGNOSTIC, mutation)
    with pytest.raises(ValueError, match=message):
        MODULE.build_evidence()


def test_diagnostic_exact_state_manifest_mutation_fails_closed(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def change_one_state(rows: list[dict[str, str]]) -> None:
        case_id = rows[0]["case_id"]
        for row in rows:
            if row["case_id"] == case_id:
                row["chiN"] = "999.0"

    mutate_source(monkeypatch, MODULE.COEFFICIENT_DIAGNOSTIC, change_one_state)
    with pytest.raises(ValueError, match="state manifest failed"):
        MODULE.build_evidence()


@pytest.mark.parametrize(
    ("mutation", "message"),
    (
        (lambda rows: rows[0].__setitem__("schema", "legacy-aba"), "schema contract failed"),
        (lambda rows: rows[0].__setitem__("model", "rogue"), "model contract failed"),
        (lambda rows: rows[0].__setitem__("accepted", "false"), "accepted contract failed"),
        (lambda rows: next(row for row in rows if row["model"] == "bvk2").__setitem__("c2", "0.12"), "c2 contract failed"),
        (lambda rows: next(row for row in rows if row["model"] == "bvk2").__setitem__("grid", "64"), "grid contract failed"),
        (lambda rows: rows.pop(), "pair cardinality failed"),
    ),
)
def test_aba_adaptive_contract_mutations_fail_closed(
    monkeypatch: pytest.MonkeyPatch, mutation, message: str
) -> None:
    mutate_source(monkeypatch, MODULE.ABA_ADAPTIVE, mutation)
    with pytest.raises(ValueError, match=message):
        MODULE.build_evidence()


def test_aba_adaptive_state_manifest_mutation_fails_closed(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def change_case(rows: list[dict[str, str]]) -> None:
        case_id = rows[0]["case_id"]
        for row in rows:
            if row["case_id"] == case_id:
                row["chiN"] = "999.0"

    mutate_source(monkeypatch, MODULE.ABA_ADAPTIVE, change_case)
    with pytest.raises(ValueError, match="state manifest failed"):
        MODULE.build_evidence()


@pytest.mark.parametrize(
    ("field", "value", "message"),
    (
        ("schema", "legacy-fixed-aba", "schema contract failed"),
        ("model", "bvk2", "model contract failed"),
        ("status", "provisional", "status contract failed"),
        ("adaptive", "true", "adaptive contract failed"),
        ("c2", "0.12", "c2 contract failed"),
        ("nx", "64", "nx contract failed"),
    ),
)
def test_aba_fixed_contract_mutations_fail_closed(
    monkeypatch: pytest.MonkeyPatch, field: str, value: str, message: str
) -> None:
    mutate_source(monkeypatch, MODULE.ABA_FIXED, lambda rows: rows[0].__setitem__(field, value))
    with pytest.raises(ValueError, match=message):
        MODULE.build_evidence()


def test_aba_fixed_state_manifest_mutation_fails_closed(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    mutate_source(
        monkeypatch,
        MODULE.ABA_FIXED,
        lambda rows: rows[0].__setitem__("chiN", "999.0"),
    )
    with pytest.raises(ValueError, match="state manifest failed"):
        MODULE.build_evidence()


@pytest.mark.parametrize(
    ("mutation", "message"),
    (
        (lambda rows: rows[0].__setitem__("schema", "legacy-stress"), "schema contract failed"),
        (lambda rows: rows[0].__setitem__("model", "rogue"), "model contract failed"),
        (lambda rows: rows[0].__setitem__("accepted", "false"), "accepted contract failed"),
        (lambda rows: rows[0].__setitem__("nx", "128"), "model-resolution contract failed"),
        (
            lambda rows: next(row for row in rows if row["model"] == "bvk2").__setitem__(
                "source_fingerprint",
                next(row for row in rows if row["model"] == "bvk2")[
                    "source_fingerprint"
                ].replace("adaptive=true", "adaptive=false"),
            ),
            "adaptive/c2 provenance contract failed",
        ),
        (
            lambda rows: next(
                row for row in rows if row["model"] == "bvk2_fixed"
            ).__setitem__(
                "source_fingerprint",
                next(row for row in rows if row["model"] == "bvk2_fixed")[
                    "source_fingerprint"
                ].replace("c2=0.16", "c2=0.12"),
            ),
            "adaptive/c2 provenance contract failed",
        ),
        (lambda rows: rows.pop(), "pair cardinality failed"),
    ),
)
def test_cell_stress_contract_mutations_fail_closed(
    monkeypatch: pytest.MonkeyPatch, mutation, message: str
) -> None:
    mutate_source(monkeypatch, MODULE.CELL_STRESS, mutation)
    with pytest.raises(ValueError, match=message):
        MODULE.build_evidence()


def test_cell_stress_state_manifest_mutation_fails_closed(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def change_complete_pair(rows: list[dict[str, str]]) -> None:
        case_id = rows[0]["case_id"]
        for row in rows:
            if row["case_id"] == case_id:
                row["chiN"] = "999.0"

    mutate_source(monkeypatch, MODULE.CELL_STRESS, change_complete_pair)
    with pytest.raises(ValueError, match="state manifest failed"):
        MODULE.build_evidence()


def test_claim_referential_integrity_rejects_unknown_references() -> None:
    report = MODULE.build_evidence()
    claims = [row.copy() for row in report["claims"]]
    claims[0]["cohort_id"] = "unknown"
    with pytest.raises(ValueError, match="unknown cohort"):
        MODULE.require_claim_integrity(claims, report["cohorts"], report["sources"])

    claims = [row.copy() for row in report["claims"]]
    claims[0]["sources"] = ["results/not-canonical.csv"]
    with pytest.raises(ValueError, match="unknown sources"):
        MODULE.require_claim_integrity(claims, report["cohorts"], report["sources"])


def test_claim_id_contract_rejects_extra_claim() -> None:
    report = MODULE.build_evidence()
    claims = [row.copy() for row in report["claims"]]
    extra = claims[0].copy()
    extra["claim_id"] = "unexpected.claim"
    claims.append(extra)
    with pytest.raises(ValueError, match="claim-ID contract failed"):
        MODULE.require_claim_integrity(claims, report["cohorts"], report["sources"])


def test_claim_schema_rejects_missing_metric() -> None:
    report = MODULE.build_evidence()
    claims = [row.copy() for row in report["claims"]]
    claims[0].pop("metric")
    with pytest.raises(ValueError, match="field schema failed"):
        MODULE.require_claim_integrity(claims, report["cohorts"], report["sources"])


def test_diagnostic_claim_cannot_be_relabelled_as_production() -> None:
    report = MODULE.build_evidence()
    claims = [row.copy() for row in report["claims"]]
    diagnostic = next(
        row for row in claims if row["evidence_class"] == "oracle_diagnostic"
    )
    diagnostic["evidence_class"] = "production"
    with pytest.raises(ValueError, match="incompatible"):
        MODULE.require_claim_integrity(claims, report["cohorts"], report["sources"])


def test_absolute_output_outside_project_is_supported(tmp_path: Path) -> None:
    output = tmp_path / "nested" / "claim_evidence.json"
    MODULE.write_or_check_output(output, check=False)
    assert output.read_text(encoding="utf-8") == MODULE.serialized_evidence()
    MODULE.write_or_check_output(output, check=True)
