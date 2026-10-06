#!/usr/bin/env python3
"""Checks for the canonical ABA Figures 7 and 8 and their data contract."""

import copy
import importlib.util
import shutil
import tempfile
from pathlib import Path
from xml.etree import ElementTree as ET


PROJECT = Path(__file__).resolve().parents[1]
renderer_source = (
    PROJECT / "scripts/render_macromolecules_figure8.py"
).read_text(encoding="utf-8")
assert "fig.add_subplot(grid[1, 0:2])" in renderer_source
assert "fig.add_subplot(grid[1, 2:4])" in renderer_source
assert "fig.add_subplot(grid[1, 4:6])" in renderer_source
assert "period_ax = fig.add_subplot(grid[2, 0:3])" in renderer_source
assert "rms_ax = fig.add_subplot(grid[2, 3:6])" in renderer_source
output6 = (
    PROJECT
    / "results/bvk2_aba_comprehensive_validation/macromolecules_figure7.svg"
)
output7 = (
    PROJECT
    / "results/bvk2_aba_comprehensive_validation/macromolecules_figure8.svg"
)
for output in (output6, output7):
    root = ET.parse(output).getroot()
    assert root.tag.endswith("svg")
    text = output.read_text(encoding="utf-8")
    assert "Matplotlib + mpltex ACS via uv" in text
    assert "<title>" not in text

text6 = output6.read_text(encoding="utf-8")
assert "BVK1" not in text6
assert text6.count("QCE") >= 3
assert text6.count("AQCE") >= 3
assert text6.count("BURP") >= 3
for label in ("(a)", "(b)", "(c)", "(d)", "(e)"):
    assert label in text6, (output6, label)
for chi_n in (20, 30, 45):
    assert f"$\\chi N={chi_n}$" in text6

text7 = output7.read_text(encoding="utf-8")
for label in ("(a)", "(b)", "(c)", "(d)", "(e)"):
    assert label not in text7, (output7, label)
for model in ("SCFT", "AQCE"):
    assert model in text7
for old_label in ("PRM", "APRM"):
    assert old_label not in text6
    assert old_label not in text7

spec = importlib.util.spec_from_file_location(
    "render_macromolecules_figure8",
    PROJECT / "scripts/render_macromolecules_figure8.py",
)
assert spec is not None and spec.loader is not None
figure8 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(figure8)

figure8.validate_fingerprints()
summary = figure8.rows("benchmark_summary.csv")
profiles = figure8.rows("benchmark_profiles.csv")
crossing = figure8.rows("ldis_crossing.csv")
crossing_summary = figure8.rows("ldis_crossing_summary.csv")
figure8.validate(summary, profiles, crossing, crossing_summary)
assert {float(row["chiN"]) for row in crossing} == {
    18.25,
    18.5,
    19.0,
    20.0,
    21.0,
    22.0,
}
assert figure8.FIGURE_MODELS == (
    "scft",
    "uneyama_doi",
    "burp_ti",
    "bvk2_fixed",
    "bvk2",
)

fixed_summary = figure8.rows("summary.csv", figure8.FIXED_DATA_DIR)
fixed_profiles = figure8.rows("profiles.csv", figure8.FIXED_DATA_DIR)
figure8.validate_fixed(fixed_summary, fixed_profiles)


def expect_fixed_rejection(test_summary=fixed_summary, test_profiles=fixed_profiles):
    try:
        figure8.validate_fixed(test_summary, test_profiles)
    except ValueError:
        return
    raise AssertionError("malformed fixed-stiffness Figure 6 input passed validation")


duplicate_fixed = copy.deepcopy(fixed_summary)
duplicate_fixed[-1] = copy.deepcopy(duplicate_fixed[0])
expect_fixed_rejection(test_summary=duplicate_fixed)

unaccepted_fixed = copy.deepcopy(fixed_summary)
unaccepted_fixed[0]["accepted"] = "false"
expect_fixed_rejection(test_summary=unaccepted_fixed)

expect_fixed_rejection(test_profiles=fixed_profiles[:-1])


def expect_rejection(
    test_summary=summary,
    test_profiles=profiles,
    test_crossing=crossing,
    test_crossing_summary=crossing_summary,
) -> None:
    try:
        figure8.validate(
            test_summary,
            test_profiles,
            test_crossing,
            test_crossing_summary,
        )
    except ValueError:
        return
    raise AssertionError("malformed Figure 6 input passed validation")


duplicate_summary = copy.deepcopy(summary)
duplicate_summary[-1] = copy.deepcopy(duplicate_summary[0])
expect_rejection(test_summary=duplicate_summary)

nonfinite_summary = copy.deepcopy(summary)
nonfinite_summary[0]["period_rg"] = "NaN"
expect_rejection(test_summary=nonfinite_summary)

unaccepted_summary = copy.deepcopy(summary)
unaccepted_summary[0]["accepted"] = "false"
expect_rejection(test_summary=unaccepted_summary)

expect_rejection(test_profiles=profiles[:-1])

duplicate_profile = copy.deepcopy(profiles)
duplicate_profile[-1] = copy.deepcopy(duplicate_profile[-2])
expect_rejection(test_profiles=duplicate_profile)

unaccepted_crossing = copy.deepcopy(crossing)
unaccepted_crossing[0]["cell_status"] = "failed"
expect_rejection(test_crossing=unaccepted_crossing)

duplicate_bifurcation = copy.deepcopy(crossing_summary)
duplicate_bifurcation[-1] = copy.deepcopy(duplicate_bifurcation[0])
expect_rejection(test_crossing_summary=duplicate_bifurcation)

with tempfile.TemporaryDirectory() as directory:
    tampered = Path(directory)
    for name in figure8.SOURCE_SHA256:
        shutil.copy2(figure8.DATA_DIR / name, tampered / name)
    with (tampered / "benchmark_summary.csv").open("a", encoding="utf-8") as handle:
        handle.write("\n")
    try:
        figure8.validate_fingerprints(tampered)
    except ValueError:
        pass
    else:
        raise AssertionError("tampered Figure 6 source passed fingerprint validation")

with tempfile.TemporaryDirectory() as directory:
    tampered_fixed = Path(directory)
    for name in figure8.FIXED_SOURCE_SHA256:
        shutil.copy2(figure8.FIXED_DATA_DIR / name, tampered_fixed / name)
    with (tampered_fixed / "summary.csv").open("a", encoding="utf-8") as handle:
        handle.write("\n")
    try:
        figure8.validate_fingerprints(fixed_data_dir=tampered_fixed)
    except ValueError:
        pass
    else:
        raise AssertionError(
            "tampered fixed-stiffness Figure 6 source passed fingerprint validation"
        )

print("ABA Figures 7 and 8, data gates, and source fingerprints passed")
