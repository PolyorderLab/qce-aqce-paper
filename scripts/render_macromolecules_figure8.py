#!/usr/bin/env python3
"""Render the canonical ABA transfer figures for Macromolecules.

Figure 7 contains the profile and error comparisons, and Figure 8 contains
the LAM--DIS bifurcation comparison. Both are assembled directly from the
accepted comprehensive-validation tables. The separate cell-stress scan
remains a Supporting Information artifact.
"""

from __future__ import annotations

import csv
import hashlib
import math
from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams.update(
    {
        "svg.fonttype": "none",
        "svg.hashsalt": "macromolecules-aba-figures7-8-v2",
    }
)
import matplotlib.pyplot as plt
import mpltex
import mpltex.acs as mpltex_acs


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "results" / "bvk2_aba_comprehensive_validation"
FIXED_DATA_DIR = ROOT / "results" / "bvk2_aba_fixed_stiffness_validation"
OUTPUT = DATA_DIR / "macromolecules_figure7.svg"
OUTPUT_BIFURCATION = DATA_DIR / "macromolecules_figure8.svg"
SVG_METADATA = {"Date": None, "Creator": "Matplotlib + mpltex ACS via uv"}

SOURCE_MODELS = ("scft", "uneyama_doi", "burp_ti", "bvk1", "bvk2")
FIGURE_MODELS = ("scft", "uneyama_doi", "burp_ti", "bvk2_fixed", "bvk2")
COMPARISON_MODELS = FIGURE_MODELS[1:]
CASE_CHI_N = {
    "f0.50_chiN20": 20.0,
    "f0.50_chiN22": 22.0,
    "f0.50_chiN25": 25.0,
    "f0.50_chiN30": 30.0,
    "f0.50_chiN35": 35.0,
    "f0.50_chiN40": 40.0,
    "f0.50_chiN45": 45.0,
}
PROFILE_SAMPLES_PER_STATE = 256
BENCHMARK_SCHEMA = "bvk2-aba-comprehensive-validation-v2"
CROSSING_SCHEMA = "bvk2-aba-comprehensive-validation-v1"
FIXED_SCHEMA = "bvk2-aba-fixed-stiffness-validation-v1"
SOURCE_SHA256 = {
    "benchmark_summary.csv": "fe890e94ed44ee5a1cbe5dd8608d583e3f9f63e4066d08f2a1c213c83f636e4b",
    "benchmark_profiles.csv": "0f206c11c4664a719a5fc4ab9ad5442a9e0d47bc34703a3347be0cbdf8debf74",
    "ldis_crossing.csv": "8d73746ba9dc434fdf110fe859a4255e14cd49bf2f38ccb83f2f5fa912c30476",
    "ldis_crossing_summary.csv": "ef23212c5e9bbd33a4e9d12952b258135b060e927719c447df4533df79c345b3",
}
FIXED_SOURCE_SHA256 = {
    "summary.csv": "e32430a11351d713672aa3fd5615d7b46985554f7551c2a72509299be474aaec",
    "profiles.csv": "c58ba77ad042a2a5e208cc0c8bc5be58ff8e46e3fd1bea3b91301df08b0f5096",
}
COLORS = {
    "scft": mpltex.tableau_10[7],
    "uneyama_doi": mpltex.tableau_10[2],
    "burp_ti": mpltex.tableau_10[3],
    "bvk1": mpltex.tableau_10[9],
    "bvk2": mpltex.tableau_10[1],
    "bvk2_fixed": mpltex.tableau_10[4],
}
LABELS = {
    "scft": "SCFT",
    "uneyama_doi": "UD",
    "burp_ti": "BURP",
    "bvk1": "nonsmooth precursor",
    "bvk2": "AQCE",
    "bvk2_fixed": "QCE",
}
MARKERS = {
    "uneyama_doi": "o",
    "burp_ti": "^",
    "bvk2_fixed": "s",
    "bvk2": "D",
}


def rows(name: str, data_dir: Path = DATA_DIR) -> list[dict[str, str]]:
    with (data_dir / name).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def validate_fingerprints(
    data_dir: Path = DATA_DIR, fixed_data_dir: Path = FIXED_DATA_DIR
) -> None:
    for name, expected in SOURCE_SHA256.items():
        path = data_dir / name
        if not path.is_file():
            raise ValueError(f"Figure 6 source is missing: {path}")
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
        if actual != expected:
            raise ValueError(
                f"Figure 6 source fingerprint mismatch for {name}: "
                f"expected {expected}, found {actual}"
            )
    for name, expected in FIXED_SOURCE_SHA256.items():
        path = fixed_data_dir / name
        if not path.is_file():
            raise ValueError(f"Figure 6 fixed-stiffness source is missing: {path}")
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
        if actual != expected:
            raise ValueError(
                f"Figure 6 fixed-stiffness fingerprint mismatch for {name}: "
                f"expected {expected}, found {actual}"
            )


def validate_fixed(
    summary: list[dict[str, str]], profiles: list[dict[str, str]]
) -> None:
    expected_cases = set(CASE_CHI_N)
    if len(summary) != len(expected_cases):
        raise ValueError("Figure 6 requires seven fixed-stiffness ABA rows")
    if {row.get("schema") for row in summary} != {FIXED_SCHEMA}:
        raise ValueError("Figure 6 fixed-stiffness rows have the wrong schema")
    if {row.get("case_id") for row in summary} != expected_cases:
        raise ValueError("Figure 6 fixed-stiffness state inventory is incomplete")
    if any(row.get("model") != "bvk2_fixed" for row in summary):
        raise ValueError("Figure 6 fixed-stiffness rows have the wrong model")
    if any(row.get("accepted", "").lower() != "true" for row in summary):
        raise ValueError("Figure 6 fixed-stiffness rows contain an unaccepted state")
    numeric = (
        "fA",
        "chiN",
        "period_rg",
        "scft_period_rg",
        "period_relative_error",
        "profile_rms",
        "profile_correlation",
        "energy_density",
    )
    for row in summary:
        try:
            values = [float(row[field]) for field in numeric]
        except (KeyError, TypeError, ValueError) as error:
            raise ValueError("Figure 6 fixed-stiffness data are malformed") from error
        if not all(math.isfinite(value) for value in values):
            raise ValueError("Figure 6 fixed-stiffness data contain nonfinite values")
        if float(row["fA"]) != 0.5 or float(row["chiN"]) != CASE_CHI_N[row["case_id"]]:
            raise ValueError("Figure 6 fixed-stiffness case metadata are inconsistent")

    expected_profiles = len(expected_cases) * PROFILE_SAMPLES_PER_STATE
    if len(profiles) != expected_profiles:
        raise ValueError(
            f"Figure 6 requires exactly {expected_profiles} fixed-stiffness samples"
        )
    counts = Counter(row.get("case_id") for row in profiles)
    if set(counts) != expected_cases or any(
        count != PROFILE_SAMPLES_PER_STATE for count in counts.values()
    ):
        raise ValueError("Figure 6 fixed-stiffness profiles are incomplete")
    keys: list[tuple[str | None, float]] = []
    for row in profiles:
        if row.get("schema") != FIXED_SCHEMA or row.get("model") != "bvk2_fixed":
            raise ValueError("Figure 6 fixed-stiffness profile metadata are inconsistent")
        try:
            s = float(row["s"])
            phi_a = float(row["phi_a"])
        except (KeyError, TypeError, ValueError) as error:
            raise ValueError("Figure 6 fixed-stiffness profiles are malformed") from error
        if not all(math.isfinite(value) for value in (s, phi_a)) or not 0.0 <= s < 1.0:
            raise ValueError("Figure 6 fixed-stiffness profiles contain invalid data")
        keys.append((row.get("case_id"), s))
    if len(set(keys)) != len(keys):
        raise ValueError("Figure 6 fixed-stiffness profiles contain duplicate samples")


def panel(ax, label: str) -> None:
    ax.text(
        -0.12,
        1.03,
        label,
        transform=ax.transAxes,
        fontsize=10,
        fontweight="bold",
        va="bottom",
    )


def validate(
    summary: list[dict[str, str]],
    profiles: list[dict[str, str]],
    crossing: list[dict[str, str]],
    crossing_summary: list[dict[str, str]],
) -> None:
    expected_cases = set(CASE_CHI_N)
    expected_summary_keys = {
        (case_id, model) for case_id in expected_cases for model in SOURCE_MODELS
    }

    if len(summary) != len(expected_summary_keys):
        raise ValueError("Figure 6 requires exactly 7 x 5 benchmark rows")
    if {row.get("schema") for row in summary} != {BENCHMARK_SCHEMA}:
        raise ValueError("Figure 6 benchmark schema is not the accepted v2 schema")
    summary_keys = [(row.get("case_id"), row.get("model")) for row in summary]
    if len(set(summary_keys)) != len(summary_keys):
        raise ValueError("Figure 6 benchmark contains duplicate state/model rows")
    if set(summary_keys) != expected_summary_keys:
        raise ValueError("Figure 6 benchmark state/model grid is incomplete")
    if any(row.get("accepted", "").lower() != "true" for row in summary):
        raise ValueError("Figure 6 benchmark contains an unaccepted state")
    summary_numeric = (
        "fA",
        "chiN",
        "period_rg",
        "scft_period_rg",
        "period_relative_error",
        "profile_rms",
        "profile_correlation",
        "energy_density",
        "grid",
    )
    for row in summary:
        try:
            values = {field: float(row[field]) for field in summary_numeric}
        except (KeyError, TypeError, ValueError) as error:
            raise ValueError("Figure 6 benchmark has malformed numeric data") from error
        if not all(math.isfinite(value) for value in values.values()):
            raise ValueError("Figure 6 benchmark contains nonfinite plotted data")
        if values["fA"] != 0.5 or values["chiN"] != CASE_CHI_N[row["case_id"]]:
            raise ValueError("Figure 6 benchmark case metadata is inconsistent")
        if values["grid"] <= 0 or not values["grid"].is_integer():
            raise ValueError("Figure 6 benchmark grid must be a positive integer")

    expected_profile_rows = len(expected_summary_keys) * PROFILE_SAMPLES_PER_STATE
    if len(profiles) != expected_profile_rows:
        raise ValueError(
            f"Figure 6 requires exactly {expected_profile_rows} profile samples"
        )
    if {row.get("schema") for row in profiles} != {BENCHMARK_SCHEMA}:
        raise ValueError("Figure 6 profile schema is not the accepted v2 schema")
    profile_counts = Counter(
        (row.get("case_id"), row.get("model")) for row in profiles
    )
    if set(profile_counts) != expected_summary_keys or any(
        count != PROFILE_SAMPLES_PER_STATE for count in profile_counts.values()
    ):
        raise ValueError("Figure 6 profiles do not form a complete 7 x 5 grid")
    profile_keys: list[tuple[str | None, str | None, float]] = []
    for row in profiles:
        try:
            s = float(row["s"])
            phi_a = float(row["phi_a"])
            f_a = float(row["fA"])
            chi_n = float(row["chiN"])
        except (KeyError, TypeError, ValueError) as error:
            raise ValueError("Figure 6 profiles have malformed numeric data") from error
        if not all(math.isfinite(value) for value in (s, phi_a, f_a, chi_n)):
            raise ValueError("Figure 6 profiles contain nonfinite data")
        if not 0.0 <= s < 1.0:
            raise ValueError("Figure 6 profile coordinates must lie in [0, 1)")
        if f_a != 0.5 or chi_n != CASE_CHI_N[row["case_id"]]:
            raise ValueError("Figure 6 profile case metadata is inconsistent")
        profile_keys.append((row.get("case_id"), row.get("model"), s))
    if len(set(profile_keys)) != len(profile_keys):
        raise ValueError("Figure 6 profiles contain duplicate samples")

    crossing_models = {"scft", "bvk2"}
    crossing_chi_n = {18.25, 18.5, 19.0, 20.0, 21.0, 22.0}
    expected_crossing_keys = {
        (model, chi_n) for model in crossing_models for chi_n in crossing_chi_n
    }
    if len(crossing) != len(expected_crossing_keys):
        raise ValueError("Figure 7 L/DIS scan requires exactly twelve accepted rows")
    if {row.get("schema") for row in crossing} != {CROSSING_SCHEMA}:
        raise ValueError("Figure 7 L/DIS scan schema is not the accepted v1 schema")
    crossing_keys: list[tuple[str | None, float]] = []
    for row in crossing:
        try:
            chi_n = float(row["chiN"])
            values = (
                chi_n,
                float(row["delta_free_energy"]),
                float(row["amplitude"]),
                float(row["period_rg"]),
            )
        except (KeyError, TypeError, ValueError) as error:
            raise ValueError("Figure 7 L/DIS scan has malformed numeric data") from error
        if not all(math.isfinite(value) for value in values):
            raise ValueError("Figure 7 L/DIS scan contains nonfinite plotted data")
        if row.get("field_status") not in {"true", "Polyorder.Successful()"}:
            raise ValueError("Figure 7 L/DIS scan contains an unaccepted field solve")
        if row.get("cell_status") not in {"true", "Polyorder.Successful()"}:
            raise ValueError("Figure 7 L/DIS scan contains an unaccepted cell solve")
        crossing_keys.append((row.get("model"), chi_n))
    if len(set(crossing_keys)) != len(crossing_keys):
        raise ValueError("Figure 7 L/DIS scan contains duplicate accepted rows")
    if set(crossing_keys) != expected_crossing_keys:
        raise ValueError("Figure 7 L/DIS accepted scan grid is incomplete")

    if len(crossing_summary) != len(crossing_models):
        raise ValueError("Figure 7 requires exactly two L/DIS bifurcation rows")
    if {row.get("schema") for row in crossing_summary} != {CROSSING_SCHEMA}:
        raise ValueError("Figure 7 L/DIS summary schema is not the accepted v1 schema")
    summary_models = [row.get("model") for row in crossing_summary]
    if len(set(summary_models)) != len(summary_models) or set(summary_models) != crossing_models:
        raise ValueError("Figure 7 L/DIS bifurcations must be unique for SCFT and AQCE")
    for row in crossing_summary:
        if row.get("transition") != "L/DIS" or row.get("nonlinear_refit") != "false":
            raise ValueError("Figure 7 L/DIS summary has incompatible semantics")
        try:
            values = (
                float(row["chiN_root"]),
                float(row["period_rg_at_bifurcation"]),
            )
        except (KeyError, TypeError, ValueError) as error:
            raise ValueError("Figure 7 L/DIS summary has malformed numeric data") from error
        if not all(math.isfinite(value) for value in values):
            raise ValueError("Figure 7 L/DIS summary contains nonfinite data")


def validated_sources(data_dir: Path = DATA_DIR):
    validate_fingerprints(data_dir)
    summary = rows("benchmark_summary.csv", data_dir)
    profiles = rows("benchmark_profiles.csv", data_dir)
    crossing = rows("ldis_crossing.csv", data_dir)
    crossing_summary = rows("ldis_crossing_summary.csv", data_dir)
    validate(summary, profiles, crossing, crossing_summary)
    fixed_summary = rows("summary.csv", FIXED_DATA_DIR)
    fixed_profiles = rows("profiles.csv", FIXED_DATA_DIR)
    validate_fixed(fixed_summary, fixed_profiles)
    return (
        summary + fixed_summary,
        profiles + fixed_profiles,
        crossing,
        crossing_summary,
    )


def save_figure(fig, output: Path) -> Path:
    output.parent.mkdir(parents=True, exist_ok=True)
    output_format = output.suffix.lower().lstrip(".") or "svg"
    metadata = SVG_METADATA if output_format == "svg" else None
    fig.savefig(output, format=output_format, dpi=180, metadata=metadata)
    plt.close(fig)
    if output_format == "svg":
        output.write_text(
            "\n".join(
                line.rstrip()
                for line in output.read_text(encoding="utf-8").splitlines()
            )
            + "\n",
            encoding="utf-8",
        )
    return output


@mpltex.acs_decorator
def render_transfer(data_dir: Path = DATA_DIR, output: Path = OUTPUT) -> Path:
    summary, profiles, _, _ = validated_sources(data_dir)

    fig = plt.figure(
        figsize=(mpltex_acs.width_double_column, 5.65), constrained_layout=True
    )
    grid = fig.add_gridspec(3, 6, height_ratios=(0.10, 1.0, 1.05))
    legend_ax = fig.add_subplot(grid[0, :])
    legend_ax.axis("off")
    profile_axes = [
        fig.add_subplot(grid[1, 0:2]),
        fig.add_subplot(grid[1, 2:4]),
        fig.add_subplot(grid[1, 4:6]),
    ]
    period_ax = fig.add_subplot(grid[2, 0:3])
    rms_ax = fig.add_subplot(grid[2, 3:6])

    profile_cases = (
        ("f0.50_chiN20", 20),
        ("f0.50_chiN30", 30),
        ("f0.50_chiN45", 45),
    )
    legend_handles = None
    legend_labels = None
    for profile_ax, (profile_case, chi_n) in zip(profile_axes, profile_cases):
        for model in FIGURE_MODELS:
            subset = sorted(
                (
                    row
                    for row in profiles
                    if row["case_id"] == profile_case and row["model"] == model
                ),
                key=lambda row: float(row["s"]),
            )
            profile_ax.plot(
                [float(row["s"]) for row in subset],
                [float(row["phi_a"]) for row in subset],
                color=COLORS[model],
                label=LABELS[model],
            )
        profile_ax.set(
            xlabel=r"$x/D$",
            xlim=(0.0, 1.0),
            ylim=(-0.03, 1.03),
        )
        profile_ax.text(
            0.95,
            0.93,
            rf"$\chi N={chi_n}$",
            transform=profile_ax.transAxes,
            ha="right",
            va="top",
        )
        if legend_handles is None:
            legend_handles, legend_labels = profile_ax.get_legend_handles_labels()
    profile_axes[0].set_ylabel(r"$\phi_A$")
    for profile_ax in profile_axes[1:]:
        profile_ax.tick_params(labelleft=False)

    for model in COMPARISON_MODELS:
        subset = sorted(
            (row for row in summary if row["model"] == model),
            key=lambda row: float(row["chiN"]),
        )
        chi_n = [float(row["chiN"]) for row in subset]
        period_ax.plot(
            chi_n,
            [100.0 * float(row["period_relative_error"]) for row in subset],
            color=COLORS[model],
            marker=MARKERS[model],
            label=LABELS[model],
        )
        rms_ax.plot(
            chi_n,
            [float(row["profile_rms"]) for row in subset],
            color=COLORS[model],
            marker=MARKERS[model],
            label=LABELS[model],
        )
    period_ax.axhline(0.0, color=mpltex.almost_black, linewidth=0.7)
    period_ax.set(xlabel=r"$\chi N$", ylabel=r"period error (\%)")
    rms_ax.set(xlabel=r"$\chi N$", ylabel="profile RMS")

    for label, ax in zip(("(a)", "(b)", "(c)"), profile_axes):
        ax.text(
            0.03,
            0.97,
            label,
            transform=ax.transAxes,
            fontsize=10,
            fontweight="bold",
            va="top",
        )
    for label, ax in zip(("(d)", "(e)"), (period_ax, rms_ax)):
        panel(ax, label)
    legend_ax.legend(
        legend_handles,
        legend_labels,
        ncol=5,
        loc="center",
    )
    period_ax.legend(ncol=2)
    rms_ax.legend(ncol=2)
    return save_figure(fig, output)


@mpltex.acs_decorator
def render_bifurcation(
    data_dir: Path = DATA_DIR, output: Path = OUTPUT_BIFURCATION
) -> Path:
    _, _, crossing, crossing_summary = validated_sources(data_dir)
    fig, ax = plt.subplots(
        figsize=(mpltex_acs.width_single_column, 2.65),
        constrained_layout=True,
    )

    for model in ("scft", "bvk2"):
        subset = sorted(
            (row for row in crossing if row["model"] == model),
            key=lambda row: float(row["chiN"]),
        )
        ax.plot(
            [float(row["chiN"]) for row in subset],
            [float(row["delta_free_energy"]) for row in subset],
            color=COLORS[model],
            marker="o",
            label=LABELS[model],
        )
        root = next(
            float(row["chiN_root"])
            for row in crossing_summary
            if row["model"] == model
        )
        ax.axvline(
            root,
            color=COLORS[model],
            linestyle=(0, (2, 2)),
            linewidth=0.9,
        )
    ax.axhline(0.0, color=mpltex.almost_black, linewidth=0.7)
    ax.set(
        xlabel=r"$\chi N$",
        ylabel=r"$F_{\mathrm{LAM}}/V-F_{\mathrm{DIS}}/V$",
    )

    ax.legend()
    return save_figure(fig, output)


def render(data_dir: Path = DATA_DIR) -> tuple[Path, Path]:
    return render_transfer(data_dir), render_bifurcation(data_dir)


if __name__ == "__main__":
    for output in render():
        print(output.relative_to(ROOT))
