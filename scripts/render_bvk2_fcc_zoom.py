#!/usr/bin/env python3
"""Render the BVK2 FCC-region phase-diagram zoom."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams.update(
    {
        "svg.fonttype": "none",
        "svg.hashsalt": "bvk2-fcc-pocket-v2",
    }
)
import matplotlib.pyplot as plt
import mpltex
import mpltex.acs as mpltex_acs
from matplotlib.lines import Line2D

from assemble_bvk2_publication_phase_diagram import load_rpa_stability_limit


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
DEFAULT_BOUNDARIES = (
    RESULTS / "bvk2_publication_phase_diagram" / "accepted_phase_boundaries.csv"
)
DEFAULT_LEGACY_FCC = (
    RESULTS / "bvk2_fcc_o70_boundary_campaign" / "boundary_brackets.csv"
)
DEFAULT_STRICT_FCC = (
    RESULTS / "bvk2_fcc_o70_boundary_campaign"
    / "matched_dx_fcc_boundary_brackets.csv"
)
DEFAULT_STRICT_POINTS = (
    RESULTS / "bvk2_fcc_o70_boundary_campaign"
    / "matched_dx_fcc_boundary_points.csv"
)
DEFAULT_SCFT = (
    RESULTS / "bvk2_publication_phase_diagram" / "scft_reference_boundaries.csv"
)
DEFAULT_RPA_STABILITY = (
    RESULTS / "bvk2_publication_phase_diagram" / "rpa_stability_limit.csv"
)
DEFAULT_OUTPUT_PREFIX = (
    RESULTS / "bvk2_publication_phase_diagram" / "fcc_region_zoom"
)

SVG_METADATA = {"Date": None, "Creator": "Matplotlib + mpltex ACS via uv"}
X_LIMITS = (0.10, 0.22)
Y_LIMITS = (20.0, 50.0)
CONTEXT_TRANSITIONS = ("S/DIS", "S/C")
STRICT_TRANSITIONS = ("FCC/DIS", "FCC/BCC")
COLORS = {
    "S/DIS": "#1f77b4",
    "S/C": "#d62728",
    "FCC/DIS": "#2ca02c",
    "FCC/BCC": "#9467bd",
}
DISPLAY_LABELS = {
    "S/DIS": "BCC/DIS",
    "S/C": "BCC/C",
}
MARKERS = {
    "S/DIS": "s",
    "S/C": "o",
    "FCC/DIS": "^",
    "FCC/BCC": "D",
}
SCFT_COLOR = mpltex.almost_black


@dataclass(frozen=True)
class TriplePointEstimate:
    f_a: float
    chi_n: float
    f_a_uncertainty: float
    chi_n_uncertainty: float
    chi_n_low: float
    chi_n_high: float
    fcc_dis_low: float
    bcc_dis_low: float
    fcc_dis_high: float
    fcc_bcc_high: float
    width_low: float
    width_high: float
    interpolation_fraction: float


def read_rows(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def number(value: str | None) -> float | None:
    try:
        result = float(value or "")
    except ValueError:
        return None
    return result if result == result else None


def normalize_transition(value: str) -> str | None:
    if value in {"DIS/S", "DIS/S_cp", "S_cp/S"}:
        return "S/DIS"
    return value if value in CONTEXT_TRANSITIONS else None


def accepted_context(path: Path) -> dict[str, list[tuple[float, float]]]:
    """Read only non-FCC accepted boundaries used as phase-diagram context."""
    grouped: dict[str, list[tuple[float, float]]] = defaultdict(list)
    for row in read_rows(path):
        transition = normalize_transition(row.get("transition", ""))
        f_a = number(row.get("fA"))
        chi_n = number(row.get("chiN"))
        if (
            transition is None
            or row.get("status") not in {"accepted", "critical_endpoint"}
            or f_a is None
            or chi_n is None
        ):
            continue
        grouped[transition].append((f_a, chi_n))
    return {
        transition: sorted(set(points), key=lambda point: point[1])
        for transition, points in grouped.items()
    }


def strict_fcc_boundaries(
    path: Path,
) -> dict[str, list[tuple[float, float, float]]]:
    """Read accepted matched-dx FCC boundary points and their fA uncertainty."""
    grouped: dict[str, list[tuple[float, float, float]]] = defaultdict(list)
    for row in read_rows(path):
        transition = row.get("boundary", "")
        f_a = number(row.get("fA_estimate"))
        chi_n = number(row.get("chiN"))
        uncertainty = number(row.get("coordinate_uncertainty"))
        if (
            transition not in STRICT_TRANSITIONS
            or row.get("status") != "accepted"
            or row.get("equilibrium_boundary", "").lower() != "true"
            or row.get("resolution_uncertainty_pass", "").lower() != "true"
            or f_a is None
            or chi_n is None
            or uncertainty is None
        ):
            continue
        grouped[transition].append((f_a, chi_n, uncertainty))
    return {
        transition: sorted(set(points), key=lambda point: point[1])
        for transition, points in grouped.items()
    }


def linear_zero(
    rows: list[dict[str, str]],
    field: str,
    *,
    require_bracket: bool,
) -> float:
    samples = sorted(
        (
            (number(row.get("fA")), number(row.get(field)))
            for row in rows
        ),
        key=lambda sample: sample[0] if sample[0] is not None else float("inf"),
    )
    samples = [
        (f_a, delta)
        for f_a, delta in samples
        if f_a is not None and delta is not None
    ]
    if require_bracket:
        pairs = [
            (left, right)
            for left, right in zip(samples, samples[1:])
            if left[1] * right[1] <= 0.0
        ]
        if not pairs:
            raise RuntimeError(f"missing signed bracket for {field}")
        left, right = min(
            pairs,
            key=lambda pair: abs(pair[0][1]) + abs(pair[1][1]),
        )
    else:
        if len(samples) < 2:
            raise RuntimeError(f"need two samples to extrapolate {field}")
        left, right = sorted(
            sorted(samples, key=lambda sample: abs(sample[1]))[:2],
            key=lambda sample: sample[0],
        )
    if right[1] == left[1]:
        raise RuntimeError(f"degenerate linear root for {field}")
    return left[0] - left[1] * (right[0] - left[0]) / (right[1] - left[1])


def estimate_triple_point(
    points_path: Path,
    strict: dict[str, list[tuple[float, float, float]]],
) -> TriplePointEstimate:
    """Estimate FCC/BCC/DIS closure from the nearest resolved slices."""
    common_chi = sorted(
        {point[1] for point in strict["FCC/DIS"]}
        & {point[1] for point in strict["FCC/BCC"]}
    )
    if not common_chi:
        raise RuntimeError("missing common accepted FCC pocket slice")
    chi_high = common_chi[0]
    strict_at_high = {
        transition: next(point for point in strict[transition] if point[1] == chi_high)
        for transition in STRICT_TRANSITIONS
    }

    rows = read_rows(points_path)
    candidate_chi = sorted(
        {
            chi_n
            for row in rows
            if (chi_n := number(row.get("chiN"))) is not None and chi_n < chi_high
        }
    )
    if not candidate_chi:
        raise RuntimeError("missing lower matched-dx triple-point slice")
    chi_low = candidate_chi[-1]
    low_rows = [row for row in rows if number(row.get("chiN")) == chi_low]

    # At chi_low BCC/DIS is a true signed bracket. FCC/DIS lies only a few
    # 1e-6 beyond the last sample, so its two closest-to-zero values provide
    # the bounded locator extrapolation used solely for triple-point closure.
    bcc_dis_low = linear_zero(low_rows, "delta_BCC_DIS", require_bracket=True)
    fcc_dis_low = linear_zero(low_rows, "delta_FCC_DIS", require_bracket=False)
    fcc_dis_high, _, fcc_dis_uncertainty = strict_at_high["FCC/DIS"]
    fcc_bcc_high, _, fcc_bcc_uncertainty = strict_at_high["FCC/BCC"]

    width_low = bcc_dis_low - fcc_dis_low
    width_high = fcc_bcc_high - fcc_dis_high
    if not width_low < 0.0 < width_high:
        raise RuntimeError("FCC pocket width does not bracket closure")
    fraction = -width_low / (width_high - width_low)
    chi_n = chi_low + fraction * (chi_high - chi_low)
    f_a = fcc_dis_low + fraction * (fcc_dis_high - fcc_dis_low)

    coordinate_floor = max(fcc_dis_uncertainty, fcc_bcc_uncertainty)
    width_slope = abs((width_high - width_low) / (chi_high - chi_low))
    chi_n_uncertainty = 2.0 * coordinate_floor / width_slope
    fcc_dis_slope = abs((fcc_dis_high - fcc_dis_low) / (chi_high - chi_low))
    f_a_uncertainty = coordinate_floor + fcc_dis_slope * chi_n_uncertainty
    return TriplePointEstimate(
        f_a=f_a,
        chi_n=chi_n,
        f_a_uncertainty=f_a_uncertainty,
        chi_n_uncertainty=chi_n_uncertainty,
        chi_n_low=chi_low,
        chi_n_high=chi_high,
        fcc_dis_low=fcc_dis_low,
        bcc_dis_low=bcc_dis_low,
        fcc_dis_high=fcc_dis_high,
        fcc_bcc_high=fcc_bcc_high,
        width_low=width_low,
        width_high=width_high,
        interpolation_fraction=fraction,
    )


def legacy_fcc_dis_locators(
    path: Path,
    strict_chi_n: set[float],
) -> list[tuple[float, float]]:
    """Retain old FCC/DIS roots only as resolution-unmatched locator evidence."""
    points = []
    for row in read_rows(path):
        f_a = number(row.get("fA_estimate"))
        chi_n = number(row.get("chiN"))
        if (
            row.get("boundary") == "FCC/DIS"
            and row.get("status") == "accepted"
            and f_a is not None
            and chi_n is not None
            and chi_n not in strict_chi_n
        ):
            points.append((f_a, chi_n))
    return sorted(set(points), key=lambda point: point[1])


def scft_curves(path: Path) -> list[list[tuple[float, float]]]:
    grouped: dict[tuple[str, str], list[tuple[int, float, float]]] = defaultdict(list)
    for row in read_rows(path):
        transition = normalize_transition(row.get("transition", ""))
        if transition is None or row.get("side") != "left":
            continue
        f_a = number(row.get("fA"))
        chi_n = number(row.get("chiN"))
        if f_a is None or chi_n is None:
            continue
        grouped[(transition, row.get("curve_id", ""))].append(
            (int(row.get("point_index", "0")), f_a, chi_n)
        )
    return [
        [(f_a, chi_n) for _, f_a, chi_n in sorted(points)]
        for points in grouped.values()
    ]


def write_plot_data(
    path: Path,
    context: dict[str, list[tuple[float, float]]],
    strict: dict[str, list[tuple[float, float, float]]],
    locators: list[tuple[float, float]],
    triple: TriplePointEstimate,
) -> None:
    rows: list[dict[str, object]] = []
    for transition, points in context.items():
        rows.extend(
            {
                "transition": transition,
                "fA": f_a,
                "chiN": chi_n,
                "coordinate_uncertainty": "",
                "chiN_uncertainty": "",
                "display_status": "accepted_context",
            }
            for f_a, chi_n in points
            if transition != "S/DIS" or chi_n < triple.chi_n
        )
    for transition, points in strict.items():
        rows.extend(
            {
                "transition": transition,
                "fA": f_a,
                "chiN": chi_n,
                "coordinate_uncertainty": uncertainty,
                "chiN_uncertainty": "",
                "display_status": "accepted_strict_matched_dx",
            }
            for f_a, chi_n, uncertainty in points
        )
    rows.extend(
        {
            "transition": "FCC/DIS",
            "fA": f_a,
            "chiN": chi_n,
            "coordinate_uncertainty": "",
            "chiN_uncertainty": "",
            "display_status": "provisional_legacy_locator",
        }
        for f_a, chi_n in locators
    )
    rows.append(
        {
            "transition": "FCC/BCC/DIS",
            "fA": triple.f_a,
            "chiN": triple.chi_n,
            "coordinate_uncertainty": triple.f_a_uncertainty,
            "chiN_uncertainty": triple.chi_n_uncertainty,
            "display_status": "estimated_triple_point",
        }
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def write_triple_point_data(path: Path, triple: TriplePointEstimate) -> None:
    row = {
        "status": "estimated_linear_pocket_closure",
        "chiN_low": triple.chi_n_low,
        "FCC_DIS_low": triple.fcc_dis_low,
        "BCC_DIS_low": triple.bcc_dis_low,
        "signed_width_low": triple.width_low,
        "chiN_high": triple.chi_n_high,
        "FCC_DIS_high": triple.fcc_dis_high,
        "FCC_BCC_high": triple.fcc_bcc_high,
        "signed_width_high": triple.width_high,
        "interpolation_fraction": triple.interpolation_fraction,
        "chiN_estimate": triple.chi_n,
        "fA_estimate": triple.f_a,
        "chiN_uncertainty": triple.chi_n_uncertainty,
        "fA_uncertainty": triple.f_a_uncertainty,
        "derivation": (
            "linear interpolation of signed FCC pocket width from "
            "BCC/DIS-FCC/DIS at chiN_low to FCC/BCC-FCC/DIS at chiN_high"
        ),
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(row), lineterminator="\n"
        )
        writer.writeheader()
        writer.writerow(row)


@mpltex.acs_decorator
def render(
    boundaries: Path,
    legacy_fcc: Path,
    strict_fcc: Path,
    strict_points: Path,
    scft: Path,
    output_prefix: Path,
    rpa_stability: Path = DEFAULT_RPA_STABILITY,
) -> None:
    context = accepted_context(boundaries)
    strict = strict_fcc_boundaries(strict_fcc)
    if not all(strict.get(transition) for transition in STRICT_TRANSITIONS):
        raise RuntimeError("missing accepted strict matched-dx FCC boundaries")
    triple = estimate_triple_point(strict_points, strict)
    strict_fcc_dis_chi = {point[1] for point in strict["FCC/DIS"]}
    locators = [
        point
        for point in legacy_fcc_dis_locators(legacy_fcc, strict_fcc_dis_chi)
        if point[1] > triple.chi_n
    ]

    fig, ax = plt.subplots(
        figsize=(mpltex_acs.width_double_column, 3.8),
        constrained_layout=True,
    )

    rpa_curve = [
        point for point in load_rpa_stability_limit(rpa_stability)
        if X_LIMITS[0] <= point[0] <= X_LIMITS[1]
        and Y_LIMITS[0] <= point[1] <= Y_LIMITS[1]
    ]
    ax.plot(
        [point[0] for point in rpa_curve],
        [point[1] for point in rpa_curve],
        color=SCFT_COLOR,
        linewidth=1.0,
        linestyle=":",
        alpha=0.75,
        zorder=2,
    )

    for curve in scft_curves(scft):
        ax.plot(
            [point[0] for point in curve],
            [point[1] for point in curve],
            color=SCFT_COLOR,
            linewidth=0.8,
            linestyle=(0, (4, 3)),
            alpha=0.75,
            zorder=1,
        )

    for transition in CONTEXT_TRANSITIONS:
        points = context.get(transition, [])
        if transition == "S/DIS":
            stable_points = [point for point in points if point[1] < triple.chi_n]
            stable_points.append((triple.f_a, triple.chi_n))
            points = stable_points
        ax.plot(
            [point[0] for point in points],
            [point[1] for point in points],
            color=COLORS[transition],
            marker=MARKERS[transition],
            markerfacecolor="white",
            markeredgewidth=0.9,
            markersize=4.2,
            linewidth=1.2,
            zorder=3,
        )

    locator_curve = sorted(locators, key=lambda point: point[1])
    ax.plot(
        [point[0] for point in locator_curve],
        [point[1] for point in locator_curve],
        color=COLORS["FCC/DIS"],
        marker=MARKERS["FCC/DIS"],
        markerfacecolor="white",
        markeredgewidth=0.9,
        markersize=4.2,
        linewidth=1.0,
        linestyle=(0, (2, 2)),
        alpha=0.8,
        zorder=2,
    )

    for transition in STRICT_TRANSITIONS:
        points = strict[transition]
        connected = [(triple.f_a, triple.chi_n)] + [
            (f_a, chi_n) for f_a, chi_n, _ in points
        ]
        ax.plot(
            [point[0] for point in connected],
            [point[1] for point in connected],
            color=COLORS[transition],
            linewidth=1.2,
            zorder=4,
        )
        ax.errorbar(
            [point[0] for point in points],
            [point[1] for point in points],
            xerr=[point[2] for point in points],
            color=COLORS[transition],
            marker=MARKERS[transition],
            markerfacecolor="white",
            markeredgewidth=1.1,
            markersize=5.2,
            linewidth=0.0,
            elinewidth=1.0,
            capsize=2.2,
            zorder=5,
        )

    common_chi = sorted(
        {point[1] for point in strict["FCC/DIS"]}
        & {point[1] for point in strict["FCC/BCC"]}
    )
    if common_chi:
        fcc_dis_by_chi = {point[1]: point[0] for point in strict["FCC/DIS"]}
        fcc_bcc_by_chi = {point[1]: point[0] for point in strict["FCC/BCC"]}
        fill_chi = [triple.chi_n] + common_chi
        fill_dis = [triple.f_a] + [fcc_dis_by_chi[chi_n] for chi_n in common_chi]
        fill_bcc = [triple.f_a] + [fcc_bcc_by_chi[chi_n] for chi_n in common_chi]
        ax.fill_betweenx(
            fill_chi,
            fill_dis,
            fill_bcc,
            color="#f1c40f",
            alpha=0.18,
            linewidth=0.0,
            zorder=2,
        )
        ax.annotate(
            "FCC",
            xy=(
                0.5 * (fcc_dis_by_chi[common_chi[-1]] + fcc_bcc_by_chi[common_chi[-1]]),
                common_chi[-1],
            ),
            xytext=(0.171, 32.2),
            color=COLORS["FCC/BCC"],
            fontsize=8,
            ha="center",
            arrowprops={"arrowstyle": "->", "linewidth": 0.7},
            zorder=6,
        )

    ax.errorbar(
        [triple.f_a],
        [triple.chi_n],
        xerr=[triple.f_a_uncertainty],
        yerr=[triple.chi_n_uncertainty],
        color=SCFT_COLOR,
        marker="*",
        markerfacecolor="white",
        markeredgewidth=1.0,
        markersize=8.5,
        linewidth=0.0,
        elinewidth=0.8,
        capsize=2.0,
        zorder=7,
    )
    ax.annotate(
        "estimated\ntriple point",
        xy=(triple.f_a, triple.chi_n),
        xytext=(triple.f_a - 0.012, triple.chi_n + 1.6),
        fontsize=7,
        ha="center",
        arrowprops={"arrowstyle": "->", "linewidth": 0.7},
        zorder=7,
    )

    ax.set_xlim(*X_LIMITS)
    ax.set_ylim(*Y_LIMITS)
    ax.set_xticks([0.10, 0.12, 0.14, 0.16, 0.18, 0.20, 0.22])
    ax.set_yticks([20, 25, 30, 35, 40, 45, 50])
    ax.set_xlabel(r"$f_A$")
    ax.set_ylabel(r"$\chi N$")
    ax.grid(False)
    for spine in ax.spines.values():
        spine.set_visible(True)

    handles = [
        Line2D(
            [0], [0], color=COLORS[transition], marker=MARKERS[transition],
            markerfacecolor="white", markeredgewidth=0.9, markersize=4.2,
            linewidth=1.2, label=DISPLAY_LABELS[transition],
        )
        for transition in CONTEXT_TRANSITIONS
    ]
    handles.extend(
        Line2D(
            [0], [0], color=COLORS[transition], marker=MARKERS[transition],
            markerfacecolor="white", markeredgewidth=1.1, markersize=5.0,
            linewidth=1.2, label=f"{transition} (strict)",
        )
        for transition in STRICT_TRANSITIONS
    )
    handles.extend(
        [
            Line2D(
                [0], [0], color=SCFT_COLOR, marker="*",
                markerfacecolor="white", linewidth=0.0, markersize=8.0,
                label="estimated FCC/BCC/DIS triple point",
            ),
            Line2D(
                [0], [0], color=COLORS["FCC/DIS"],
                marker=MARKERS["FCC/DIS"], markerfacecolor="white",
                linewidth=1.0, linestyle=(0, (2, 2)),
                label="legacy FCC/DIS locator",
            ),
            Line2D(
                [0], [0], color=SCFT_COLOR, linewidth=1.0,
                linestyle=":", label="RPA stability limit",
            ),
            Line2D(
                [0], [0], color=SCFT_COLOR, linewidth=0.8,
                linestyle=(0, (4, 3)), label="SCFT",
            ),
        ]
    )
    ax.legend(
        handles=handles,
        loc="center left",
        bbox_to_anchor=(1.015, 0.5),
        frameon=False,
        borderaxespad=0.0,
    )

    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    write_plot_data(
        output_prefix.with_name(f"{output_prefix.name}_data.csv"),
        context, strict, locators, triple,
    )
    write_triple_point_data(
        output_prefix.with_name(f"{output_prefix.name}_triple_point.csv"),
        triple,
    )
    svg_path = output_prefix.with_suffix(".svg")
    fig.savefig(
        svg_path,
        format="svg",
        metadata=SVG_METADATA,
    )
    svg_path.write_text(
        "\n".join(
            line.rstrip()
            for line in svg_path.read_text(encoding="utf-8").splitlines()
        )
        + "\n",
        encoding="utf-8",
    )
    fig.savefig(
        output_prefix.with_suffix(".png"),
        format="png",
        dpi=240,
        metadata={"Software": "Matplotlib + mpltex ACS via uv"},
    )
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--boundaries", type=Path, default=DEFAULT_BOUNDARIES)
    parser.add_argument("--legacy-fcc", type=Path, default=DEFAULT_LEGACY_FCC)
    parser.add_argument("--strict-fcc", type=Path, default=DEFAULT_STRICT_FCC)
    parser.add_argument("--strict-points", type=Path, default=DEFAULT_STRICT_POINTS)
    parser.add_argument("--scft", type=Path, default=DEFAULT_SCFT)
    parser.add_argument(
        "--rpa-stability", type=Path, default=DEFAULT_RPA_STABILITY
    )
    parser.add_argument("--output-prefix", type=Path, default=DEFAULT_OUTPUT_PREFIX)
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    render(
        args.boundaries.resolve(),
        args.legacy_fcc.resolve(),
        args.strict_fcc.resolve(),
        args.strict_points.resolve(),
        args.scft.resolve(),
        args.output_prefix.resolve(),
        args.rpa_stability.resolve(),
    )
