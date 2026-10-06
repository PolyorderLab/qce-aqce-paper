"""Endpoint-preserving display smoothing for quantized SCFT reference paths."""

from __future__ import annotations


def _monotone(values: list[float]) -> bool:
    differences = [right - left for left, right in zip(values, values[1:])]
    return all(value >= 0.0 for value in differences) or all(
        value <= 0.0 for value in differences
    )


def _pchip_slopes(values: list[float]) -> list[float]:
    """Fritsch--Carlson slopes for equally spaced monotone samples."""
    if len(values) == 2:
        slope = values[1] - values[0]
        return [slope, slope]

    secants = [right - left for left, right in zip(values, values[1:])]
    slopes = [0.0] * len(values)
    for index in range(1, len(values) - 1):
        before, after = secants[index - 1], secants[index]
        if before * after > 0.0:
            slopes[index] = 2.0 * before * after / (before + after)

    def endpoint(first: float, second: float) -> float:
        candidate = 0.5 * (3.0 * first - second)
        if candidate * first <= 0.0:
            return 0.0
        if first * second < 0.0 and abs(candidate) > 3.0 * abs(first):
            return 3.0 * first
        return candidate

    slopes[0] = endpoint(secants[0], secants[1])
    slopes[-1] = endpoint(secants[-1], secants[-2])
    return slopes


def _dequantize_monotone_values(
    values: list[float], passes: int = 10
) -> list[float]:
    """Remove one-pixel stair steps with a positive binomial display filter."""
    filtered = list(values)
    for _ in range(passes):
        filtered = [
            filtered[0],
            *(
                0.25 * filtered[index - 1]
                + 0.5 * filtered[index]
                + 0.25 * filtered[index + 1]
                for index in range(1, len(filtered) - 1)
            ),
            filtered[-1],
        ]
    return filtered


def smooth_scft_display_curve(
    points: list[tuple[float, float]], samples_per_segment: int = 5,
    *, smooth_chi: bool = True,
) -> list[tuple[float, float]]:
    """Densify a monotone digitized path without changing its data contract.

    The SCFT CSV retains the original quantized vertices. A short positive
    binomial filter removes sub-pixel stair steps, then PCHIP densifies the
    display path. Both endpoints stay exact and monotonicity is preserved.
    Non-monotone paths are returned unchanged rather than risking a topology
    change. Set smooth_chi=False to retain the raw segregation coordinates
    at each display anchor (useful for steep GYR/LAM branches).
    For the O70 curves, the maximum filter displacement is about one
    digitization increment in composition and below one increment in chiN.
    """
    if len(points) < 3 or samples_per_segment < 2:
        return list(points)
    f_values = [point[0] for point in points]
    chi_values = [point[1] for point in points]
    if not (_monotone(f_values) and _monotone(chi_values)):
        return list(points)

    display_f = _dequantize_monotone_values(f_values)
    display_chi = (
        _dequantize_monotone_values(chi_values) if smooth_chi else chi_values
    )
    f_slopes = _pchip_slopes(display_f)
    chi_slopes = _pchip_slopes(display_chi)
    result: list[tuple[float, float]] = []
    for index in range(len(points) - 1):
        for sample in range(samples_per_segment):
            t = sample / samples_per_segment
            t2, t3 = t * t, t * t * t
            h00 = 2.0 * t3 - 3.0 * t2 + 1.0
            h10 = t3 - 2.0 * t2 + t
            h01 = -2.0 * t3 + 3.0 * t2
            h11 = t3 - t2
            f_a = (
                h00 * display_f[index]
                + h10 * f_slopes[index]
                + h01 * display_f[index + 1]
                + h11 * f_slopes[index + 1]
            )
            chi_n = (
                h00 * display_chi[index]
                + h10 * chi_slopes[index]
                + h01 * display_chi[index + 1]
                + h11 * chi_slopes[index + 1]
            )
            result.append((
                min(max(f_a, min(display_f[index:index + 2])),
                    max(display_f[index:index + 2])),
                min(max(chi_n, min(display_chi[index:index + 2])),
                    max(display_chi[index:index + 2])),
            ))
    result.append(points[-1])
    return result
