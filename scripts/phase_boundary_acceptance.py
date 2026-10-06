#!/usr/bin/env python3
"""Canonical acceptance contract for phase-boundary campaigns.

Acceptance has two explicit scopes:

``phase_endpoint``
    One phase at one thermodynamic coordinate.  Acceptance certifies phase
    identity, composition, cell stress, and free-energy stationarity only.

``boundary_bracket``
    A complete two-phase boundary locator.  Acceptance additionally requires
    accepted endpoints for both competitors on both sides, a signed
    free-energy bracket, and the requested coordinate uncertainty.

An accepted endpoint is never, by itself, an accepted boundary point.
"""

from __future__ import annotations

from dataclasses import dataclass
import math


ACCEPTED = "accepted"
PROVISIONAL = "provisional"
REJECTED = "rejected"

PHASE_ENDPOINT = "phase_endpoint"
BOUNDARY_BRACKET = "boundary_bracket"


@dataclass(frozen=True)
class AcceptanceContract:
    """Numerical gates shared by all current BVK2 boundary workflows."""

    mean_error_max: float = 1.0e-8
    cell_stress_target: float = 1.0e-3
    cell_stress_max: float = 2.0e-3
    energy_change_max: float = 1.0e-6
    coordinate_uncertainty_max: float = 1.0e-3

    @property
    def bracket_width_max(self) -> float:
        return 2.0 * self.coordinate_uncertainty_max


DEFAULT_CONTRACT = AcceptanceContract()


@dataclass(frozen=True)
class AcceptanceDecision:
    status: str
    claim_kind: str
    reason: str
    stress_tier: str = ""

    @property
    def accepted(self) -> bool:
        return self.status == ACCEPTED


def _finite(value: float) -> bool:
    return math.isfinite(float(value))


def classify_phase_endpoint(
    *,
    identity_ok: bool,
    mean_error: float,
    cell_stress: float,
    energy_change: float,
    audit_energy: float,
    contract: AcceptanceContract = DEFAULT_CONTRACT,
) -> AcceptanceDecision:
    """Classify one phase endpoint under the common physical contract."""
    values = (mean_error, cell_stress, energy_change, audit_energy)
    if not all(_finite(value) for value in values):
        return AcceptanceDecision(
            PROVISIONAL, PHASE_ENDPOINT, "nonfinite_endpoint_evidence"
        )
    if not identity_ok:
        return AcceptanceDecision(
            REJECTED, PHASE_ENDPOINT, "phase_identity_failed"
        )
    if abs(mean_error) > contract.mean_error_max:
        return AcceptanceDecision(
            PROVISIONAL, PHASE_ENDPOINT, "composition_gate_failed"
        )
    stress = abs(cell_stress)
    if stress > contract.cell_stress_max:
        return AcceptanceDecision(
            PROVISIONAL, PHASE_ENDPOINT, "cell_stress_gate_failed"
        )
    if abs(energy_change) > contract.energy_change_max:
        return AcceptanceDecision(
            PROVISIONAL, PHASE_ENDPOINT, "energy_change_gate_failed"
        )
    tier = (
        "target"
        if stress <= contract.cell_stress_target
        else "fallback"
    )
    return AcceptanceDecision(
        ACCEPTED,
        PHASE_ENDPOINT,
        "physical_endpoint_gates_pass",
        stress_tier=tier,
    )


def classify_boundary_bracket(
    *,
    lower_endpoints_accepted: bool,
    upper_endpoints_accepted: bool,
    coordinate_lower: float,
    coordinate_upper: float,
    delta_lower: float,
    delta_upper: float,
    contract: AcceptanceContract = DEFAULT_CONTRACT,
) -> tuple[AcceptanceDecision, float]:
    """Classify a two-phase bracket and return its secant coordinate."""
    values = (
        coordinate_lower,
        coordinate_upper,
        delta_lower,
        delta_upper,
    )
    if not all(_finite(value) for value in values):
        return (
            AcceptanceDecision(
                PROVISIONAL, BOUNDARY_BRACKET, "nonfinite_boundary_evidence"
            ),
            math.nan,
        )
    lower = float(coordinate_lower)
    upper = float(coordinate_upper)
    d_lower = float(delta_lower)
    d_upper = float(delta_upper)
    if upper < lower:
        lower, upper = upper, lower
        d_lower, d_upper = d_upper, d_lower
    if not lower_endpoints_accepted or not upper_endpoints_accepted:
        return (
            AcceptanceDecision(
                PROVISIONAL,
                BOUNDARY_BRACKET,
                "phase_endpoint_gate_failed",
            ),
            math.nan,
        )
    if d_lower * d_upper > 0.0:
        return (
            AcceptanceDecision(
                PROVISIONAL,
                BOUNDARY_BRACKET,
                "signed_bracket_missing",
            ),
            math.nan,
        )
    width = upper - lower
    if width > contract.bracket_width_max + 1.0e-12:
        return (
            AcceptanceDecision(
                PROVISIONAL,
                BOUNDARY_BRACKET,
                "coordinate_uncertainty_too_large",
            ),
            _secant(lower, upper, d_lower, d_upper),
        )
    root = _secant(lower, upper, d_lower, d_upper)
    if not lower - 1.0e-12 <= root <= upper + 1.0e-12:
        return (
            AcceptanceDecision(
                REJECTED,
                BOUNDARY_BRACKET,
                "secant_outside_signed_bracket",
            ),
            root,
        )
    return (
        AcceptanceDecision(
            ACCEPTED,
            BOUNDARY_BRACKET,
            "accepted_endpoints_signed_bracket_within_uncertainty",
        ),
        root,
    )


def _secant(lower: float, upper: float, d_lower: float, d_upper: float) -> float:
    denominator = d_upper - d_lower
    if denominator == 0.0:
        return lower if d_lower == 0.0 else math.nan
    return lower - d_lower * (upper - lower) / denominator


def require_boundary_claim(row: dict[str, object]) -> None:
    """Reject a canonical ledger row that is not a boundary-level claim."""
    if row.get("claim_kind") != BOUNDARY_BRACKET:
        raise ValueError(
            f"claim_kind={row.get('claim_kind')!r}; expected "
            f"{BOUNDARY_BRACKET!r}"
        )
    if row.get("status") not in {ACCEPTED, PROVISIONAL, REJECTED}:
        raise ValueError(f"unsupported canonical status={row.get('status')!r}")
