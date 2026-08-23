"""Direct shared-onset and marginal-subtraction comparisons."""

from __future__ import annotations

from dataclasses import dataclass
import math

from .bounds import single_policy_bounds, valid_at


@dataclass(frozen=True)
class PairedBounds:
    lower: int
    upper: int


def _breakpoints(lower: float, upper: float, alert_a: float, alert_b: float, lead_min: float, lead_max: float) -> list[float]:
    points = [lower, upper]
    for alert in (alert_a, alert_b):
        if math.isfinite(alert):
            points.extend([alert + lead_min, alert + lead_max])
    return sorted({min(upper, max(lower, x)) for x in points})


def paired_valid_bounds(
    lower: float,
    upper: float,
    alert_a: float,
    alert_b: float,
    lead_min: float = 6.0,
    lead_max: float = 24.0,
) -> PairedBounds:
    """Sharp bounds for ``valid_a(T) - valid_b(T)`` using the same latent ``T``."""
    single_policy_bounds(lower, upper, alert_a, lead_min, lead_max)
    single_policy_bounds(lower, upper, alert_b, lead_min, lead_max)
    points = _breakpoints(lower, upper, alert_a, alert_b, lead_min, lead_max)
    probes = list(points)
    probes.extend((a + b) / 2 for a, b in zip(points[:-1], points[1:]) if a < b)
    values = [
        valid_at(t, alert_a, lead_min, lead_max) - valid_at(t, alert_b, lead_min, lead_max)
        for t in probes
    ]
    return PairedBounds(min(values), max(values))


def marginal_subtraction(lower_a: float, upper_a: float, lower_b: float, upper_b: float) -> tuple[float, float]:
    """Difference interval obtained by subtracting separate marginal bounds."""
    return lower_a - upper_b, upper_a - lower_b
