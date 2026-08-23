"""Small exact checks for analytic interval geometry."""

from __future__ import annotations

import math

from .bounds import single_policy_bounds, valid_at


def enumerate_valid_extrema(lower: float, upper: float, alert: float, lead_min: float = 6.0, lead_max: float = 24.0) -> tuple[int, int]:
    points = [lower, upper]
    if math.isfinite(alert):
        points.extend([max(lower, min(upper, alert + lead_min)), max(lower, min(upper, alert + lead_max))])
    points = sorted(set(points))
    probes = points + [(a + b) / 2 for a, b in zip(points[:-1], points[1:]) if a < b]
    values = [valid_at(x, alert, lead_min, lead_max) for x in probes]
    return min(values), max(values)


def analytic_matches_enumeration(lower: float, upper: float, alert: float) -> bool:
    bound = single_policy_bounds(lower, upper, alert)
    return (bound.valid_lower, bound.valid_upper) == enumerate_valid_extrema(lower, upper, alert)
