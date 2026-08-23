"""Bounded utility extrema over an onset interval."""

from __future__ import annotations

import math


def utility_at(onset: float, alert: float) -> float:
    """Frozen utility: no alert/late event -1, 0--6 h 0, 6--24 h 1, >24 h -0.25."""
    if math.isinf(alert):
        return -1.0
    gap = onset - alert
    if gap < 0:
        return -1.0
    if gap < 6:
        return 0.0
    if gap <= 24:
        return 1.0
    return -0.25


def utility_bounds(lower: float, upper: float, alert: float) -> tuple[float, float]:
    if lower > upper or not (math.isfinite(lower) and math.isfinite(upper)):
        raise ValueError("invalid onset interval")
    if math.isinf(alert):
        return -1.0, -1.0
    lo, hi = lower - alert, upper - alert
    values: list[float] = []
    if lo < 0:
        values.append(-1.0)
    if hi >= 0 and lo < 6:
        values.append(0.0)
    if max(lo, 6) <= min(hi, 24):
        values.append(1.0)
    if hi > 24:
        values.append(-0.25)
    if not values:
        raise AssertionError("utility partition did not intersect onset interval")
    return min(values), max(values)
