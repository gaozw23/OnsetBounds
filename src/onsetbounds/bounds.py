"""Sharp single-policy bounds under interval-uncertain event onset."""

from __future__ import annotations

from dataclasses import dataclass
import math


@dataclass(frozen=True)
class BinaryBounds:
    valid_lower: int
    valid_upper: int
    premature_lower: int
    premature_upper: int
    late_missed_lower: int
    late_missed_upper: int


def _validate_interval(lower: float, upper: float) -> None:
    if not math.isfinite(lower) or not math.isfinite(upper) or lower > upper:
        raise ValueError("onset interval must be finite and satisfy lower <= upper")


def valid_at(onset: float, alert: float, lead_min: float, lead_max: float) -> int:
    """Closed-boundary valid indicator: alert + lead_min <= onset <= alert + lead_max."""
    if math.isinf(alert):
        return 0
    return int(alert + lead_min <= onset <= alert + lead_max)


def single_policy_bounds(
    lower: float,
    upper: float,
    alert: float,
    lead_min: float = 6.0,
    lead_max: float = 24.0,
) -> BinaryBounds:
    """Return sharp binary-event bounds over all onsets in ``[lower, upper]``.

    ``alert=+inf`` represents no alert. Boundary equality is valid, matching the
    frozen manuscript and real-data implementations.
    """
    _validate_interval(lower, upper)
    if lead_min > lead_max:
        raise ValueError("lead_min must not exceed lead_max")
    if math.isinf(alert):
        if alert < 0:
            raise ValueError("only +inf is a valid no-alert sentinel")
        return BinaryBounds(0, 0, 0, 0, 1, 1)
    if not math.isfinite(alert):
        raise ValueError("alert must be finite or +inf")
    window_left = alert + lead_min
    window_right = alert + lead_max
    valid_lower = int(lower >= window_left and upper <= window_right)
    valid_upper = int(max(lower, window_left) <= min(upper, window_right))
    premature_lower = int(lower > window_right)
    premature_upper = int(upper > window_right)
    late_lower = int(upper < window_left)
    late_upper = int(lower < window_left)
    return BinaryBounds(
        valid_lower,
        valid_upper,
        premature_lower,
        premature_upper,
        late_lower,
        late_upper,
    )
