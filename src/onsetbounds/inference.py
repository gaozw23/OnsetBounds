"""Patient contributions and endpoint aggregation."""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Iterable

import numpy as np

from .bounds import single_policy_bounds
from .utility import utility_bounds


@dataclass(frozen=True)
class PatientContribution:
    event: int
    evaluable: int
    alerted: int
    valid_lower: float
    valid_upper: float
    premature_lower: float
    premature_upper: float
    late_missed_lower: float
    late_missed_upper: float
    ppv_numerator_lower: float
    ppv_numerator_upper: float
    utility_lower: float
    utility_upper: float


def patient_contribution(
    event: bool,
    evaluable: bool,
    lower: float | None,
    upper: float | None,
    alert: float,
) -> PatientContribution:
    alerted = int(math.isfinite(alert))
    if event:
        if lower is None or upper is None:
            raise ValueError("event patients require an onset interval")
        bound = single_policy_bounds(lower, upper, alert)
        q_lower, q_upper = utility_bounds(lower, upper, alert)
        return PatientContribution(
            1,
            int(evaluable),
            alerted,
            bound.valid_lower,
            bound.valid_upper,
            bound.premature_lower,
            bound.premature_upper,
            bound.late_missed_lower,
            bound.late_missed_upper,
            bound.valid_lower,
            bound.valid_upper,
            q_lower,
            q_upper,
        )
    q = -float(alerted)
    return PatientContribution(0, 0, alerted, 0, 0, 0, 0, 0, 0, 0, 0, q, q)


def aggregate_contributions(contributions: Iterable[PatientContribution]) -> dict[str, float | int | str | None]:
    rows = list(contributions)
    if not rows:
        raise ValueError("at least one patient contribution is required")
    evaluable = [x for x in rows if x.evaluable]
    alert_denominator = sum(x.alerted for x in rows)
    result: dict[str, float | int | str | None] = {
        "patients": len(rows),
        "primary_evaluable_events": len(evaluable),
        "alerts": alert_denominator,
        "status": "DEFINED",
    }
    if not evaluable:
        result.update({k: None for k in (
            "valid_lower", "valid_upper", "premature_lower", "premature_upper",
            "late_missed_lower", "late_missed_upper",
        )})
        result["status"] = "MISSING_EVALUABLE_DENOMINATOR"
    else:
        for stem in ("valid", "premature", "late_missed"):
            result[f"{stem}_lower"] = float(np.mean([getattr(x, f"{stem}_lower") for x in evaluable]))
            result[f"{stem}_upper"] = float(np.mean([getattr(x, f"{stem}_upper") for x in evaluable]))
    if alert_denominator:
        result["ppv_lower"] = float(sum(x.ppv_numerator_lower for x in rows) / alert_denominator)
        result["ppv_upper"] = float(sum(x.ppv_numerator_upper for x in rows) / alert_denominator)
    else:
        result["ppv_lower"] = result["ppv_upper"] = None
        result["status"] = "MISSING_ALERT_DENOMINATOR" if evaluable else "MISSING_BOTH_DENOMINATORS"
    result["utility_lower"] = float(np.mean([x.utility_lower for x in rows]))
    result["utility_upper"] = float(np.mean([x.utility_upper for x in rows]))
    return result
