"""Public feature helpers; no frozen patient-derived constants are distributed."""

from __future__ import annotations

import pandas as pd


def hourly_last_value(events: pd.DataFrame, value_column: str, origin_column: str = "intime") -> pd.DataFrame:
    frame = events.copy()
    frame["charttime"] = pd.to_datetime(frame["charttime"], errors="raise")
    frame[origin_column] = pd.to_datetime(frame[origin_column], errors="raise")
    frame["hour_index"] = ((frame["charttime"] - frame[origin_column]).dt.total_seconds() // 3600).astype(int)
    return frame.sort_values(["stay_id", "charttime"], kind="stable").groupby(["stay_id", "hour_index"], as_index=False).last()[["stay_id", "hour_index", value_column]]


def apply_external_train_medians(features: pd.DataFrame, medians: dict[str, float]) -> pd.DataFrame:
    """Apply user-supplied TRAIN-only medians; this repository includes no fitted arrays."""
    result = features.copy()
    for column, value in medians.items():
        if column in result:
            result[column] = result[column].fillna(float(value))
    return result
