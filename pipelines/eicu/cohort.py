"""eICU v2.0 analysis-unit helpers for authorized local data."""

from __future__ import annotations

import pandas as pd


REQUIRED_COLUMNS = {"patientunitstayid", "uniquepid", "hospitalid", "unitadmittime24"}


def select_first_unit_stay(stays: pd.DataFrame) -> pd.DataFrame:
    missing = REQUIRED_COLUMNS.difference(stays.columns)
    if missing:
        raise ValueError(f"missing required columns: {sorted(missing)}")
    return stays.sort_values(["uniquepid", "unitadmittime24", "patientunitstayid"], kind="stable").drop_duplicates("uniquepid", keep="first")
