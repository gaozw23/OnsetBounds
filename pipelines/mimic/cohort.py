"""MIMIC-IV cohort interface; requires an authorized local database extract."""

from __future__ import annotations

import pandas as pd


REQUIRED_COLUMNS = {"subject_id", "hadm_id", "stay_id", "intime", "outtime"}


def select_analysis_stays(stays: pd.DataFrame) -> pd.DataFrame:
    """Deterministically select the first adult ICU stay per hospitalization."""
    missing = REQUIRED_COLUMNS.difference(stays.columns)
    if missing:
        raise ValueError(f"missing required columns: {sorted(missing)}")
    frame = stays.copy()
    frame["intime"] = pd.to_datetime(frame["intime"], errors="raise")
    return frame.sort_values(["subject_id", "hadm_id", "intime", "stay_id"], kind="stable").drop_duplicates("hadm_id", keep="first")
