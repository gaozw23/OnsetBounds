"""Evaluate externally frozen alerts without fitting or threshold selection."""

from __future__ import annotations

from onsetbounds.inference import aggregate_contributions, patient_contribution


def evaluate_records(records: list[dict]) -> dict:
    contributions = [
        patient_contribution(
            bool(row["event"]), bool(row.get("evaluable", row["event"])),
            row.get("lower"), row.get("upper"), float(row["alert"]),
        )
        for row in records
    ]
    return aggregate_contributions(contributions)
