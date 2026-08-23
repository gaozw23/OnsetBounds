"""Frozen patient and two-stage hospital/patient bootstrap utilities."""

from __future__ import annotations

from collections.abc import Callable, Sequence

import numpy as np


def patient_bootstrap(
    records: Sequence[object],
    statistic: Callable[[list[object]], np.ndarray],
    b: int,
    seed: int,
) -> np.ndarray:
    """Resample complete patient records with a deterministic generator stream."""
    n = len(records)
    if n == 0 or b <= 0:
        raise ValueError("nonempty records and positive b required")
    rng = np.random.default_rng(seed)
    outputs = []
    for _ in range(b):
        indices = rng.integers(0, n, size=n)
        outputs.append(np.asarray(statistic([records[i] for i in indices]), dtype=float))
    return np.vstack(outputs)


def two_stage_hospital_patient_bootstrap(
    hospital_records: dict[object, Sequence[object]],
    statistic: Callable[[list[object]], np.ndarray],
    b: int,
    seed: int,
) -> np.ndarray:
    """Resample hospitals, then patients within every sampled hospital cluster."""
    hospitals = list(hospital_records)
    if not hospitals or b <= 0 or any(not hospital_records[h] for h in hospitals):
        raise ValueError("nonempty hospital clusters and positive b required")
    rng = np.random.default_rng(seed)
    outputs = []
    for _ in range(b):
        sampled_hospitals = rng.integers(0, len(hospitals), size=len(hospitals))
        sample: list[object] = []
        for position in sampled_hospitals:
            members = hospital_records[hospitals[int(position)]]
            indices = rng.integers(0, len(members), size=len(members))
            sample.extend(members[int(i)] for i in indices)
        outputs.append(np.asarray(statistic(sample), dtype=float))
    return np.vstack(outputs)


def outer_confidence_set(
    estimate_lower: float,
    estimate_upper: float,
    bootstrap_lower: np.ndarray,
    bootstrap_upper: np.ndarray,
    alpha: float = 0.05,
    natural_space: tuple[float, float] = (0.0, 1.0),
) -> dict[str, float]:
    """Unstudentized outward expansion used for the frozen final analyses."""
    blo = np.asarray(bootstrap_lower, dtype=float)
    bhi = np.asarray(bootstrap_upper, dtype=float)
    if not len(blo) or len(blo) != len(bhi):
        raise ValueError("aligned bootstrap endpoints required")
    error = np.maximum(estimate_lower - blo, bhi - estimate_upper)
    critical = max(0.0, float(np.quantile(error, 1 - alpha)))
    lo, hi = natural_space
    return {
        "estimate_lower": estimate_lower,
        "estimate_upper": estimate_upper,
        "outer_lower": max(lo, estimate_lower - critical),
        "outer_upper": min(hi, estimate_upper + critical),
        "critical": critical,
    }


def simultaneous_paired_outer_sets(
    estimates: np.ndarray,
    bootstrap_endpoints: np.ndarray,
    alpha: float = 0.05,
) -> tuple[np.ndarray, float]:
    """Shared max-statistic expansion across paired comparisons."""
    estimate = np.asarray(estimates, dtype=float)
    boot = np.asarray(bootstrap_endpoints, dtype=float)
    if estimate.ndim != 2 or estimate.shape[1] != 2 or boot.ndim != 3 or boot.shape[1:] != estimate.shape:
        raise ValueError("expected estimates Kx2 and bootstrap BxKx2")
    errors = np.maximum(estimate[None, :, 0] - boot[:, :, 0], boot[:, :, 1] - estimate[None, :, 1])
    critical = max(0.0, float(np.quantile(errors.max(axis=1), 1 - alpha)))
    outer = np.column_stack((np.maximum(-1, estimate[:, 0] - critical), np.minimum(1, estimate[:, 1] + critical)))
    return outer, critical
