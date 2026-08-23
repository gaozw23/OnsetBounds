"""Frozen model-selection and first-alert metrics."""

from __future__ import annotations

import numpy as np


def weighted_average_precision(y_true: np.ndarray, scores: np.ndarray, weights: np.ndarray) -> float:
    """Patient-weighted AP with exact grouping of tied scores."""
    y = np.asarray(y_true, dtype=float)
    s = np.asarray(scores, dtype=float)
    w = np.asarray(weights, dtype=float)
    if y.ndim != 1 or not (len(y) == len(s) == len(w)):
        raise ValueError("inputs must be aligned one-dimensional arrays")
    if len(y) == 0 or not np.isfinite(s).all() or not np.isfinite(w).all() or (w < 0).any():
        raise ValueError("invalid weighted AP inputs")
    positives = float(np.dot(w, y))
    if positives <= 0:
        return float("nan")
    order = np.argsort(-s, kind="stable")
    y, s, w = y[order], s[order], w[order]
    ends = np.r_[np.flatnonzero(s[:-1] != s[1:]), len(s) - 1]
    weighted_positive = w * y
    tp = np.cumsum(weighted_positive)[ends]
    predicted = np.cumsum(w)[ends]
    recall = tp / positives
    precision = tp / predicted
    delta_recall = np.diff(np.r_[0.0, recall])
    return float(np.dot(delta_recall, precision))


def select_strict_threshold(patient_max_scores: np.ndarray, all_scores: np.ndarray, target_rate: float) -> dict[str, float | int | str]:
    """Select ``score > c`` threshold closest to target; break ties upward."""
    maxima = np.sort(np.asarray(patient_max_scores, dtype=float))
    scores = np.asarray(all_scores, dtype=float)
    if not len(maxima) or not (np.isfinite(maxima).all() and np.isfinite(scores).all()):
        raise ValueError("finite nonempty scores required")
    candidates = np.unique(np.r_[scores, -np.inf, np.inf])
    rates = (len(maxima) - np.searchsorted(maxima, candidates, side="right")) / len(maxima)
    distance = np.abs(rates - target_rate)
    minimum = float(distance.min())
    tied = np.flatnonzero(np.isclose(distance, minimum, rtol=0, atol=1e-15))
    chosen = int(tied[np.argmax(candidates[tied])])
    threshold = float(candidates[chosen])
    if not np.isfinite(threshold):
        raise ValueError("target implies an infinite threshold")
    observed = float((maxima > threshold).mean())
    return {
        "threshold": threshold,
        "observed_rate": observed,
        "absolute_target_difference": minimum,
        "candidate_thresholds": int(len(candidates)),
        "alert_rule": "score > c",
        "tie_rule": "higher threshold",
    }


def first_alert_hour(scores: np.ndarray, threshold: float, hours: np.ndarray | None = None) -> float:
    """Return the first hour with ``score > threshold`` or ``+inf``."""
    score = np.asarray(scores, dtype=float)
    time = np.arange(len(score), dtype=float) if hours is None else np.asarray(hours, dtype=float)
    if len(score) != len(time):
        raise ValueError("scores and hours must align")
    hit = np.flatnonzero(score > threshold)
    return float(time[hit[0]]) if len(hit) else float("inf")
