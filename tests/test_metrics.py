import numpy as np

from onsetbounds.metrics import first_alert_hour, select_strict_threshold, weighted_average_precision


def test_weighted_ap_tie_grouping():
    value = weighted_average_precision(np.array([1, 0, 1]), np.array([.9, .9, .5]), np.ones(3))
    assert np.isclose(value, 7 / 12)


def test_strict_threshold_and_first_alert():
    selected = select_strict_threshold(np.array([.1, .2, .3]), np.array([.1, .2, .3]), 1 / 3)
    assert selected["threshold"] == .2
    assert selected["observed_rate"] == 1 / 3
    assert first_alert_hour(np.array([.2, .21]), .2) == 1
