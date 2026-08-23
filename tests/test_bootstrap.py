import numpy as np

from onsetbounds.bootstrap import outer_confidence_set, patient_bootstrap, simultaneous_paired_outer_sets


def test_deterministic_patient_bootstrap():
    records = [1.0, 2.0, 4.0]
    stat = lambda x: np.array([np.mean(x), np.mean(x)])
    a = patient_bootstrap(records, stat, 20, 123)
    b = patient_bootstrap(records, stat, 20, 123)
    assert np.array_equal(a, b)


def test_natural_space_clipping():
    result = outer_confidence_set(.01, .99, np.array([-.2, 0]), np.array([1.2, 1.0]))
    assert result["outer_lower"] == 0
    assert result["outer_upper"] == 1


def test_simultaneous_paired_clipping():
    estimates = np.array([[-.95, .95], [0, .1]])
    boot = np.array([[[-1.2, 1.2], [-.2, .4]], [[-.9, .9], [0, .1]]])
    outer, critical = simultaneous_paired_outer_sets(estimates, boot, alpha=.5)
    assert critical >= 0
    assert outer.min() >= -1 and outer.max() <= 1
