import math

from onsetbounds.bounds import single_policy_bounds


def test_no_alert_case():
    b = single_policy_bounds(10, 12, math.inf)
    assert (b.valid_lower, b.valid_upper, b.premature_lower, b.premature_upper, b.late_missed_lower, b.late_missed_upper) == (0, 0, 0, 0, 1, 1)


def test_certain_valid_and_boundary_equality():
    b = single_policy_bounds(6, 24, 0)
    assert (b.valid_lower, b.valid_upper) == (1, 1)


def test_partial_and_no_intersection():
    assert single_policy_bounds(2, 8, 0).valid_lower == 0
    assert single_policy_bounds(2, 8, 0).valid_upper == 1
    assert (single_policy_bounds(25, 30, 0).valid_lower, single_policy_bounds(25, 30, 0).valid_upper) == (0, 0)


def test_refinement_monotonicity():
    coarse = single_policy_bounds(0, 10, -2, 4, 10)
    refined = single_policy_bounds(3, 7, -2, 4, 10)
    assert coarse.valid_lower <= refined.valid_lower <= refined.valid_upper <= coarse.valid_upper
