from onsetbounds.bounds import single_policy_bounds
from onsetbounds.paired import marginal_subtraction, paired_valid_bounds


def test_identical_policy_direct_zero():
    b = paired_valid_bounds(0, 10, 0, 0, 0, 5)
    assert (b.lower, b.upper) == (0, 0)


def test_direct_full_minus_one_to_one():
    b = paired_valid_bounds(0, 4, 0, 2, 0, 2)
    assert (b.lower, b.upper) == (-1, 1)


def test_certain_plus_and_minus_one():
    plus = paired_valid_bounds(10, 12, 4, float("inf"))
    minus = paired_valid_bounds(10, 12, float("inf"), 4)
    assert (plus.lower, plus.upper) == (1, 1)
    assert (minus.lower, minus.upper) == (-1, -1)


def test_direct_inside_marginal_subtraction():
    lower, upper, a, b = 0, 15, 0, 5
    direct = paired_valid_bounds(lower, upper, a, b, 2, 8)
    ba = single_policy_bounds(lower, upper, a, 2, 8)
    bb = single_policy_bounds(lower, upper, b, 2, 8)
    marginal = marginal_subtraction(ba.valid_lower, ba.valid_upper, bb.valid_lower, bb.valid_upper)
    assert marginal[0] <= direct.lower <= direct.upper <= marginal[1]
