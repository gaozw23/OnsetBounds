import itertools

from onsetbounds.validation import analytic_matches_enumeration


def test_analytic_vs_cell_enumeration():
    for lower, width, alert in itertools.product([0, 5, 10], [0, 3, 12], [-4, 0, 8, float("inf")]):
        assert analytic_matches_enumeration(lower, lower + width, alert)
