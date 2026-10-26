"""C1/C2/D3 aggregate method core on authorized, already frozen local inputs.

No data loading, model fitting, threshold selection, clinical TEST execution,
audit export or automatic analysis is performed by this module.
"""
import itertools
import math
from fractions import Fraction as F
import numpy as np
from onsetbounds.bounds import single_policy_bounds, valid_at
from onsetbounds.paired import paired_valid_bounds, marginal_subtraction
from onsetbounds.utility import utility_at

MODELS = ('ELASTIC_NET', 'XGBOOST', 'TREND')
WINDOWS = ((0, 6), (6, 12), (6, 24), (6, 48))
GRID = ('0', '.05', '.10', '.20', '.30', '.40', '.50', '.60', '.75', '.90', '1.00')


def frac(x):
    return F(str(float(x)))


def rank(lo, hi):
    return 'first_superior' if lo > 0 else 'second_superior' if hi < 0 else 'unresolved'


def restrict(L, U, d, direction):
    if direction not in ('early', 'late') or not 0 <= d <= 1:
        raise ValueError('invalid onset restriction')
    return (L + d * (U - L), U) if direction == 'late' else (L, U - d * (U - L))


def calc(rows, d, direction, paired=True):
    lo = hi = 0
    for L, U, a, b in rows:
        l, u = restrict(L, U, d, direction)
        if paired:
            v = paired_valid_bounds(l, u, a, b, F(6), F(24))
            z = (v.lower, v.upper)
        else:
            v = single_policy_bounds(l, u, a, F(6), F(24))
            z = (v.valid_lower, v.valid_upper)
        lo += z[0]
        hi += z[1]
    return (lo, hi)


def changes(rows, direction):
    cuts = {F(0), F(1)}
    for L, U, a, b in rows:
        if U == L:
            continue
        for alert in [a, b]:
            if not math.isfinite(alert):
                continue
            for q in [alert + 6, alert + 24]:
                d = (q - L) / (U - L) if direction == 'late' else (U - q) / (U - L)
                if 0 <= d <= 1:
                    cuts.add(d)
    return sorted(cuts)


def breakdown(rows, direction):
    if not rows:
        raise ValueError('nonempty fixed event population required')
    cuts = changes(rows, direction)

    def at(j):
        return cuts[j // 2] if j % 2 == 0 else (cuts[j // 2] + cuts[j // 2 + 1]) / 2
    cache = {}

    def state(j):
        if j not in cache:
            cache[j] = calc(rows, at(j), direction)
        return cache[j]
    n = 2 * len(cuts) - 1
    if rank(*state(n - 1)) == 'unresolved':
        return (dict(delta_star=None, delta_star_lower_bound=None, delta_star_upper_bound=None, threshold_attained=False, ranking='unresolved', at_threshold_lower=None, at_threshold_upper=None, first_resolved_lower=None, first_resolved_upper=None, endpoint_lower=state(n - 1)[0] / len(rows), endpoint_upper=state(n - 1)[1] / len(rows), first_prespecified_grid_delta=None), None)
    left, right = (0, n - 1)
    while left < right:
        mid = (left + right) // 2
        if rank(*state(mid)) != 'unresolved':
            right = mid
        else:
            left = mid + 1
    j = left
    t = cuts[j // 2]
    z = calc(rows, t, direction)
    zz = state(j)
    scale = 10 ** 12
    floor = t.numerator * scale // t.denominator
    ceil = -(-t.numerator * scale // t.denominator)
    return (dict(delta_star=float(t), delta_star_lower_bound=floor / scale, delta_star_upper_bound=ceil / scale, threshold_attained=j % 2 == 0, ranking=rank(*zz), at_threshold_lower=z[0] / len(rows), at_threshold_upper=z[1] / len(rows), first_resolved_lower=zz[0] / len(rows), first_resolved_upper=zz[1] / len(rows), endpoint_lower=state(n - 1)[0] / len(rows), endpoint_upper=state(n - 1)[1] / len(rows), first_prespecified_grid_delta=next((float(F(d)) for d in GRID if rank(*calc(rows, F(d), direction)) != 'unresolved'))), t)


def _fixed_events(lower, upper, alerts, names):
    lower, upper = np.asarray(lower), np.asarray(upper)
    alerts = np.asarray(alerts)
    if not len(lower) or upper.shape != lower.shape or alerts.shape != (len(lower), len(names)):
        raise ValueError('nonempty aligned fixed-E population required')
    if (not np.isfinite(lower).all() or not np.isfinite(upper).all() or
            np.any(lower > upper) or np.isnan(alerts).any() or np.isneginf(alerts).any()):
        raise ValueError('finite brackets and finite/+inf alerts required')
    return lower, upper, alerts


def _endpoints(lower, upper, alerts, lead_min=6, lead_max=24):
    singles = []
    for column in alerts.T:
        bounds = [single_policy_bounds(l, u, a, lead_min, lead_max)
                  for l, u, a in zip(lower, upper, column)]
        singles.append(np.array([(b.valid_lower, b.valid_upper) for b in bounds]))
    pairs = {}
    for a, b in itertools.combinations(range(alerts.shape[1]), 2):
        bounds = [paired_valid_bounds(l, u, aa, bb, lead_min, lead_max)
                  for l, u, aa, bb in zip(lower, upper, alerts[:, a], alerts[:, b])]
        pairs[a, b] = np.array([(v.lower, v.upper) for v in bounds])
    return singles, pairs


def _ratio(value, denominator):
    return float(value / denominator) if denominator else None


def lead_windows(lower, upper, alerts, names=MODELS, windows=WINDOWS):
    """D3: same original E, brackets and first alerts in every closed window.

    Supply the original primary-evaluable event subset, not a newly selected
    population for each window. Results describe the fixed finite cohort.
    """
    lower, upper, alerts = _fixed_events(lower, upper, alerts, names)
    n = len(lower)
    primary_single, primary_pair = _endpoints(lower, upper, alerts)
    singles, pairs = [], []
    for lo, hi in windows:
        ss, pp = _endpoints(lower, upper, alerts, lo, hi)
        for j, values in enumerate(ss):
            low, high = values.sum(axis=0) / n
            base = primary_single[j].sum(axis=0) / n
            singles.append(dict(model=names[j], lead_min_h=lo, lead_max_h=hi, fixed_E=n,
                lower=low, upper=high, width=high-low, ambiguous_event_fraction=high-low,
                delta_lower_vs_primary=low-base[0], delta_upper_vs_primary=high-base[1],
                delta_width_vs_primary=(high-low)-(base[1]-base[0])))
        for (a, b), values in pp.items():
            low, high = values.sum(axis=0) / n
            ml, mu = marginal_subtraction(*ss[a].mean(axis=0), *ss[b].mean(axis=0))
            base = primary_pair[a, b].mean(axis=0)
            pairs.append(dict(contrast=names[a]+'-'+names[b], lead_min_h=lo, lead_max_h=hi,
                fixed_E=n, direct_lower=low, direct_upper=high, direct_width=high-low,
                marginal_lower=ml, marginal_upper=mu, marginal_width=mu-ml,
                relative_width_reduction=_ratio((mu-ml)-(high-low), mu-ml),
                contains_zero=bool(low <= 0 <= high), ranking=rank(low, high),
                delta_lower_vs_primary=low-base[0], delta_upper_vs_primary=high-base[1],
                delta_width_vs_primary=(high-low)-(base[1]-base[0])))
    return singles, pairs


def alert_concordance(first, second):
    """C1 aggregate alert concordance for one caller-selected frozen population.

    Positive cells <=10 and complementary disclosures use the accepted C1
    withholding rules; patient alert values are never returned.
    """
    first, second = np.asarray(first), np.asarray(second)
    if first.ndim != 1 or first.shape != second.shape or any(
            np.isnan(x).any() or np.isneginf(x).any() for x in (first, second)):
        raise ValueError('aligned finite/+inf alert arrays required')
    fa, fb = np.isfinite(first), np.isfinite(second)
    both = fa & fb
    counts = [int(both.sum()), int((fa & ~fb).sum()),
              int((~fa & fb).sum()), int((~fa & ~fb).sum())]
    delta = first[both] - second[both]
    absolute = np.abs(delta)
    protected = 0 < counts[0] <= 10
    pooled = any(0 < z <= 10 for z in [counts[3], counts[1]+counts[2]])
    directional = any(0 < z <= 10 for z in counts[1:3])
    values = dict(zip(('both_alert', 'only_first', 'only_second', 'neither_alert'), counts))
    values.update(both_alert_jaccard=_ratio(counts[0], sum(counts[:3])),
        jaccard_denominator_either_alerted=sum(counts[:3]), only_one_alert_pooled=counts[1]+counts[2],
        signed_median_h=float(np.median(delta)) if len(delta) else None,
        signed_Q1_h=float(np.quantile(delta, .25)) if len(delta) else None,
        signed_Q3_h=float(np.quantile(delta, .75)) if len(delta) else None,
        median_absolute_h=float(np.median(absolute)) if len(delta) else None)
    for h, key in [(0, 'same_hour'), (1, 'within_1h'), (3, 'within_3h'), (6, 'within_6h')]:
        values[key+'_proportion'] = _ratio(int((absolute <= h).sum()), len(delta))
    complementary = {'only_first', 'only_second', 'neither_alert', 'only_one_alert_pooled',
                     'both_alert_jaccard', 'jaccard_denominator_either_alerted'}
    return dict(N=len(first), status=('withheld_partition_small_cell' if protected else
        'complementary_partition_cells_withheld' if pooled else
        'directional_cells_withheld' if directional else 'released'), **{
            key: 'withheld' if protected or (pooled and key in complementary) or
            (directional and key in {'only_first', 'only_second'}) else value
            for key, value in values.items()})


def shared_onset(lower, upper, alerts, names=MODELS):
    """C1 paired width decomposition and coarse/fine gain strata on fixed E."""
    lower, upper, alerts = _fixed_events(lower, upper, alerts, names)
    ss, pp = _endpoints(lower, upper, alerts)
    results, strata = [], []
    for (a, b), paired in pp.items():
        contrast = names[a]+'-'+names[b]
        wm = ss[a][:, 1]-ss[a][:, 0]+ss[b][:, 1]-ss[b][:, 0]
        wd = paired[:, 1]-paired[:, 0]
        gain = wm-wd
        ml, mu = marginal_subtraction(*ss[a].mean(axis=0), *ss[b].mean(axis=0))
        ua = np.array([valid_at(u, alert, 6, 24) for u, alert in zip(upper, alerts[:, a])])
        ub = np.array([valid_at(u, alert, 6, 24) for u, alert in zip(upper, alerts[:, b])])
        results.append(dict(contrast=contrast, fixed_E=len(lower),
            direct_lower=paired[:, 0].mean(), direct_upper=paired[:, 1].mean(),
            marginal_lower=ml, marginal_upper=mu, mean_marginal_width=wm.mean(),
            mean_direct_width=wd.mean(), mean_absolute_gain=gain.mean(), sum_gain=int(gain.sum()),
            relative_reduction=_ratio(gain.mean(), wm.mean()), u_paired_point=float((ua-ub).mean())))
        fa, fb = np.isfinite(alerts[:, a]), np.isfinite(alerts[:, b])
        joint = fa & fb
        delta = np.full(len(lower), np.inf)
        delta[joint] = np.abs(alerts[joint, a]-alerts[joint, b])
        bracket = upper-lower
        ambiguous = (ss[a][:, 1]-ss[a][:, 0]).astype(bool) & (ss[b][:, 1]-ss[b][:, 0]).astype(bool)
        families = {
            'alert_timing': np.where(joint, np.where(delta == 0, 'both_same_hour',
                np.where(delta <= 3, 'both_difference_(0,3]h', 'both_difference_>3h')),
                np.where(fa | fb, 'only_one_alert', 'neither_alert')),
            'onset_bracket_width': np.where(bracket <= 6, '<=6h', np.where(bracket <= 12,
                '(6,12]h', np.where(bracket <= 24, '(12,24]h', '>24h'))),
            'alert_timing_coarse': np.where(joint & (delta <= 3), 'both_difference_<=3h', 'all_other_alert_patterns'),
            'onset_bracket_width_coarse': np.where(bracket <= 12, '<=12h', '>12h'),
            'joint_valid_ambiguity': np.where(ambiguous, 'both_valid_classifications_ambiguous', 'at_least_one_classification_fixed'),
            'joint_ambiguity_geometry': np.where(ambiguous, np.where(wd == 0,
                'both_ambiguous_identical_valid_sets', np.where(wd == 1,
                'both_ambiguous_nested_valid_sets', 'both_ambiguous_both_exclusive_regions')),
                'at_least_one_classification_fixed')}
        for family, labels in families.items():
            categories = sorted(set(labels))
            protected = any(0 < int((labels == c).sum()) <= 10 for c in categories)
            for category in categories:
                mask = labels == category
                values = dict(n=int(mask.sum()), mean_marginal_width=wm[mask].mean(),
                    mean_direct_width=wd[mask].mean(), mean_absolute_gain=gain[mask].mean(),
                    sum_gain=int(gain[mask].sum()), share_total_gain=_ratio(gain[mask].sum(), gain.sum()))
                strata.append(dict(contrast=contrast, family=family, stratum=category,
                    family_denominator=len(lower), status='withheld_family_small_cell' if protected else 'released',
                    **{key: 'withheld' if protected else value for key, value in values.items()}))
    return results, strata


def right_endpoint_reference(D, E, upper, alerts):
    """C1 U reference for a single policy, retaining all three denominators.

    Pass the complete frozen cohort; event bracket upper endpoints and alert
    times are hours since each patient's fixed monitoring origin. +inf is no
    alert. E is the original primary-evaluable mask, a subset of D.
    """
    D, E = np.asarray(D, bool), np.asarray(E, bool)
    upper, alerts = np.asarray(upper), np.asarray(alerts)
    if (D.ndim != 1 or not len(D) or any(x.shape != D.shape for x in (E, upper, alerts))
            or np.any(E & ~D) or not np.isfinite(upper[D]).all()
            or np.isnan(alerts).any() or np.isneginf(alerts).any()):
        raise ValueError('aligned frozen event/evaluability and alert arrays required')
    finite = np.isfinite(alerts)
    valid = np.zeros(len(D), dtype=int)
    premature = np.zeros(len(D), dtype=int)
    missed = np.zeros(len(D), dtype=int)
    utility = -finite.astype(float)
    for i in np.flatnonzero(D):
        b = single_policy_bounds(upper[i], upper[i], alerts[i])
        valid[i], premature[i], missed[i] = b.valid_lower, b.premature_lower, b.late_missed_lower
        utility[i] = utility_at(upper[i], alerts[i])
    return dict(valid=_ratio(valid[E].sum(), E.sum()), premature=_ratio(premature[E].sum(), E.sum()),
                late_missed=_ratio(missed[E].sum(), E.sum()), ppv=_ratio(valid.sum(), finite.sum()),
                utility=float(utility.mean()), E=int(E.sum()), alerted=int(finite.sum()), N=len(D))


def onset_location(lower, upper, alerts, names=MODELS):
    """C2 single/paired restricted-support curves and strict-ranking thresholds.

    Original E and alerts stay fixed. Decimal hour values are converted to exact
    rationals for the closed-point/open-cell change-point search.
    """
    lower, upper, alerts = _fixed_events(lower, upper, alerts, names)
    columns = [[(frac(l), frac(u), frac(a) if math.isfinite(a) else math.inf, math.inf)
                for l, u, a in zip(lower, upper, column)] for column in alerts.T]
    singles, pairs, thresholds = [], [], []
    n = len(lower)
    for direction in ('early', 'late'):
        for j, rows in enumerate(columns):
            base = calc(rows, F(0), direction, paired=False)
            basewidth = (base[1]-base[0])/n
            for ds in GRID:
                low, high = np.array(calc(rows, F(ds), direction, paired=False))/n
                singles.append(dict(model=names[j], direction=direction, delta=float(F(ds)),
                    fixed_E=n, lower=low, upper=high, width=high-low,
                    relative_width_reduction=_ratio(basewidth-(high-low), basewidth)))
        for a, b in itertools.combinations(range(len(names)), 2):
            rows = [(l, u, aa, brow[2]) for (l, u, aa, _), brow in zip(columns[a], columns[b])]
            base = calc(rows, F(0), direction)
            basewidth = (base[1]-base[0])/n
            for ds in GRID:
                low, high = np.array(calc(rows, F(ds), direction))/n
                sa = next(s for s in singles if s['model'] == names[a] and s['direction'] == direction and s['delta'] == float(F(ds)))
                sb = next(s for s in singles if s['model'] == names[b] and s['direction'] == direction and s['delta'] == float(F(ds)))
                ml, mu = marginal_subtraction(sa['lower'], sa['upper'], sb['lower'], sb['upper'])
                pairs.append(dict(contrast=names[a]+'-'+names[b], direction=direction, delta=float(F(ds)),
                    fixed_E=n, direct_lower=low, direct_upper=high, direct_width=high-low,
                    marginal_lower=ml, marginal_upper=mu, ranking=rank(low, high),
                    relative_width_reduction=_ratio(basewidth-(high-low), basewidth)))
            summary, _ = breakdown(rows, direction)
            thresholds.append(dict(contrast=names[a]+'-'+names[b], direction=direction, fixed_E=n, **summary))
    return singles, pairs, thresholds
