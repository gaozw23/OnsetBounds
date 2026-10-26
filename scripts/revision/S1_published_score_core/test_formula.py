"""Fictitious inputs only; independent transcription, not author predictions."""
from decimal import Decimal, localcontext
import math, unittest
import numpy as np
from ugwuowo_2020_model import COEFFICIENTS, UNITS, score, score_scalar, expit, above_published_threshold, first_alert
import pandas as pd
FICTIONAL = {'age': 60.0, 'bun': 20.0, 'creatinine': 1.0, 'potassium': 4.0, 'bicarbonate': 24.0, 'chloride': 103.0, 'hemoglobin': 12.0, 'sodium': 138.0, 'male': 1, 'black': 0, 'icu': 0, 'surgical': 0}

def independent(x):
    with localcontext() as ctx:
        ctx.prec = 70
        d = lambda k: Decimal(str(x[k]))
        y = Decimal('.0090205') * d('age') + Decimal('.0028409') * d('bun') + Decimal('.3149055') * d('creatinine') + Decimal('.2605861') * d('potassium') - Decimal('.0768276') * d('bicarbonate') - Decimal('.0606976') * d('chloride') - Decimal('.0115591') * d('hemoglobin') + Decimal('.0442439') * d('sodium') + Decimal('.0274777') * d('male') + Decimal('.155564') * d('black') + Decimal('.09760579') * d('icu') + Decimal('.2331921') * d('surgical') - Decimal('3.333648')
        return (float(y), float(y.exp() / (1 + y.exp())))

class ModelTests(unittest.TestCase):

    def test_independent_decimal(self):
        for f in [FICTIONAL, dict(FICTIONAL, age=75, bun=40, creatinine=2.0, black=1, icu=1, surgical=1), {k: 0.0 for k in FICTIONAL}]:
            y, p = independent(f)
            r = score_scalar(f, UNITS)
            self.assertAlmostEqual(r['linear_predictor'], y, places=13)
            self.assertAlmostEqual(r['risk'], p, places=14)

    def test_all_coefficients_intercept_signs_no_scaling(self):
        independent_coeffs = [0.0090205, 0.0028409, 0.3149055, 0.2605861, -0.0768276, -0.0606976, -0.0115591, 0.0442439, 0.0274777, 0.155564, 0.09760579, 0.2331921]
        base = {k: 0.0 for k in FICTIONAL}
        r0 = score_scalar(base, UNITS)
        self.assertEqual(r0['linear_predictor'], -3.333648)
        for v, b in zip(FICTIONAL, independent_coeffs):
            with self.subTest(variable=v):
                q = score_scalar(dict(base, **{v: 1.0}), UNITS)
                self.assertAlmostEqual(q['linear_predictor'] - r0['linear_predictor'], b, places=13)
                self.assertEqual(q['risk'] > r0['risk'], b > 0)

    def test_missing_units_nonfinite(self):
        for v in FICTIONAL:
            f = dict(FICTIONAL)
            del f[v]
            self.assertEqual(score_scalar(f, UNITS)['status'], 'UNSCORABLE')
            for value in [None, float('nan'), float('inf'), 'unknown']:
                r = score_scalar(dict(FICTIONAL, **{v: value}), UNITS)
                self.assertEqual(r['status'], 'UNSCORABLE')
                self.assertTrue(math.isnan(r['risk']))
            u = dict(UNITS)
            del u[v]
            self.assertEqual(score_scalar(FICTIONAL, u)['status'], 'UNSCORABLE')
            self.assertEqual(score_scalar(FICTIONAL, dict(UNITS, **{v: 'wrong'}))['status'], 'UNSCORABLE')

    def test_binary_domain(self):
        for v in ['male', 'black', 'icu', 'surgical']:
            for value in [-1, 0.5, 2]:
                self.assertEqual(score_scalar(dict(FICTIONAL, **{v: value}), UNITS)['status'], 'UNSCORABLE')

    def test_does_not_accept_2019_extra_inputs(self):
        for v in ['wbc', 'platelets', 'creatinine_change_48h', 'creatinine_baseline']:
            self.assertEqual(score_scalar(dict(FICTIONAL, **{v: 1}), UNITS)['status'], 'UNSCORABLE')

    def test_extreme_probability_range(self):
        self.assertTrue(np.array_equal(expit([-1e+300, -1000, 0, 1000, 1e+300]), [0, 0, 0.5, 1, 1]))
        for value in [100000000.0, 1e+300]:
            r = score_scalar(dict(FICTIONAL, bun=value), UNITS)
            self.assertTrue(0 <= r['risk'] <= 1)

    def test_batch_scalar(self):
        fs = [dict(FICTIONAL, bun=v) for v in [0, 10, 20, 100]]
        b = score(fs, UNITS)
        for i, f in enumerate(fs):
            r = score_scalar(f, UNITS)
            self.assertEqual(b.iloc[i].status, r['status'])
            self.assertEqual(b.iloc[i].reason, r['reason'])
            self.assertAlmostEqual(b.iloc[i].linear_predictor, r['linear_predictor'], places=13)
            self.assertAlmostEqual(b.iloc[i].risk, r['risk'], places=14)

    def test_threshold_strict_boundary(self):
        self.assertFalse(above_published_threshold(0.15))
        self.assertFalse(above_published_threshold(np.nextafter(0.15, 0)))
        self.assertTrue(above_published_threshold(np.nextafter(0.15, 1)))
        self.assertIsNone(above_published_threshold(None))
        self.assertIsNone(above_published_threshold(float('nan')))

    def test_first_alert_no_reset_and_stop(self):
        x = [(0, None), (1, 0.15), (2, 0.151), (3, 0.9)]
        self.assertEqual(first_alert(x, 0, 4), {'status': 'ALERTED', 'alert_time': 2})
        self.assertEqual(first_alert(x, 0, 4, aki_stop=2)['status'], 'POLICY_STOPPED_NO_ALERT')
        self.assertEqual(first_alert([(0, None)], 0, 4)['status'], 'NEVER_SCORABLE')
        self.assertEqual(first_alert([(0, 0.1)], 0, 4)['status'], 'NO_THRESHOLD_CROSSING')
        self.assertEqual(first_alert([(0, None)], 0, 4, interrupted=True)['status'], 'MONITORING_INTERRUPTED_NO_ALERT')
        self.assertEqual(first_alert(x, 0, 4, aki_stop=-1)['status'], 'AKI_BEFORE_MONITORING')

    def test_threshold_logit_arithmetic(self):
        y = math.log(0.15 / 0.85)
        self.assertAlmostEqual(float(expit(y)), 0.15, places=15)
if __name__ == '__main__':
    unittest.main(verbosity=2)
