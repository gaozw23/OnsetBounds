"""Synthetic tests; never access clinical records."""
import inspect,unittest
import numpy as np
import pandas as pd
import timeline as T
from policy_math import (choose_threshold,first_landmark_alert,event_metrics,
    scoring_opportunities,IndependentFormulaAudit,threshold_independent_member)

class PolicyTests(unittest.TestCase):
    def test_calibration_patient_level_strict_and_rule(self):
        a=[np.array([.001,(i+1)/100]) for i in range(10)]
        r=choose_threshold(a)
        self.assertEqual(r['threshold'],.09);self.assertEqual(r['alerts'],1)
        self.assertEqual(r['achieved_rate'],.1)
        self.assertEqual(list(inspect.signature(choose_threshold).parameters),['non_event_scores'])
    def test_independent_brute_force(self):
        arrays=[np.array(x) for x in [[.1,.2],[.2],[.15,.3],[],[.05],[.07],[.08],[.09],[.1],[.11],[.12]]]
        candidates=sorted({-np.inf,1.}|{float(x) for a in arrays for x in a})
        trial=[]
        for c in candidates:
            k=sum(any(x>c for x in a) for a in arrays)
            trial.append((abs(10*k-len(arrays)),-c,c,k))
        _,_,c,k=min(trial)
        r=choose_threshold(arrays)
        self.assertEqual((r['threshold'],r['alerts']),(c,k))
    def test_ties_choose_higher_threshold(self):
        # N5: rates0 and.2 equally close to.1; higher no-alert extreme wins.
        r=choose_threshold([np.array([x]) for x in [.1,.2,.3,.4,.5]])
        self.assertEqual(r['threshold'],1.);self.assertEqual(r['alerts'],0)
    def test_missing_scores_retained_in_denominator(self):
        r=choose_threshold([np.array([]) for _ in range(8)]+[np.array([.1]),np.array([.2])])
        self.assertEqual(r['non_events'],10);self.assertEqual(r['alerts'],1)
        self.assertEqual(r['threshold'],.1)
    def test_episode_origin_strict_boundary(self):
        self.assertEqual(first_landmark_alert([10,20,30],[.15,.15001,.9],.15),20)
        self.assertIsNone(first_landmark_alert([10],[.15],.15))
    def test_opportunity_not_first_alert_validity(self):
        e=dict(D=1,finite=True,E=True,L=T.hr(36),U=T.hr(52))
        times=[T.hr(0),T.hr(20),T.hr(30)]
        self.assertTrue(scoring_opportunities(e,times)['possible'])
        a=first_landmark_alert(times,[.2,.1,.3],.15)
        self.assertEqual(a,0);self.assertEqual(event_metrics(e,a)['valid'],(0,0,0,0))
    def test_points_and_partitions(self):
        e=dict(D=1,finite=True,E=True,L=T.hr(36),U=T.hr(52))
        m=event_metrics(e,T.hr(20))
        self.assertEqual(m['valid'],(0,1,0,1))
        for x in m.values():self.assertTrue(x[0]<=x[2]<=x[1] and x[0]<=x[3]<=x[1])
        self.assertEqual(event_metrics(e,None)['late_or_missed'],(1,1,1,1))
    def test_absolute_timestamp_offset_preserves_opportunity(self):
        offset=6_000_000_000_000_000_000
        a=dict(D=1,finite=True,E=True,L=T.hr(36),U=T.hr(52))
        b=dict(a,L=a['L']+offset,U=a['U']+offset)
        self.assertEqual(scoring_opportunities(a,[T.hr(20),T.hr(30)]),
                         scoring_opportunities(b,[offset+T.hr(20),offset+T.hr(30)]))
    def test_new_cohort_independent_of_preorigin_alert(self):
        priority=['age','known_aki','first_alert_before_t0']
        self.assertTrue(threshold_independent_member(dict(age=False,known_aki=False,first_alert_before_t0=True),priority))
        self.assertFalse(threshold_independent_member(dict(age=True,known_aki=False,first_alert_before_t0=False),priority))
    def test_real_formula_auditor_on_fictitious_rows(self):
        row=dict(age=60.,bun=20.,creatinine=1.,potassium=4.,bicarbonate=24.,chloride=103.,
                 hemoglobin=12.,sodium=138.,male=1.,black=0.,icu=1.,surgical=0.)
        a=IndependentFormulaAudit();r=a(pd.DataFrame([row]),T.FORMULA.UNITS)
        self.assertEqual(a.rows,1);self.assertLess(a.max_p_error,1e-14)
        self.assertEqual(r.status.iloc[0],'SCORED_FIXED_FORMULA')

if __name__=='__main__':unittest.main(verbosity=2)
