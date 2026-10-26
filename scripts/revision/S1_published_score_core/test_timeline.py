"""Entirely synthetic temporal tests for the locked retrospective policy."""
import math,unittest
import numpy as np
from timeline import (Lab,hr,HOUR,FORMULA,run_policy,baseline,outcome,sensitivity,
                      active_location,active_service,detect_any_aki,latest_inputs)

STATIC=dict(age=60.,male=1.,black=0.)
BASE=dict(sodium=138.,potassium=4.,chloride=103.,bicarbonate=24.,bun=20.,creatinine=1.,hemoglobin=12.)
def labs(at=0,**change):return [Lab(v,hr(at),hr(at),x) for v,x in dict(BASE,**change).items()]
def replay(rs,services=None,intervals=None,**kwargs):
    return run_policy(rs,intervals or [(hr(0),hr(80),1.)],services or [(hr(0),0.)],STATIC,hr(0),hr(24),hr(72),**kwargs)

class TemporalTests(unittest.TestCase):
    def test_missing_means_no_score_not_low_risk(self):
        rs=[r for r in labs() if r.variable!='hemoglobin']+[Lab('hemoglobin',hr(26),hr(27),12.)]
        p=replay(rs)
        self.assertFalse(p['landmark_complete']);self.assertEqual(p['times'][0],hr(27))
        self.assertGreater(p['missing_pre'],0);self.assertIsNone(p['alert'])
    def test_atomic_location_service_updates(self):
        rs=labs();p=replay(rs,services=[(0,0.),(hr(10),1.)],intervals=[(0,hr(10),0.),(hr(10),hr(80),1.)])
        indices={t:i for i,t in enumerate(p['times'])}
        a=p['risks'][indices[0]];b=p['risks'][indices[hr(10)]]
        self.assertAlmostEqual(math.log(b/(1-b))-math.log(a/(1-a)),.09760579+.2331921,places=12)
        self.assertTrue(np.isnan(active_location([(0,hr(1),1.),(hr(2),hr(3),1.)],hr(1.5))))
    def test_positive_batch_stops_before_score(self):
        rs=labs()+[Lab('creatinine',hr(30),hr(30),1.4),Lab('bun',hr(30),hr(30),1000.)]
        p=replay(rs)
        self.assertEqual(p['stop'],hr(30));self.assertNotIn(hr(30),p['times']);self.assertIsNone(p['alert'])
    def test_preorigin_first_alert_never_resets(self):
        rs=labs()+[Lab('bun',hr(5),hr(5),1000.),Lab('bun',hr(25),hr(25),20.)]
        p=replay(rs)
        self.assertEqual(p['alert'],hr(5));self.assertTrue(p['pre_scorable'])
    def test_stage1_stop_does_not_censor_stage2_outcome(self):
        rs=labs()+[Lab('creatinine',hr(30),hr(30),1.4),Lab('creatinine',hr(40),hr(42),2.1)]
        p=replay(rs);c0=baseline(rs,0,hr(24));e=outcome(rs,c0,hr(24),hr(72))
        self.assertEqual(p['stop'],hr(30));self.assertEqual(e['U'],hr(40))
        self.assertEqual(e['L'],hr(30));self.assertTrue(e['E'])
    def test_measurement_and_availability_clocks_alert_after_U(self):
        rs=labs()+[Lab('creatinine',hr(30),hr(30),1.1),Lab('creatinine',hr(40),hr(44),2.1),Lab('bun',hr(42),hr(42),1000.)]
        p=replay(rs);e=outcome(rs,1.,hr(24),hr(72))
        self.assertEqual(p['alert'],hr(42));self.assertEqual(p['stop'],hr(44))
        self.assertEqual(e['U'],hr(40));self.assertEqual(sensitivity(e,p['alert']),(0,0,0,0))
    def test_late_reference_stops_now_without_rewriting_alert(self):
        rs=labs(bun=1000.)+[Lab('creatinine',hr(-1),hr(5),.6,False)]
        p=replay(rs)
        self.assertEqual(p['alert'],0);self.assertEqual(p['stop'],hr(5))
        self.assertIn(0,p['times']);self.assertNotIn(hr(5),p['times'])
    def test_baseline_uses_only_available_first24h(self):
        rs=[Lab('creatinine',0,0,1.1),Lab('creatinine',hr(10),hr(26),.8)]
        self.assertEqual(baseline(rs,0,hr(24)),1.1)
        x,provenance=latest_inputs(rs,hr(24));self.assertEqual(x['creatinine'],1.1)
        self.assertTrue(all(o<=hr(24) and a<=hr(24) for o,a in provenance))
    def test_duplicates_only_currently_known_median(self):
        rs=[Lab('creatinine',0,0,1.),Lab('creatinine',0,hr(28),1.2)]
        self.assertEqual(baseline(rs,0,hr(24)),1.)
        self.assertEqual(latest_inputs(rs,hr(28))[0]['creatinine'],1.1)
    def test_fixed_L_U_E_and_unknown_outcome(self):
        rs=[Lab('creatinine',hr(31),hr(32),1.),Lab('creatinine',hr(40),hr(44),2.)]
        e=outcome(rs,1.,hr(24),hr(72));self.assertTrue(e['E'])
        self.assertEqual((e['L'],e['U']),(hr(31),hr(40)))
        self.assertEqual(sensitivity(e,hr(20)),(1,1,1,1))
        self.assertEqual(outcome(rs[1:],1.,hr(24),hr(72))['finite'],False)
        self.assertIsNone(outcome([],1.,hr(24),hr(72))['D'])
        self.assertEqual(outcome(rs[:1],1.,hr(24),hr(72))['D'],0)
    def test_unknown_history_window_and_exact_boundary(self):
        positive,unknown,n=detect_any_aki([Lab('creatinine',0,0,1.)],0,0)
        self.assertFalse(positive);self.assertEqual(unknown,1)
        rs=[Lab('creatinine',0,0,1.),Lab('creatinine',hr(48),hr(48),1.4)]
        self.assertTrue(detect_any_aki(rs,hr(48),0)[0])
        rs=[rs[0],Lab('creatinine',hr(48)+1,hr(48)+1,1.4)]
        self.assertFalse(detect_any_aki(rs,hr(48)+1,0)[0])
    def test_initial_Cr_strict4_clock_and_no_lookback_replacement(self):
        p=replay(labs(creatinine=4.));self.assertFalse(p['initial_gt4'])
        p=replay(labs(creatinine=4.1));self.assertTrue(p['initial_gt4']);self.assertFalse(p['ever_scorable'])
        rs=labs()+[Lab('creatinine',hr(.1),hr(10),5.)]
        p=replay(rs);self.assertEqual(p['initial_cr'],1.)
    def test_latest_service_and_source_missing_state(self):
        self.assertTrue(np.isnan(active_service([(hr(2),1.)],hr(1))))
        self.assertEqual(active_service([(0,0.),(hr(2),1.)],hr(1)),0)
        self.assertTrue(np.isnan(active_service([(0,0.),(0,1.)],0)))
    def test_outcome_not_a_policy_argument(self):
        import inspect
        params=inspect.signature(run_policy).parameters
        for k in ['U','D','C0','outcome']:self.assertNotIn(k,params)
        p=replay(labs());self.assertEqual(p['timeline_violations'],0)
        self.assertTrue(all(t<hr(72) for t in p['times']))

if __name__=='__main__':unittest.main(verbosity=2)
