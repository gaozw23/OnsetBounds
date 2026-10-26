"""Prespecified D0-only burden calibration and frozen-alert evaluation."""
import math
import numpy as np
import timeline as T
HOUR=T.HOUR

def threshold_independent_member(flags,priorities):
    assert priorities[-1]=='first_alert_before_t0'
    return all(flags[k] is False for k in priorities[:-1])

def post_scores(p,t0,H):
    mask=np.array([t0<=t<H for t in p['times']],dtype=bool)
    return np.asarray(p['times'],dtype=np.int64)[mask],p['risks'][mask]

def first_landmark_alert(times,risks,c):
    hits=np.flatnonzero(np.asarray(risks)>c)
    return int(times[int(hits[0])]) if len(hits) else None

def choose_threshold(non_event_scores):
    """Only D0 patient score arrays enter. No event/onset/time/TEST inputs.

    Enumerate distinct D0 risk values plus extreme policies. Missing-score
    patients stay in denominator with max=-inf. Exact integer rate objective.
    """
    if not non_event_scores:raise ValueError('NO_OBSERVED_VALIDATION_NON_EVENTS')
    n=len(non_event_scores)
    arrays=[np.asarray(a,dtype=float) for a in non_event_scores]
    if any(not np.isfinite(a).all() or not ((a>=0)&(a<=1)).all() for a in arrays):
        raise ValueError('NONFINITE_OR_INVALID_CALIBRATION_RISK')
    pooled=np.concatenate(arrays) if any(len(a) for a in arrays) else np.array([])
    if not len(pooled):raise ValueError('NO_SCORABLE_VALIDATION_NON_EVENTS')
    distinct=np.unique(pooled)
    candidates=np.unique(np.concatenate(([-np.inf],distinct,[1.])))
    maxima=np.sort([float(np.max(a)) if len(a) else -np.inf for a in arrays])
    counts=n-np.searchsorted(maxima,candidates,side='right')
    errors=np.abs(10*counts-n)
    best=int(np.flatnonzero(errors==errors.min())[-1])
    c=float(candidates[best]);alerts=int(counts[best])
    return dict(threshold=c,alerts=alerts,non_events=n,achieved_rate=alerts/n,
      absolute_rate_error=abs(alerts/n-.1),objective_integer=int(errors[best]),
      distinct_risk_values=len(distinct),candidate_thresholds=len(candidates),
      minimizing_candidates=int((errors==errors.min()).sum()))

def event_metrics(event,alert):
    """Sharp intervals and U/midpoints, closed6–24h, same E."""
    if not event['E']:return None
    v=T.sensitivity(event,alert)
    if alert is None:p=(0,0,0,0);m=(1,1,1,1)
    else:
        l,u=event['L'],event['U'];a=alert
        p=(int(l-a>24*HOUR),int(u-a>24*HOUR),int(u-a>24*HOUR),int(l+u-2*a>48*HOUR))
        m=(int(u-a<6*HOUR),int(l-a<6*HOUR),int(u-a<6*HOUR),int(l+u-2*a<12*HOUR))
    for x in [v,p,m]:assert x[0]<=x[2]<=x[1] and x[0]<=x[3]<=x[1]
    assert v[2]+p[2]+m[2]==1 and v[3]+p[3]+m[3]==1
    return {'valid':v,'premature':p,'late_or_missed':m}

def scoring_opportunities(event,times):
    """Actual legal score times, not first-alert validity; no threshold fit."""
    if event['D']!=1:return None
    a=np.asarray(times,dtype=np.int64);u=event['U'];l=event['L']
    out=dict(U=bool(((a>=u-24*HOUR)&(a<=u-6*HOUR)).any()),possible=None,guaranteed=None,midpoint=None)
    if not event['finite']:return out
    out['possible']=bool(((a>=l-24*HOUR)&(a<=u-6*HOUR)).any())
    out['guaranteed']=bool(((a>=u-24*HOUR)&(a<=l-6*HOUR)).any())
    # Use relative coordinates; doubling absolute deidentified timestamps can
    # overflow int64 even though all clinical time differences are small.
    delta=a-l;width=u-l
    out['midpoint']=bool(((2*delta>=width-48*HOUR)&(2*delta<=width-12*HOUR)).any())
    return out

class IndependentFormulaAudit:
    """Separately transcribed fsum/libm audit on every actual complete row.

    The wrapper never changes coefficients/values/returned predictions.
    Feature rows stay in process memory and only aggregate errors persist.
    """
    def __init__(self):
        self.original=T.FORMULA.score;self.rows=0;self.max_y_error=0.;self.max_p_error=0.;self.invalid_encoding=0
    def __call__(self,values,units):
        result=self.original(values,units)
        cols=['age','bun','creatinine','potassium','bicarbonate','chloride','hemoglobin','sodium','male','black','icu','surgical']
        coeff=[.0090205,.0028409,.3149055,.2605861,-.0768276,-.0606976,-.0115591,.0442439,.0274777,.155564,.09760579,.2331921]
        ys=[];ps=[]
        for row in values[cols].itertuples(index=False,name=None):
            self.invalid_encoding+=int(not(18<=row[0]<=89 and all(x in (0.,1.) for x in row[8:])))
            y=math.fsum([-3.333648]+[x*b for x,b in zip(row,coeff)])
            p=1/(1+math.exp(-y)) if y>=0 else math.exp(y)/(1+math.exp(y))
            ys.append(y);ps.append(p)
        self.rows+=len(ys)
        if ys:
            self.max_y_error=max(self.max_y_error,float(np.max(np.abs(result.linear_predictor.to_numpy()-ys))))
            self.max_p_error=max(self.max_p_error,float(np.max(np.abs(result.risk.to_numpy()-ps))))
        assert self.max_y_error<5e-13 and self.max_p_error<5e-14 and self.invalid_encoding==0
        return result
    def summary(self):return dict(rows=self.rows,max_y_error=self.max_y_error,max_probability_error=self.max_p_error,invalid_encoding=self.invalid_encoding)
