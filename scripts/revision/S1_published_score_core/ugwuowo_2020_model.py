"""Ugwuowo 2020 Box 1, fixed formula. Never import Simonov coefficients.

Formula scoring is distinct from EHR input and risk-set eligibility gates.
No fitting, calibration, imputation or online scoring service.
"""
import numpy as np
import pandas as pd

INTERCEPT=-3.333648
COEFFICIENTS={
 'age':.0090205,'bun':.0028409,'creatinine':.3149055,
 'potassium':.2605861,'bicarbonate':-.0768276,'chloride':-.0606976,
 'hemoglobin':-.0115591,'sodium':.0442439,'male':.0274777,
 'black':.155564,'icu':.09760579,'surgical':.2331921,
}
UNITS={'age':'years','bun':'mg/dL','creatinine':'mg/dL',
       'potassium':'mEq/L','bicarbonate':'mEq/L','chloride':'mEq/L',
       'hemoglobin':'g/dL','sodium':'mEq/L',
       'male':'binary','black':'binary','icu':'binary','surgical':'binary'}
THRESHOLD=.15

def expit(y):
    y=np.asarray(y,dtype=float);r=np.empty_like(y);pos=y>=0
    r[pos]=1/(1+np.exp(-y[pos]))
    z=np.exp(y[~pos]);r[~pos]=z/(1+z)
    return r

def score(values,units):
    """Return explicit status, fixed linear predictor and probability.

    Inputs are raw Box 1 units, without BUN/10 or platelets/50. All four
    indicator variables must be explicit 0/1, never unknown-to-zero.
    Additional columns are rejected to prevent 2019 predictors being added.
    This validates arithmetic input only, not EHR as-of or cohort semantics.
    """
    f=pd.DataFrame([values]) if isinstance(values,dict) else pd.DataFrame(values).copy()
    reasons=[[] for _ in range(len(f))];matrix=np.full((len(f),len(COEFFICIENTS)),np.nan)
    extra=sorted(set(f.columns)-set(COEFFICIENTS))
    if extra:
        for x in reasons:x.append('unexpected_columns:'+','.join(extra))
    for j,v in enumerate(COEFFICIENTS):
        if v not in f:
            for x in reasons:x.append('missing_column:'+v)
            continue
        if units.get(v)!=UNITS[v]:
            for x in reasons:x.append('unit_mismatch:'+v)
        data=pd.to_numeric(f[v],errors='coerce').to_numpy(dtype=float,na_value=np.nan)
        matrix[:,j]=data
        for i,x in enumerate(data):
            if not np.isfinite(x):reasons[i].append('missing_or_nonfinite:'+v)
            elif v in {'male','black','icu','surgical'} and x not in (0.,1.):
                reasons[i].append('invalid_binary:'+v)
            elif v=='age' and x<0:reasons[i].append('invalid_age')
    valid=np.array([not x for x in reasons],dtype=bool)
    y=np.full(len(f),np.nan);p=np.full(len(f),np.nan)
    with np.errstate(over='ignore',invalid='ignore'):
        y[valid]=INTERCEPT+matrix[valid]@np.array(list(COEFFICIENTS.values()))
    for i in np.flatnonzero(valid&~np.isfinite(y)):
        reasons[i].append('nonfinite_linear_predictor');valid[i]=False
    y[~valid]=np.nan;p[valid]=expit(y[valid])
    return pd.DataFrame(dict(status=np.where(valid,'SCORED_FIXED_FORMULA','UNSCORABLE'),
                             linear_predictor=y,risk=p,reason=[';'.join(x) for x in reasons]),index=f.index)

def score_scalar(values,units):return score(values,units).iloc[0].to_dict()

def above_published_threshold(risk):
    """Strict >15%; None denotes invalid/unscorable, not 'no alert'."""
    if risk is None:return None
    try:p=float(risk)
    except (TypeError,ValueError):return None
    if not np.isfinite(p) or p<0 or p>1:return None
    return p>THRESHOLD

def first_alert(scored_updates,origin,end,aki_stop=None,interrupted=False):
    """First qualifying update in one hospitalization; synthetic policy helper.

    Caller supplies covariate updates and a known stopping episode. Never
    computes AKI from current-stay final diagnoses or future data. A return
    status of NEVER_SCORABLE remains distinct from NO_THRESHOLD_CROSSING.
    """
    if end<origin:raise ValueError('invalid_monitoring_interval')
    if aki_stop is not None and aki_stop<origin:
        return {'status':'AKI_BEFORE_MONITORING','alert_time':None}
    seen=set();valid=0
    for t,p in sorted(scored_updates,key=lambda z:z[0]):
        if t in seen:raise ValueError('group_simultaneous_updates_before_policy')
        seen.add(t)
        if t<origin or t>end:continue
        if aki_stop is not None and t>=aki_stop:break
        hit=above_published_threshold(p)
        if hit is None:continue
        valid+=1
        if hit:return {'status':'ALERTED','alert_time':t}
    status=('POLICY_STOPPED_NO_ALERT' if aki_stop is not None and aki_stop<=end else
            'MONITORING_INTERRUPTED_NO_ALERT' if interrupted else
            'NO_THRESHOLD_CROSSING' if valid else 'NEVER_SCORABLE')
    return {'status':status,'alert_time':None}
