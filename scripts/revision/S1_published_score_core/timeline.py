"""Standalone retrospective timeline core. No clinical IO or cohort results."""
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
import json
import numpy as np
import pandas as pd
import ugwuowo_2020_model as FORMULA
HOUR=3_600_000_000_000
DAY=24*HOUR
RULES=json.loads(Path(__file__).with_name('input_mapping.json').read_text(encoding='utf-8'))
LABS=list(RULES['lab_items'])
def hr(x):return int(round(x*HOUR))

def race_code(value):
    s=str(value).strip().upper()
    if s in {'BLACK/AFRICAN AMERICAN','BLACK/CAPE VERDEAN','BLACK/CARIBBEAN ISLAND','BLACK/HAITIAN'}:return 1.
    if s=='WHITE' or s.startswith('WHITE - ') or s=='ASIAN' or s.startswith('ASIAN - '):return 0.
    if s in {'AMERICAN INDIAN/ALASKA NATIVE','NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER'}:return 0.
    return np.nan

def service_code(s):
    if s in RULES['surgical_1']:return 1.
    if s in RULES['surgical_0']:return 0.
    return np.nan

def icu_code(s):
    if s in RULES['icu_1']:return 1.
    if s in RULES['icu_0']:return 0.
    return np.nan

@dataclass(frozen=True)
class Lab:
    variable:str
    obs:int
    avail:int
    value:float
    current:bool=True

def active_location(intervals,t):
    flags=[x for start,end,x in intervals if start<=t<end]
    return flags[0] if flags and all(np.isfinite(x) for x in flags) and len(set(flags))==1 else np.nan

def active_service(changes,t):
    prior=[(s,x) for s,x in changes if s<=t]
    if not prior:return np.nan
    last=max(s for s,x in prior);flags=[x for s,x in prior if s==last]
    return flags[0] if flags and all(np.isfinite(x) for x in flags) and len(set(flags))==1 else np.nan

def medians(records,t=None,current=None):
    d=defaultdict(list)
    for r in records:
        if t is not None and (r.avail>t or r.obs>t):continue
        if current is not None and r.current!=current:continue
        d[(r.variable,r.obs)].append(r.value)
    return {k:float(np.median(v)) for k,v in d.items()}

def latest_inputs(known,t):
    m=medians(known,t,current=True);out={};provenance=[]
    for v in LABS:
        obs=[s for vv,s in m if vv==v]
        if not obs:continue
        s=max(obs);out[v]=m[(v,s)]
        provenance.extend((r.obs,r.avail) for r in known if r.current and r.variable==v and r.obs==s and r.avail<=t)
    return out,provenance

def detect_any_aki(known,t,admission):
    """Record-based stop, recomputed as-of now; UNKNOWN without reference."""
    m=medians([r for r in known if r.variable=='creatinine'],t)
    seq=sorted((s,x) for (v,s),x in m.items())
    current_obs={r.obs for r in known if r.current and r.variable=='creatinine' and r.avail<=t and r.obs<=t}
    return detect_from_sequence(seq,current_obs,admission)

def detect_from_sequence(seq,current_obs,admission):
    unknown=0;assessed=0
    for s,c in seq:
        if s<admission or s not in current_obs:continue
        a=[x for q,x in seq if s-2*DAY<=q<s]
        b=[x for q,x in seq if s-7*DAY<=q<s]
        if not a and not b:unknown+=1;continue
        assessed+=1
        # Decimal source values are compared with a tiny floating representation
        # tolerance ONLY for the .3 arithmetic subtraction (not risk threshold).
        if (a and c-min(a)>=.3-1e-12) or (b and c>=1.5*min(b)-1e-12):
            return True,unknown,assessed
    return False,unknown,assessed

def baseline(records,icu_start,t0):
    m=medians([r for r in records if r.variable=='creatinine' and r.current and icu_start<=r.obs<=t0],t0)
    return min(m.values()) if m else np.nan

def outcome(records,c0,t0,H):
    """Observation-clock measurements. No future U is passed to the policy."""
    if not np.isfinite(c0):return dict(D=None,L=None,U=None,E=False,finite=False,n_post=0)
    m=medians([r for r in records if r.current and r.variable=='creatinine' and t0<=r.obs<=H])
    seq=sorted((s,c) for (v,s),c in m.items())
    post=[(s,c) for s,c in seq if s>t0]
    positive=[s for s,c in post if c>=2*c0]
    if positive:
        u=min(positive);neg=[s for s,c in seq if t0<=s<u and c<2*c0]
        l=max(neg) if neg else None;finite=l is not None and t0<=l<u<=H
        return dict(D=1,L=l,U=u,E=bool(finite and l>=t0+6*HOUR),finite=bool(finite),n_post=len(post))
    return dict(D=0 if post else None,L=None,U=None,E=False,finite=False,n_post=len(post))

def run_policy(records,intervals,services,static,admission,t0,H,known_eskd=False,death=None):
    """Replay hospital episode; evaluate stop before risk in each update batch.

    No fixed-outcome U, D, C0 or event status enters this function.
    """
    records=sorted(records,key=lambda r:(r.avail,r.obs,r.variable,r.value))
    batch=defaultdict(list)
    known=[]
    for r in records:
        if r.avail<admission:known.append(r)
        elif admission<=r.avail<=H:batch[r.avail].append(r)
    decisions={admission,t0}
    decisions.update(s for s,x in services if admission<=s<=H)
    for s,e,x in intervals:
        if admission<=s<=H:decisions.add(s)
        if admission<=e<=H:decisions.add(e)
    decisions.update(t for t,rs in batch.items() if any(r.current for r in rs))
    times=decisions|set(batch)
    times=sorted(t for t in times if admission<=t<=H)
    values=[];valid_times=[];missing_counts=defaultdict(int)
    unknown_stop_checks=0;stopped=None;initial_cr=None;initial_gt4=False
    timeline_violations=0;stopped_score_attempts=0;stop_priority_batches=0
    pre_requests=post_requests=pre_missing=post_missing=0
    before_origin_scorable=False;landmark_complete=False;landmark_labs_complete=False;landmark_icu=landmark_service=np.nan
    first_complete=None;protocol_stopped_updates=0;esrd_blocked=0;stopped_pre=stopped_post=0
    current_values=defaultdict(list);current_medians={};current_maxavail={};last_obs={}
    cr_values=defaultdict(list);cr_medians={};current_cr_obs=set()
    def insert(rs):
        touched=set();cr_touched=set()
        for r in rs:
            if r.current:
                k=(r.variable,r.obs);current_values[k].append(r.value);touched.add(k)
                current_maxavail[k]=max(current_maxavail.get(k,r.avail),r.avail)
                last_obs[r.variable]=max(last_obs.get(r.variable,r.obs),r.obs)
            if r.variable=='creatinine':
                cr_values[r.obs].append(r.value);cr_touched.add(r.obs)
                if r.current:current_cr_obs.add(r.obs)
        for k in touched:current_medians[k]=float(np.median(current_values[k]))
        for s in cr_touched:cr_medians[s]=float(np.median(cr_values[s]))
    insert(known)
    for t in times:
        group=batch.get(t,[]);known.extend(group)
        insert(group)
        current_new=[r for r in group if r.current and r.variable=='creatinine']
        if initial_cr is None and current_new:
            first_obs=min(r.obs for r in current_new)
            initial_cr=float(np.median([r.value for r in current_new if r.obs==first_obs]))
            initial_gt4=initial_cr>4.
        if stopped is None and any(r.variable=='creatinine' for r in group):
            positive,unknown,assessed=detect_from_sequence(sorted(cr_medians.items()),current_cr_obs,admission)
            unknown_stop_checks+=unknown
            if positive:stopped=t;stop_priority_batches+=1
        icu=active_location(intervals,t);surg=active_service(services,t)
        labs={v:current_medians[(v,s)] for v,s in last_obs.items()}
        prov=[(s,current_maxavail[(v,s)]) for v,s in last_obs.items()]
        complete=(len(labs)==7 and all(np.isfinite(static.get(v,np.nan)) for v in ['age','male','black']) and np.isfinite(icu) and np.isfinite(surg))
        if t==t0:
            landmark_complete=complete;landmark_labs_complete=len(labs)==7;landmark_icu=icu;landmark_service=surg
        if t not in decisions or t>=H or (death is not None and t>=death):continue
        # H is a closed outcome boundary; a transfer-out/end state may itself
        # be unknown, so endpoint update cannot manufacture an ICU score.
        is_pre=t<t0
        pre_requests+=int(is_pre);post_requests+=int(not is_pre)
        if stopped is not None and t>=stopped:
            protocol_stopped_updates+=1;stopped_pre+=int(is_pre);stopped_post+=int(not is_pre);continue
        if known_eskd:esrd_blocked+=1;continue
        if initial_gt4:continue
        if not complete:
            pre_missing+=int(is_pre);post_missing+=int(not is_pre)
            for v in LABS:
                if v not in labs:missing_counts[v]+=1
            for v in ['age','male','black']:
                if not np.isfinite(static.get(v,np.nan)):missing_counts[v]+=1
            if not np.isfinite(icu):missing_counts['icu']+=1
            if not np.isfinite(surg):missing_counts['surgical']+=1
            continue
        timeline_violations+=sum(o>t or a>t for o,a in prov)
        row={**labs,**static,'icu':float(icu),'surgical':float(surg)}
        values.append(row);valid_times.append(t)
        if first_complete is None:first_complete=t
        before_origin_scorable |= is_pre
    if values:
        f=pd.DataFrame(values,columns=list(FORMULA.COEFFICIENTS))
        scored=FORMULA.score(f,FORMULA.UNITS)
        if not scored.status.eq('SCORED_FIXED_FORMULA').all():raise RuntimeError('FIXED_FORMULA_INPUT_OR_NUMERICAL_FAILURE')
        risks=scored.risk.to_numpy()
        if not np.isfinite(risks).all() or not ((risks>=0)&(risks<=1)).all():raise RuntimeError('INVALID_PROBABILITY')
        alerts=np.flatnonzero(risks>.15)
        alert=valid_times[int(alerts[0])] if len(alerts) else None
    else:risks=np.array([],dtype=float);alert=None
    assert all(stopped is None or t<stopped for t in valid_times),'SCORED_AFTER_ANY_AKI_STOP'
    return dict(alert=alert,stop=stopped,initial_cr=initial_cr,initial_gt4=initial_gt4,
      times=valid_times,risks=risks,requests_pre=pre_requests,requests_post=post_requests,
      missing_pre=pre_missing,missing_post=post_missing,missing_counts=dict(missing_counts),
      unknown_detector_checks=unknown_stop_checks,protocol_stopped_updates=protocol_stopped_updates,
      stopped_pre=stopped_pre,stopped_post=stopped_post,
      stop_priority_batches=stop_priority_batches,timeline_violations=timeline_violations,
      first_complete=first_complete,pre_scorable=before_origin_scorable,
      landmark_complete=landmark_complete,landmark_labs_complete=landmark_labs_complete,landmark_icu=landmark_icu,landmark_service=landmark_service,
      ever_scorable=bool(values),nonfinite_scores=0,out_of_range_probabilities=0)

def sensitivity(event,alert):
    if not event['E']:return None
    l,u=event['L'],event['U']
    if alert is None:return (0,0,0,0)
    a=alert;lo=int(l>=a+6*HOUR and u<=a+24*HOUR)
    hi=int(max(l,a+6*HOUR)<=min(u,a+24*HOUR))
    at_u=int(6*HOUR<=u-a<=24*HOUR)
    # Integer doubled coordinates avoid midpoint rounding at a boundary.
    mid=int(12*HOUR<=l+u-2*a<=48*HOUR)
    assert lo<=at_u<=hi and lo<=mid<=hi
    return lo,hi,at_u,mid
