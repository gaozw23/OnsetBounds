"""A1 contribution adapter to the existing OnsetBounds statistical core."""
import numpy as np
from onsetbounds.bounds import single_policy_bounds
from onsetbounds.paired import paired_valid_bounds
from onsetbounds.utility import utility_bounds

PAIRS = ((0, 1), (0, 2), (1, 2))
POLICIES = ('A', 'B', 'C')
METRICS = ('valid', 'premature', 'late_missed', 'ppv', 'utility')
NAMES = tuple(f'{p}:{m}' for p in POLICIES for m in METRICS) + ('A-B', 'A-C', 'B-C')
NATURAL = np.array([(-1., 1.) if n.endswith('utility') or n in NAMES[15:]
                    else (0., 1.) for n in NAMES])


def geometry(lower, upper, alerts):
    """Closed 6--24 h bounds; utility retains the submission definition."""
    lower, upper, alerts = np.asarray(lower), np.asarray(upper), np.asarray(alerts)
    n = len(lower)
    if alerts.shape != (n, 3) or upper.shape != (n,):
        raise ValueError('expected aligned intervals and three policies')
    result = {name: np.empty((n, 3, 2)) for name in
              ('valid', 'premature', 'late_missed', 'utility', 'paired')}
    for i, (l, u, row) in enumerate(zip(lower, upper, alerts)):
        for j, alert in enumerate(row):
            b = single_policy_bounds(l, u, alert)
            result['valid'][i, j] = b.valid_lower, b.valid_upper
            result['premature'][i, j] = b.premature_lower, b.premature_upper
            result['late_missed'][i, j] = b.late_missed_lower, b.late_missed_upper
            result['utility'][i, j] = utility_bounds(l, u, alert)
        for j, (a, b) in enumerate(PAIRS):
            pair = paired_valid_bounds(l, u, row[a], row[b])
            result['paired'][i, j] = pair.lower, pair.upper
    return result


def contributions(data):
    n=len(data.D); num=np.zeros((n,18,2)); den=np.zeros((n,18)); event=np.flatnonzero(data.D)
    g=geometry(data.L[event],data.U[event],data.alerts[event]) if len(event) else None
    for p in range(3):
        alerted=np.isfinite(data.alerts[:,p]).astype(float)
        for m,metric in enumerate(METRICS):
            j=p*5+m
            if metric in ('valid','premature','late_missed'):
                den[:,j]=data.E
                if len(event):num[event,j]=g[metric][:,p]*data.E[event,None]
            elif metric=='ppv':
                den[:,j]=alerted
                if len(event):num[event,j]=g['valid'][:,p]
            else:
                den[:,j]=1.; num[:,j]=-alerted[:,None]
                if len(event):num[event,j]=g['utility'][:,p]
    for j in range(3):
        den[:,15+j]=data.E
        if len(event):num[event,15+j]=g['paired'][:,j]*data.E[event,None]
    return num,den


def aggregate(num,den,weights):
    w=np.atleast_2d(np.asarray(weights,float)); sums=(w@num.reshape(len(num),-1)).reshape(len(w),18,2); ds=w@den
    estimates=np.divide(sums,ds[:,:,None],out=np.full(sums.shape,np.nan),where=ds[:,:,None]>0)
    return estimates,ds
