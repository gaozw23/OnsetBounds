"""A1 batched resampling and draw validity; existing core confidence sets."""
import numpy as np
from aligned_geometry import aggregate, NATURAL
from onsetbounds.bootstrap import outer_confidence_set, simultaneous_paired_outer_sets


def patient_weights(n,rng,batch=1):return rng.multinomial(n,np.full(n,1/n),size=batch)


def hospital_weights(hospitals,rng,selected=None,trace=False,groups=None):
    groups={h:np.flatnonzero(hospitals==h) for h in np.unique(hospitals)} if groups is None else groups
    labels=np.array(list(groups)); chosen=rng.choice(labels,len(labels),replace=True) if selected is None else np.asarray(selected)
    counts=np.zeros(len(hospitals),dtype=np.int64); instances=[]
    for label in chosen:
        idx=groups[label]; c=rng.multinomial(len(idx),np.full(len(idx),1/len(idx)))
        counts[idx]+=c
        if trace:instances.append((int(label),c.tolist()))
    return (counts,instances) if trace else counts


def bootstrap_draws(num,den,hospitals,B,rng,scheme,batch_size=16):
    groups={h:np.flatnonzero(hospitals==h) for h in np.unique(hospitals)} if scheme=='hospital' else None
    n=len(num); draws=np.empty((B,18,2)); denominators=np.empty((B,18))
    for start in range(0,B,batch_size):
        m=min(batch_size,B-start)
        weights=patient_weights(n,rng,m) if scheme=='patient' else np.stack([hospital_weights(hospitals,rng,groups=groups) for _ in range(m)])
        draws[start:start+m],denominators[start:start+m]=aggregate(num,den,weights)
    return draws,denominators


def contains(outer,oracle):return bool(np.isfinite(outer).all() and np.isfinite(oracle).all() and outer[0]<=oracle[0] and outer[1]>=oracle[1])


def infer(point, draws, mode='clinical_shared', min_valid=.95):
    if mode not in ('clinical_shared', 'target_specific'):
        raise ValueError(mode)
    valid = np.isfinite(draws).all(axis=2)
    shared = valid.all(axis=1)
    B = len(draws)
    results = []
    for indices in [*[np.array([j]) for j in range(15)], np.arange(15, 18)]:
        mask = shared if mode == 'clinical_shared' else valid[:, indices].all(axis=1)
        count = int(mask.sum())
        own = int(valid[:, indices].all(axis=1).sum())
        finite = np.isfinite(point[indices]).all()
        status = ('OK' if finite and count >= np.ceil(min_valid * B) else
                  'SAMPLE_UNDEFINED' if not finite else 'INSUFFICIENT_VALID_DRAWS')
        q = np.nan
        outer = np.full((len(indices), 2), np.nan)
        if status == 'OK':
            if len(indices) == 1:
                j = indices[0]
                result = outer_confidence_set(*point[j], draws[mask, j, 0],
                                              draws[mask, j, 1], natural_space=NATURAL[j])
                outer[0] = result['outer_lower'], result['outer_upper']
                q = result['critical']
            else:
                outer, q = simultaneous_paired_outer_sets(point[indices], draws[mask][:, indices])
        for j, interval in zip(indices, outer):
            results.append(dict(index=int(j), status=status, B_requested=B,
                                B_valid=count, B_invalid=B-count, B_own_valid=own,
                                excluded_by_other_targets=own-count, critical=q, outer=interval))
    return results
