"""A1 revision DGM v1.1; independent pilot, three policies, early/evaluable mixture."""
from dataclasses import dataclass,asdict
import numpy as np
from scipy.special import expit
VERSION='A1_RAW_ALIGNED_v1.2_RANDOM_HOSPITALS'
@dataclass
class SyntheticData:
    D:np.ndarray
    E:np.ndarray
    T:np.ndarray
    L:np.ndarray
    U:np.ndarray
    alerts:np.ndarray
    hospitals:np.ndarray
    def __len__(self):return len(self.D)

def stream_seed(base,cell,rep,stream):return np.random.SeedSequence([int(base),int(cell),int(rep),int(stream)])
def _structure(cfg,rng,n=None,hospitals=None):
    if cfg['scheme']=='patient':
        size=int(cfg['n'] if n is None else n); return np.zeros(size,dtype=int),np.zeros(size)
    h=int(cfg.get('hospital_count',len(cfg['hospital_sizes'])) if hospitals is None else hospitals)
    if cfg.get('hospital_size_mode')!='iid_discrete_support':raise ValueError('formal hospital DGM requires reviewed iid size mode')
    sizes=rng.choice(np.asarray(cfg['hospital_sizes'],int),size=h,replace=True)
    ids=np.repeat(np.arange(h),sizes); b=rng.normal(0,float(cfg['tau']),h)
    return ids,b[ids]

def calibrate(cfg,seed,pilot_n=5000,pilot_hospitals=120):
    rng=np.random.default_rng(seed); ids,b=_structure(cfg,rng,n=pilot_n,hospitals=pilot_hospitals)
    x=b+rng.normal(size=len(b)); rho=float(cfg['rho']); z=.5*x[:,None]+np.sqrt(rho)*rng.normal(size=(len(b),1))+np.sqrt(1-rho)*rng.normal(size=(len(b),3))
    low,high=-20.,20.; target=float(cfg['event_prevalence'])
    for _ in range(70):
        mid=(low+high)/2
        if expit(mid+b).mean()<target:low=mid
        else:high=mid
    intercept=(low+high)/2; pilot_nonevent=rng.random(len(b))>=expit(intercept+b)
    if not pilot_nonevent.any():raise ValueError('pilot has no non-events')
    thresholds=np.quantile(z[pilot_nonevent],1-float(cfg['non_event_alert_rate']),axis=0,method='linear')
    return {'pilot_non_events':int(pilot_nonevent.sum()),'thresholds':thresholds.tolist(),'event_intercept':(low+high)/2,'pilot_patients':len(b),'pilot_hospitals':len(np.unique(ids)) if cfg['scheme']=='hospital' else None,'source':'independent synthetic pilot; no coverage used','seed_spawn_key':list(seed.entropy) if isinstance(seed,np.random.SeedSequence) else seed}

def generate(cfg,policy,seed,n=None,hospitals=None):
    rng=np.random.default_rng(seed); ids,b=_structure(cfg,rng,n=n,hospitals=hospitals); N=len(ids)
    x=b+rng.normal(size=N); D=rng.random(N)<expit(policy['event_intercept']+b)
    early=rng.random(N)<float(cfg['early_event_fraction']); T=np.where(early,30.,78.)+42*rng.beta(2,2,N)
    width=np.clip(float(cfg['width'])*rng.gamma(4,.25,N)*np.exp(-float(cfg['kappa'])*x-float(cfg['kappa'])**2/2),.5,48)
    relative=rng.beta(1,1,N) if cfg['location']=='uniform' else rng.beta(5,2,N)
    L=np.maximum(24.,T-relative*width); U=np.minimum(168.,T+(1-relative)*width)
    rho=float(cfg['rho']); noise=np.sqrt(rho)*rng.normal(size=(N,1))+np.sqrt(1-rho)*rng.normal(size=(N,3))
    risk=.5*x[:,None]+noise+float(cfg['event_risk_shift'])*D[:,None]
    mu=np.array(cfg['lead_means']); sd=np.array(cfg['lead_sds']); lead=mu+2*x[:,None]+sd*(np.sqrt(rho)*rng.normal(size=(N,1))+np.sqrt(1-rho)*rng.normal(size=(N,3)))
    event_alert=np.maximum(24.,T[:,None]-lead)
    nonevent_alert=24.+144.*rng.beta(2,3,size=(N,3))
    alerts=np.where(D[:,None],event_alert,nonevent_alert)
    available=risk>np.asarray(policy['thresholds'])[None,:]
    available &= (~D[:,None])|(alerts<U[:,None])
    alerts=np.where(available,alerts,np.inf)
    E=D&(L>=30.)
    T=np.where(D,T,np.nan); L=np.where(D,L,np.nan); U=np.where(D,U,np.nan)
    if np.any(D&~((L<=T)&(T<=U)&(L<U))):raise AssertionError('true onset outside generated set')
    return SyntheticData(D,E,T,L,U,alerts,ids)
