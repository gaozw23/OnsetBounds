"""A1 independent-unit ratio moments and oracle uncertainty classification."""
import numpy as np
from scipy.stats import norm


def empty_state():return {'units':0,'patients':0,'batches':0,'sum_y':np.zeros((18,2)),'sum_d':np.zeros(18),'sum_yy':np.zeros((18,2)),'sum_yd':np.zeros((18,2)),'sum_dd':np.zeros(18),'events':0,'evaluable':0,'alerts':np.zeros(3),'size_min':None,'size_max':None}


def add_units(s,y,d):
 s['units']+=len(y);s['sum_y']+=y.sum(axis=0);s['sum_d']+=d.sum(axis=0);s['sum_yy']+=(y*y).sum(axis=0);s['sum_yd']+=(y*d[:,:,None]).sum(axis=0);s['sum_dd']+=(d*d).sum(axis=0)


def estimates(s):
 n=s['units'];ds=s['sum_d'];theta=np.divide(s['sum_y'],ds[:,None],out=np.full((18,2),np.nan),where=ds[:,None]>0)
 ss=np.maximum(0,s['sum_yy']-2*theta*s['sum_yd']+theta**2*s['sum_dd'][:,None]);se=np.divide(np.sqrt(n/(n-1)*ss),ds[:,None],out=np.full((18,2),np.nan),where=ds[:,None]>0) if n>1 else np.full((18,2),np.nan)
 return theta,se


def oracle_band_classification(outer,theta,se):
 z=float(norm.ppf(1-.05/(2*36)));ol=np.asarray(outer)[:,0];ou=np.asarray(outer)[:,1];tl=np.asarray(theta)[:,0];tu=np.asarray(theta)[:,1];sl=np.asarray(se)[:,0];su=np.asarray(se)[:,1]
 nominal=(ol<=tl)&(ou>=tu);definite=(ol<=tl-z*sl)&(ou>=tu+z*su);possible=(ol<=tl+z*sl)&(ou>=tu-z*su)
 return {'nominal':nominal,'definite':definite,'possible':possible,'boundary_sensitive':definite!=possible,'z':z}
