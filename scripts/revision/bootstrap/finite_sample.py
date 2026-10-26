"""A1 finite-sample conditional-binomial sensitivity for iid patient targets.

Use only for binary contribution sums with their exact integer denominators;
the conditional-binomial guard does not cover hospital-resampled targets.
"""
import math
from fractions import Fraction
from functools import lru_cache
import numpy as np
from scipy.special import betaincinv


def integer(value,name):
    if isinstance(value,(bool,np.bool_)) or not isinstance(value,(int,np.integer)):
        raise ValueError(name+' must be an integer, not a rounded proportion')
    return int(value)


@lru_cache(maxsize=None,typed=True)
def rounding_cell_count(x,n):
    """Certify a UNIQUE integer whose correctly rounded binary64 quotient is x.

    Uses exact rational midpoints between adjacent binary64 values, not an
    arbitrary tolerance on a rounded decimal proportion. Requires provenance
    establishing a binary sum divided by the exact integer denominator.
    Returns (count, exact-rounding-cell width after multiplication by n).
    """
    n=integer(n,'n')
    if n<=0 or n>2**53 or not math.isfinite(x) or not 0<=x<=1:
        raise ValueError('nonpositive/unsupported denominator or invalid probability')
    x=float(x);q=Fraction.from_float(x)
    left=(Fraction.from_float(math.nextafter(x,-math.inf))+q)*n/2
    right=(q+Fraction.from_float(math.nextafter(x,math.inf)))*n/2
    low=max(0,-((-left.numerator)//left.denominator))
    high=min(n,right.numerator//right.denominator)
    # Including ties conservatively can only make certification harder.
    if low!=high or float(low/n)!=x:
        raise ValueError('success count not uniquely certified from full-precision quotient')
    return low,float(right-left)


@lru_cache(maxsize=None,typed=True)
def cp_limits(s,n,delta=.025):
    """One-sided lower/upper each with tail error delta; n=0 is vacuous."""
    s=integer(s,'s');n=integer(n,'n')
    if n<0 or s<0 or s>n or not 0<delta<.5:
        raise ValueError('invalid count/denominator/tail probability')
    if n==0:return (0.,1.)
    lo=0. if s==0 else (math.exp(math.log(delta)/n) if s==n else float(betaincinv(s,n-s+1,delta)))
    hi=1. if s==n else (-math.expm1(math.log(delta)/n) if s==0 else float(betaincinv(s+1,n-s,1-delta)))
    return lo,hi


def exact_outer(s_minus,s_plus,n,alpha=.05):
    sm=integer(s_minus,'S_minus');sp=integer(s_plus,'S_plus');n=integer(n,'n')
    if sm>sp:raise ValueError('lower successes exceed upper successes')
    return np.array([cp_limits(sm,n,alpha/2)[0],cp_limits(sp,n,alpha/2)[1]])


def guard_outer(bootstrap,exact):
    b=np.asarray(bootstrap,float);e=np.asarray(exact,float)
    if b.shape!=(2,) or e.shape!=(2,) or not np.isfinite(e).all() or e[0]>e[1]:
        raise ValueError('invalid interval')
    if not np.isfinite(b).all():return e.copy()
    if b[0]>b[1]:raise ValueError('reversed bootstrap interval')
    return np.array([min(b[0],e[0]),max(b[1],e[1])])
