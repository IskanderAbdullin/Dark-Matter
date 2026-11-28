"""Hard-cut and quasiparticle line-resolution sensitivity of the reference."""
import json
import numpy as np
from scipy.special import roots_legendre

from early_pole_hard_gauge_pair_source_reference import calc, central
from early_pole_decay_mode_clock_reference import intervals,Q


def rate(m,pmin,smin,ct,n):
    z,w=roots_legendre(n)
    p=pmin+(z+1)*(40-pmin)/2
    fp=w*(40-pmin)*p*p/(2*np.pi**2*np.expm1(p))
    all=4*p[:,None]*p[None,:]
    low=np.minimum(1.,smin/all)
    ans=np.zeros_like(all)
    for a,b in zip(z,w,strict=True):
        u=low+.5*(1-low)*(a+1)
        ans+=b*(1-low)*u*central(all*u,m,1.,ct)
    return fp@ans@fp


def one():
    samples=[(.5,1.,.8),(.5,1.,.6),(.75,1.5,.6),(1.,2.,.5)]
    m=.003
    R=[rate(m,*k,80) for k in samples]
    M=[rate(m,*k,48) for k in samples]
    assert abs(R[0]/calc(m,1.,.8,80)-1)<2e-13
    assert max(abs(np.array(R)/M-1))<4e-4
    assert max(np.diff(R))<0.
    ks=[]
    for s in (.003,.004,.0045,.005):
        x=.05; pole=.0008672
        P=np.hypot(x,s/2);B=np.hypot(.06*x,s/2)
        _,lo,hi,_,_=intervals(P,pole,B,Q)
        d=(lo+hi)/2
        p=np.sqrt((d-P)*(d+P))
        res=np.hypot(p+Q,B)-np.hypot(p-Q,B)
        ks.append(dict(e=s,single_line_width_proxy_over_angular_split=s*s*np.log(1/s)/(4*np.pi*res)))
    return dict(cuts=[dict(p_min=k[0],s_min=k[1],cos_cm_max=k[2],
                           unblocked_number_rate_per_e4=r,relative_hard_subset=r/R[0],
                           e_reweight_if_only_normalization_changed=(R[0]/r)**.25)
                      for k,r in zip(samples,R,strict=True)],line_scale_comparison=ks)


if __name__=='__main__':
    print(json.dumps(one(),indent=2,allow_nan=False))
