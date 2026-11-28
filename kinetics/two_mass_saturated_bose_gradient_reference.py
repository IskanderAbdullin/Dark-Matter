"""Positive evaluation of the line-averaged Bose entropy variables at high occupation."""
import json
from functools import lru_cache

import numpy as np
from scipy.special import roots_legendre

DT=np.longdouble


@lru_cache(maxsize=5)
def nodes(order):
    return tuple(np.asarray(x,dtype=DT) for x in roots_legendre(order))


def bose_average(start,end,capacity,order=16):
    x,y,g=map(DT,(start,end,capacity))
    if min(x,y,g)<=0 or not np.all(np.isfinite([x,y,g])):
        raise ValueError("Positive, finite count endpoints and capacity are required.")
    if x==y:
        return np.log1p(g/x),(g/(g+x))/(2*x)
    v=(y-x)/x
    span=np.log1p(v) if abs(v)<DT(.25) else np.log(y/x)
    u,w=nodes(order)
    cuts=np.linspace(DT(0),span,max(2,1+int(np.ceil(abs(span)/DT(.6)))),dtype=DT)
    Q=DT(0)
    H=DT(0)
    for l,r in zip(cuts[:-1],cuts[1:],strict=True):
        z=l+(r-l)*(1+u)/2
        positive_weights=(r-l)/v*np.exp(z)*w/2
        t=np.expm1(z)/v
        n=x*np.exp(z)
        Q+=positive_weights@np.log1p(g/n)
        H+=positive_weights@(t*(g/(g+n))/n)
    if min(Q,H)<=0:
        raise ValueError("A positive entropy gradient or its curvature was lost.")
    return Q,H


def high_precision(a,b,g):
    import mpmath as mp
    with mp.workdps(100):
        u,v,c=(mp.mpf(n)/d for n,d in (z.as_integer_ratio() for z in (a,b,g)))
        S=lambda x:x*mp.log1p(c/x)+c*mp.log1p(x/c)
        s=S(v)-S(u)
        if u==v:
            return float(mp.log1p(c/u)),float(c/(2*u*(u+c))),0.
        q=s/(v-u)
        h=(q-mp.log1p(c/v))/(v-u)
        return float(q),float(h),float(s)


def run():
    out=[]
    for f,r,g in [(1e8,1.00000000001,1e-12),(1e24,1.,1e-5),(1e23,1.004,1.),
                  (1e22,.99,.06),(1e30,1e8,1e-9),(1e35,1e-8,1e-9),(1e21,1.3,1e-18)]:
        x=DT(f)*DT(g)
        y=x*DT(r)
        gap=DT(g)
        a,b=bose_average(x,y,gap)
        a12,b12=bose_average(x,y,gap,12)
        q,h,d=high_precision(x,y,gap)
        target=np.array([q,h,d]); got=np.array([a,b,(y-x)*a],dtype=float)
        cond=target!=0.
        err=float(np.max(abs(got[cond]/target[cond]-1)))
        loss=abs(a12/a-1)+abs(b12/b-1)
        stable_primitive=lambda z:z*np.log1p(gap/z)+gap*np.log1p(z/gap)
        raw=np.log(x+gap)-np.log(x)
        if err>5e-13 or loss>5e-13:
            raise ValueError(f"Stable high-occupation check failed at {f}, {r}: {err}, {loss}.")
        out.append(dict(starting_bose_occupation=f,end_to_start_count=r,chosen_capacity=g,
                        averaged_affinity=float(a),negative_end_gradient=float(b),
                        stable_entropy_increment=float((y-x)*a),
                        direct_float_difference_increment=float(stable_primitive(y)-stable_primitive(x)),
                        naive_start_log_difference=float(raw),high_precision_relative_error=err,
                        relative_order_error=float(loss)))
    return dict(scope="Local arithmetic for high Bose occupations; not a full saturated kinetic solve.",
                cases=out)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
