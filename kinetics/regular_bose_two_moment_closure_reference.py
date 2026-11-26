"""Maximum-entropy-like regular BE two-moment closure with explicit failure.

Finite-cell temperatures including negative slopes are representation
parameters, never a claim of a physical thermal bath or of condensation.
"""
import json
import numpy as np
from scipy.optimize import brentq
from fermion_two_moment_closure_reference import mesh


def occupation(alpha,beta,s):
    with np.errstate(over="ignore"):
        return 1/np.expm1(alpha+beta*s)


def infer(lo,hi,mass,fill,mean_s,order=96):
    if not (fill>0 and 0<mean_s<1):
        raise ValueError("strictly interior finite BE moments required")
    s,w,_,_=mesh(lo,hi,mass,order)
    base=np.log1p(1/fill)
    target=fill*mean_s
    def inner(beta):
        # Exponent shifted by its lower endpoint: s if beta>=0,
        # or 1-s if beta<0, to avoid cancellation at the pole.
        p=s if beta>=0 else 1-s
        slope=abs(beta)
        def f(a):
            return occupation(a,slope,p)
        if beta==0:
            return base,fill*float(w@s)
        if float(w@f(0.))<=fill:
            return None
        A=brentq(lambda a:float(w@f(a))/fill-1.,0.,base,
                 xtol=np.nextafter(0.,1.),rtol=1e-14)
        return A,float(w@(s*f(A)))
    e0=inner(0.)[1]
    if abs(e0-target)<1e-14*fill:
        return base,0.
    if e0>target:
        sign=1
    else:
        sign=-1
    def reached(x):
        v=inner(x*sign)
        if v is None:
            return None
        return (v[1]-target)*sign
    a,b=0.,1.
    for _ in range(20):
        ub=reached(b)
        if ub is None:
            # Bisect the resolved pole boundary. It may be a genuine
            # continuous-shell capacity or a finite quadrature limitation.
            c,d=a,b
            for _ in range(45):
                mid=(c+d)/2
                if reached(mid) is None:
                    d=mid
                else:
                    c=mid
            b=c
            ub=reached(b)
            if ub is None or ub>0:
                raise ValueError("target outside resolved normal BE family; no condensate inferred")
        if ub<=0:
            break
        a,b=b,2*b
    else:
        raise ValueError("no normal BE shape located within resolved slope range")
    val=brentq(lambda a:reached(a)/fill,a,b,xtol=2e-12)
    beta=val*sign
    A,_=inner(beta)
    return A+max(0.,-beta),beta


def test():
    rows=[]
    for lo,hi,m in [(0.,1e-7,.1),(.0001,.0007,.0004),(0.,.015,.0004),(.1,4.,.04)]:
        s,w,_,_=mesh(lo,hi,m,256)
        for a,b in [(1.,3.),(.01,3.),(.02,.2),(15.,.8),(3.1,-3.),(40.1,-40.),(.05,0.)]:
            f=occupation(a,b,s)
            n,k=float(w@f),float(w@(s*f))
            x,y=infer(lo,hi,m,n,k/n,order=192)
            g=occupation(x,y,s)
            err=max(abs(w@g-n)/n,abs(w@(s*g)-k)/k)
            assert err<2e-8
            rows.append(dict(qlo=lo,qhi=hi,mass=m,alpha=a,beta=b,
                             independent_moment_relative_error=float(err)))
    # A normal positive-temperature near-rest BE shape cannot match an
    # arbitrarily small mean kinetic energy at very high occupation.
    try:
        infer(0.,1e-7,.1,100.,.001)
    except ValueError:
        pass
    else:
        raise AssertionError("impossible normal near-rest shape accepted")
    return dict(cases=rows,max_independent_error=max(z["independent_moment_relative_error"] for z in rows),
                normal_BE_impossibility_control_rejected=True)
if __name__ == "__main__":
    print(json.dumps(test(),indent=2,allow_nan=False))
