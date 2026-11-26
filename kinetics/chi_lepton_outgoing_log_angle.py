"""Exact change of variables for finite-vector t-channel outgoing CM polar angle.

Let z=1-cos(theta), A=2*pstar²/mV². Vector propagator² is
(1+A*z)^-2. y=log(1+A*z) gives (1+A*z)^-2 dz=exp(-y)dy/A.
This precisely clusters angular nodes near the physically relevant
forward region. This is not an analytic solution of the full collision.
"""
from __future__ import annotations
import math
import json
import numpy as np
from numpy.polynomial.legendre import leggauss
from kinetics.independent_chi_lepton_5d_quadrature import evaluate


def exact_propagator_moments(A):
    """I0=∫ dc/[1+A(1-c)]²; I1=∫(1-c)dc/[1+A(1-c)]²."""
    if not(math.isfinite(A) and A>0):raise ValueError("physical A>0")
    i0=2/(1+2*A)
    if A<1e-4:
        i1=2-(16/3)*A+12*A*A-(128/5)*A**3
    else:
        i1=(math.log1p(2*A)-2*A/(1+2*A))/(A*A)
    return i0,i1


def cm_angle_moments(A,n,scheme):
    if n<2:raise ValueError
    x,w=leggauss(n)
    if scheme=="uniform_cos":
        c=x
        W=w/(1+A*(1-x))**2
    elif scheme=="propagator_compensated":
        c=1-(1-x)/(1+A*(1+x))
        W=w/(1+2*A)
    elif scheme=="log_transfer":
        ymax=math.log1p(2*A)
        y=(x+1)*ymax/2
        c=1-np.expm1(y)/A
        W=w*(ymax/(2*A))*np.exp(-y)
    else:raise ValueError
    return float(sum(W)),float(np.dot(W,1-c))


def audit(full=False):
    rows=[]
    for A in (.1,1,10,100,1000,10000):
        e0,e1=exact_propagator_moments(A)
        row={"A":A,"analytic_I0":e0,"analytic_I1":e1}
        for scheme in ("uniform_cos","propagator_compensated","log_transfer"):
            for n in ((8,16) if not full else (4,8,16,32)):
                a,b=cm_angle_moments(A,n,scheme)
                row[f"{scheme}_n{n}"]={"I0_relative_error":a/e0-1,
                    "I1_relative_error":b/e1-1}
        rows.append(row)
    physical=[]
    for mv,shape in ((.1,"oscillatory"),(.01,"smooth")):
        for scheme in ("uniform_cos","log_transfer"):
            for nc in ((8,16,32) if not full else (8,16,32,64)):
                kw=dict(radial_scheme="segmented_p",segments=12 if mv==.1 else 8,
                        nr=5 if mv==.1 else 4,nmu=8 if mv==.1 else 6,
                        nc=nc,nphi=16 if mv==.1 else 12,mv=mv,
                        angular_scheme=scheme)
                q=evaluate(shape,**kw)
                physical.append({"mv":mv,"shape":shape,"scheme":scheme,
                    "nc":nc,"Qchi":q["Qchi"],"entropy":q["entropy"]})
    return {"schema_version":1,"analytic_propagator":rows,
            "synthetic_5D_control":physical,
            "complete_physical_chi_gamma_or_86eV_dm":False}


if __name__=="__main__":
    import argparse
    p=argparse.ArgumentParser()
    p.add_argument("--full",action="store_true")
    a=p.parse_args()
    print(json.dumps(audit(full=a.full),indent=2))
