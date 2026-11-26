"""Recheck finite-q native CLASS-like input moments and dilation response."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline


def audit(path,ns,qmax):
    q,F=np.loadtxt(path).T
    assert np.all(F>0.)
    f=CubicSpline(q,F)
    hi=-((F[-1]-F[-2])/(q[-1]-q[-2]))/F[-1]
    assert hi>0.
    def evaluate(p):
        x=np.clip(p,q[0],q[-1])
        dst=f(x)
        derivative=f(x,1)
        tail=F[-1]*np.exp(-hi*np.maximum(0.,p-q[-1]))
        return np.where(p>q[-1],tail,dst),np.where(p<q[0],0.,np.where(p>q[-1],-hi*tail,derivative))
    def moment(p,t):
        z,w=evaluate(p)
        n=[]; r=[]
        for k in (0,1,2,3,4):
            value=t@(p**(2+k)*z)
            dot=t@(p**(3+k)*w)
            n.append(value)
            r.append(dot/(-(3+k)*value)-1)
        return np.array(n),np.array(r)
    za=np.linspace(0,qmax,64001); wa=np.ones(len(za))*qmax/(len(za)-1);wa[[0,-1]]*=.5
    ref,dref=moment(za,wa)
    assert np.max(abs(dref))<2e-6
    rows=[]
    for n in ns:
        p=np.arange(1,n+1)*qmax/n
        w=np.full(n,qmax/n);w[-1]*=.5
        med,err=moment(p,w)
        rows.append(dict(momentum_bins=n,physical_normalization_residual=med[0]/ref[0]-1.,
                  normalized_q_moment_errors=(med[1:]/med[0]/(ref[1:]/ref[0])-1.).tolist(),
                  dilation_integral_by_parts_residual=err.tolist()))
    return dict(file=Path(path).name,qmax=qmax,results=rows,
       convention="Native scheme-3 trapezoidal nodes and weights, with analytic derivative of an independent cubic spline and matching external-tails convention. Dilation identity residual is a quadrature diagnostic, not a CLASS transfer error or proof.")

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("standard");p.add_argument("fine")
    a=p.parse_args()
    print(json.dumps(dict(standard=audit(a.standard,[64,128,256,512,1040,4160],20.),
                         fine=audit(a.fine,[64,128,256,512,1040,4160],20.)),
                     indent=2,allow_nan=False))
