"""Conservative fixture: within-shell shape changes preserving all three native moments."""
import argparse
import json

import numpy as np
from scipy.special import roots_legendre
from early_pole_frozen_dynamic_vlasov_reference import load

def shape(edges,occ,mass,sign,active=(1.39,1.48),degree=3):
    if not abs(sign)<=.85 or int(degree) not in (3,4,5):
        raise ValueError("Choose a bounded perturbation and a cubic, quartic or quintic probe")
    z,w=roots_legendre(24)
    L,H=np.asarray(edges[:-1]),np.asarray(edges[1:])
    mid,h=(L+H)/2,(H-L)/2
    anchor=(np.hypot(L,mass)+np.hypot(H,mass))/2
    de=(np.hypot(H,mass)-np.hypot(L,mass))/2
    def basis(i,t):
        p=mid[i]+h[i]*t
        e=np.hypot(p,mass)
        x=(e-anchor[i])/de[i]
        # x² Ec/E represents the exact 1/E constraint after subtracting
        # its constant and linear-in-energy terms, without cancellation.
        return p,e,np.array([np.ones_like(t),x,x*x*anchor[i]/e])
    coefficients={}
    errors=[]
    for i in np.flatnonzero((H>active[0])&(L<active[1])):
        p,e,A=basis(i,z)
        bulk=w*p*p*(1-z*z)
        high=z**degree
        c=np.linalg.solve((A*bulk)@A.T,(A*bulk)@high)
        t=np.linspace(-1,1,1025)
        _,_,D=basis(i,t)
        norm=np.max(np.abs((1-t*t)*(t**degree-c@D)))
        residual=bulk*(high-c@A)/norm
        err=(A@residual)/np.sum(abs(residual))
        assert np.max(abs(err))<2e-12
        errors.append(err.tolist())
        coefficients[i]=(c,norm)
    def evaluate(p):
        k=np.searchsorted(edges,p)-1
        out=np.r_[0.,occ,0.][1+k].copy()
        for i,(a,n) in coefficients.items():
            selected=(k==i)
            if not np.any(selected):
                continue
            t=(p[selected]-mid[i])/h[i]
            _,_,A=basis(i,t)
            out[selected]+=sign*occ[i]*(1-occ[i])*(1-t*t)*(t**degree-a@A)/n
        assert np.min(out)>=0. and np.max(out)<=1.
        return out
    return evaluate,errors

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--degree",type=int,choices=(3,4,5),default=3)
    p.add_argument("--amplitude",type=float,default=.75)
    a=p.parse_args()
    out=[]
    for grid in ('fine','third'):
        d=load(f'results/early-pole-spectrum-{grid}-m06-2026-11-26.json.gz')
        edge=np.asarray(d['endpoint_fermion_momentum_edges'])
        f=np.asarray(d['endpoint_fermion_occupations'])
        p=np.linspace(.45,5.,5000)
        fun,errors=shape(edge,f,.003,a.amplitude,degree=a.degree)
        lo=fun(p)
        fa=f[np.searchsorted(edge,p)-1]
        z,w=roots_legendre(40)
        L,H=edge[:-1],edge[1:]
        keep=(H>1.39)&(L<1.48)
        h=.5*(H[keep]-L[keep])
        q=(.5*(H[keep]+L[keep]))[:,None]+h[:,None]*z
        E=np.hypot(q,.003)
        d=(fun(q.ravel()).reshape(q.shape)-f[keep,None])*(h[:,None]*w)*q*q
        vals=np.stack((d,d*E,d/E),axis=1)
        den=np.maximum(np.sum(abs(vals),axis=2),1e-30)
        independent=float(np.max(abs(vals.sum(axis=2))/den))
        assert independent<2e-11
        out.append(dict(grid=grid,degree=a.degree,amplitude=a.amplitude,max_delta=float(np.max(abs(lo-fa))),
                        max_relative_moment_residual=float(np.max(abs(np.asarray(errors)))),
                        independent_relative_residual=independent))
    print(json.dumps(out,indent=2))
