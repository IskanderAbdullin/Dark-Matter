"""Adaptive parent-momentum integral on separately located daughter packet bounds."""
import argparse
import gzip
import json
from math import fsum

import numpy as np
from scipy.integrate import quad
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq
from scipy.special import expit

from early_pole_decay_mode_clock_reference import G,KAPPA
from early_pole_daughter_local_fraction_curvature_attribution import initial
from early_pole_frozen_dynamic_vlasov_reference import load


def make_interpolant(grid):
    d=load(f'results/early-pole-spectrum-{grid}-m06-2026-11-26.json.gz')
    p=np.asarray(d['endpoint_fermion_momentum_edges'])
    f=np.asarray(d['endpoint_fermion_occupations'])
    v=p**3/3.
    h=PchipInterpolator(v,np.r_[0.,np.cumsum(f*np.diff(v))],extrapolate=False).derivative()
    def gap(x):
        j=np.searchsorted(p,x)-1
        if not 0<=j<len(f):
            raise ValueError('Daughter outside native archive.')
        return max(float(h(x**3/3.)),0.)-f[j]
    return p,gap


def bounds(A,B,va,vb,E,q,scan):
    # Each dressed energy is affine in its underlying square-root bare energy.
    ap,bp=va[0]-va[1]*A.ec/A.de,va[1]/A.de
    ad,bd=vb[0]-vb[1]*B.ec/B.de,vb[1]/B.de
    lo=A.lo
    hi=A.hi
    # Intersect the parent momentum interval with the exact daughter energy range.
    wl=ad+bd*np.hypot(B.lo,B.m)+E
    wh=ad+bd*np.hypot(B.hi,B.m)+E
    low=(wl-ap)/bp
    high=(wh-ap)/bp
    if high<=A.m:
        return (),None
    to_p=lambda e:np.sqrt(max(0.,(e-A.m)*(e+A.m)))
    lo=max(lo,to_p(low))
    hi=min(hi,to_p(high))
    if hi<=lo:
        return (),None
    def geom(p):
        w=ap+bp*np.hypot(p,A.m)
        wb=w-E
        eb=(wb-ad)/bd
        k=np.sqrt(max(0.,(eb-B.m)*(eb+B.m)))
        return k,w,wb,bd*k/eb
    f=lambda p:geom(p)[0]-p
    grid=np.linspace(lo,hi,scan)
    cuts=[lo,hi]
    for side in (-q,q):
        g=lambda w:f(w)-side
        v=np.array([g(w) for w in grid])
        for j in range(scan-1):
            if v[j]==0.:
                cuts.append(grid[j])
            elif v[j]*v[j+1]<0.:
                cuts.append(brentq(g,grid[j],grid[j+1],xtol=8e-16))
    cut=np.unique(cuts)
    return [(a,b) for a,b in zip(cut[:-1],cut[1:],strict=True) if abs(f((a+b)/2))<q],geom


def compare(bundle,scan):
    C,n,(x,_,_),_,depth=initial(bundle)
    M=C.M
    with gzip.open(bundle,'rt') as f:
        grid=json.load(f)['branches']['steps']['grid']
    native_edges,gap=make_interpolant(grid)
    e,F=M.weights(x,True)[:2]
    ee,_,_,fs,m=C.scalar(n)
    w=.05*np.exp(getattr(M,'eta',0.))
    scale=C.sw*G*KAPPA*w*((1.06*w)**2-m*m)/(C.q*ee)
    daughter=[j for j in range(M.r[2].start,M.r[2].stop)
              if M.cells[j].lo>=1.4138888888888888-2e-14 and M.cells[j].hi<=1.4666666666666666+2e-14]
    count=0
    mode=[]
    maxerr=0.
    for q,ek,fc,wk in zip(C.q,ee,fs,scale,strict=True):
        sums=[]
        for i in range(M.r[1].start,M.r[1].stop):
            A=M.cells[i]
            for j in daughter:
                B=M.cells[j]
                intervals,geom=bounds(A,B,e[i],e[j],ek,q,scan)
                for a,b in intervals:
                    count+=1
                    ka,kb=[geom(p)[0] for p in (a,b)]
                    cuts=[brentq(lambda t:geom(t)[0]-p,a,b,xtol=8e-16)
                          for p in native_edges[(native_edges>ka)&(native_edges<kb)]]
                    def rhs(p):
                        k,wp,wb,d=geom(p)
                        fp=expit(F[i,0]+F[i,1]*A.basis(p))
                        return -wk*p*k/(wp*wb*d)*(fp+fc)*gap(k)
                    v,err=quad(rhs,a,b,points=cuts,epsabs=2e-14,epsrel=1e-11)
                    sums.append(v)
                    maxerr=max(maxerr,err)
        mode.append(fsum(sums))
    return dict(hard_depth=depth,cold_modes=len(n),adaptive_direct_rate=fsum(mode),adaptive_direct_by_mode=mode,
                independently_located_packet_count=count,root_scan=scan,max_quad_reported_error=maxerr)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('archives',nargs='+')
    p.add_argument('--root-scan',type=int,default=129)
    v=p.parse_args()
    if v.root_scan<33:
        p.error('Require at least 33 independent root scan points.')
    print(json.dumps([compare(f,v.root_scan) for f in v.archives],indent=2,allow_nan=False))
