"""Independently integrate the initial bare source against archived scalar shells."""
import argparse
import gzip
import json

import numpy as np
from scipy.special import roots_legendre

from early_pole_native_aperture_witness_reference import GRID, PROFILES, refine, mode
from early_pole_decay_mode_clock_reference import G, KAPPA


def case(profile):
    grid=profile['grid']
    with gzip.open(GRID/f'early-pole-spectrum-{grid}-m06-2026-11-26.json.gz','rt') as f:
        a=json.load(f)
    M=a['x_end']
    b=a['mass_ratio']*M
    m=a['checkpoints'][-1]['scalar_pole_over_T']
    P0=np.asarray(a['endpoint_parent_momentum_edges'])
    F=np.asarray(a['endpoint_parent_occupations'])
    pe,base,new,_=refine(P0,F,profile['parent'],M)
    DE=np.asarray(a['endpoint_fermion_momentum_edges'])
    DF=np.asarray(a['endpoint_fermion_occupations'])
    old=np.asarray(a['scalar_momentum_edges'])**3
    counts=np.asarray(a['endpoint_scalar_signed_shell_number_over_T3'])
    C=G*KAPPA*M*((M+b)**2-m*m)
    A=np.array([0.,.85,.97])
    lo,hi=1e-18,1e-15

    def compute(q,N,V):
        out=np.zeros(len(A))
        for x,n,v in zip(q,N,V,strict=True):
            E,I0,_=mode(M,m,b,x,pe,base,new,DE,DF,0.)
            _,I1,_=mode(M,m,b,x,pe,base,new,DE,DF,1.)
            I=I0+(I1-I0)*A[:,None]
            out+=C*(n*(I[:,0]-I[:,1])+v*I[:,2])/(x*E)
        return out

    def regrid(n):
        edges=np.geomspace(lo,hi,n+1)
        ov=np.maximum(0.,np.minimum(old[1:,None],edges[None,1:])-
                         np.maximum(old[:-1,None],edges[None,:-1]))
        N=(counts/np.diff(old))@ov
        V=np.diff(edges)/(6*np.pi**2)
        q=np.cbrt((edges[1:]+edges[:-1])/2)
        return compute(q,N,V)

    def native(order):
        # Resolved independent scalar integration in q cubed, splitting
        # at original native density discontinuities, unlike regrid above.
        e=np.unique(np.r_[np.geomspace(lo,hi,65),old[(lo<old)&(old<hi)]])
        h=np.diff(e)[:,None]/2
        z,w=roots_legendre(order)
        x=((e[:-1,None]+e[1:,None])/2+h*z).ravel()
        u=(h*w).ravel()
        density=counts/np.diff(old)
        return compute(np.cbrt(x),u*density[np.searchsorted(old,x)-1],u/(6*np.pi**2))

    ns=(4,8,16,32,64,128,256,512)
    rates=np.array([regrid(n) for n in ns])
    r24, r64=native(16),native(32)
    assert np.max(abs(r64-r24))<2e-5
    assert np.max(abs(r64-rates[-1]))<2e-5
    assert np.all(rates[2:,0]>0.) and np.all(rates[2:,1:]<0.)
    return dict(grid=grid,witness_fractions=A.tolist(),geometric_bins=list(ns),representative_rates=rates.tolist(),
                panel_gauss16=r24.tolist(),panel_gauss32=r64.tolist())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--profiles',type=str,default=str(PROFILES))
    a=p.parse_args()
    with open(a.profiles) as f:
        profiles=json.load(f)
    print(json.dumps([case(row) for row in profiles],indent=2,allow_nan=False))
