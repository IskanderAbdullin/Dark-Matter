"""Vacuum neutral-scalar plus photon crossed channel, with charged antiquasiparticle."""
import gzip
import json
from pathlib import Path

import numpy as np
from scipy.special import expit

from early_pole_decay_mode_clock_reference import Q
from early_pole_hard_yukawa_photoproduction_reference import (
    I, G, adj, dot, invariant_form, slash,
)


def com(w,M,m,ms,c):
    r=np.sqrt((w-M-m)*(w+M+m)*(w-M+m)*(w+M-m))/(2*w)
    k=(w-ms)*(w+ms)/(2*w)
    EP=np.hypot(M,r)
    EB=np.hypot(m,r)
    vx=r*np.sqrt((1-c)*(1+c))
    p=np.array([EP,vx,0.,r*c])
    b=np.array([EB,-vx,0.,-r*c])
    photon=np.array([k,0.,0.,k])
    scalar=p+b-photon
    return p,b,photon,scalar


def direct(p,b,k,M,m,pol):
    P,B,K,E=map(slash,(p,b,k,pol))
    # The internal lines are P-K and K-B. The final fermion
    # spin sum of the uncrossed process becomes a v spin sum.
    v=(P-K+M*I)@E/(-2*dot(p,k))+E@(K-B+m*I)/(-2*dot(b,k))
    trace=np.trace((B-m*I)@v@(P+M*I)@adj(v))
    if abs(trace.imag)>1e-10*max(1.,abs(trace.real)):
        raise ArithmeticError('Complex spin sum')
    return float(trace.real)


def thermal(p,b,k,q,T,y):
    # Boost along the outgoing photon for a moving bath observer.
    energy=lambda v:np.cosh(y)*v[0]+np.sinh(y)*v[3]
    ep,eb,eγ,es=map(energy,(p,b,k,q))
    f,g=expit(-ep/T),expit(-eb/T)
    nγ,ns=1/np.expm1(eγ/T),1/np.expm1(es/T)
    forward=f*g*(1+nγ)*(1+ns)
    backward=(1-f)*(1-g)*nγ*ns
    return abs(forward-backward)/max(forward,backward)


def run():
    result=[]
    for grid in ('fine','third'):
        with gzip.open(Path(__file__).parent/'results'/f'early-pole-spectrum-{grid}-m06-2026-11-26.json.gz') as f:
            d=json.load(f)
        M=d['x_end']
        m=M*d['mass_ratio']
        ms=d['checkpoints'][-1]['scalar_pole_over_T']
        W0=M+m
        controls=[]
        ex,ey=np.eye(4)[1:3]
        for h in (1e-6,.01,.08,.5,2.,10.):
            w=W0*(1+h)
            for c in (-.9999,-.7,0.,.7,.9999):
                p,b,k,q=com(w,M,m,ms,c)
                tot=direct(p,b,k,M,m,ex)+direct(p,b,k,M,m,ey)
                if tot<=0.:
                    raise ArithmeticError('Negative crossed spin sum')
                crossed=-invariant_form(-dot(p,k),dot(b,k),-dot(p,b),M,m)
                iw=abs(crossed/tot-1.)
                gauge=abs((direct(p,b,k,M,m,ex+1.1*k)+direct(p,b,k,M,m,ey-.3*k))/tot-1.)
                ward=abs(direct(p,b,k,M,m,k)/tot)
                balance=max(thermal(p,b,k,q,T,y) for T in (.2,1.) for y in (-.8,0.,.9))
                if max(iw,gauge,ward,balance)>1e-8:
                    raise ArithmeticError(('crossing, gauge, or thermal check failed',iw,gauge,ward,balance))
                controls.append(dict(initial_pair_invariant_mass=float(w),cos_theta=c,crossed_spin_sum=tot,
                                     ward_relative_residual=ward,gauge_shift_relative_residual=gauge,
                                     crossed_invariant_relative_error=iw,pointwise_FD_BE_balance_defect=balance))
        hard=[]
        for qmag in (1e-7,Q,1e-5,1e-4,.01):
            E=np.hypot(qmag,ms)
            d0=(W0-ms)*(W0+ms)/2
            vals=[d0/(E-qmag*c) for c in (-1.,0.,1.)]
            # Use exact on-shell s from four-vector reconstruction at both
            # sides of each angular threshold rather than sampling only fit points.
            for c,k0 in zip((-1.,0.,1.),vals,strict=True):
                px=qmag*np.sqrt((1-c)*(1+c))
                qs=np.array([E,px,0.,qmag*c])
                for factor,want in ((.995,False),(1.005,True)):
                    ks=factor*k0*np.array([1.,0.,0.,1.])
                    if (dot(qs+ks,qs+ks)>W0*W0)!=want:
                        raise ArithmeticError('Angular threshold did not bracket')
            hard.append(dict(q=qmag,scalar_E=E,photon_threshold_at_cos_minus_zero_plus=vals,
                             total_energy_threshold=E+vals[0],equilibrium_photon_at_T1_at_minimum=1/np.expm1(vals[0])))
        result.append(dict(grid=grid,M=M,m=m,scalar_mass=ms,minimum_com_photon_energy=(W0-ms)*(W0+ms)/(2*W0),
                           spin_controls=controls,lab_fixed_q_thresholds=hard))
    return dict(scope='One vacuum crossed charge assignment P+anti-b <-> scalar+photon; unaveraged spin and polarization sum, no couplings, integrated density, antiparticle archive assumption, or medium correction.',data=result)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
