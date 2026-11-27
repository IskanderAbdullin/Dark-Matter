"""Certify cap-free late count floors by summing conservative passage times."""
import argparse
import gzip
import json
from pathlib import Path
import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq
from archived_collisionless_streaming_reference import Background
from conditional_relic_normalization_reference import conversion
from sm_dark_entropy_exchange_reference import MPL
from three_species_decay_onsager_reference import m, mp, M, y
from warm_late_visible_recycling_envelope_reference import (
    energy_rate_coefficient, COUP2, VECTOR_MASS, VECTOR_GAMMA,
)
from warm_prompt_reheating_expanding_slow_reference import visible
from warm_prompt_reheating_cold_band_loss_reference import vector_energy_ceiling


def coefficients(T, n, rho, split):
    h = np.pi*T*T/(np.sqrt(45)*MPL)
    ka = 4*y**4/(M*M)
    k0 = y**4/(8*np.pi*M*M)
    kb = ka/(3*np.sqrt(3)*np.pi*m*m)
    cp = (3*np.pi**2)**(2/3)/(2*np.pi**2)
    ke, _ = vector_energy_ceiling()
    p = VECTOR_MASS**2
    w2 = (VECTOR_MASS*VECTOR_GAMMA)**2
    kvac = 3*COUP2/(4*np.pi)*(np.sqrt(p*p+w2)+p)/(2*w2)
    ell = np.sqrt(1-4*split*split/(m*m))
    emin = m*m/(4*split)
    r2 = 1-(m/emin)**2
    return dict(r=rho/(m*n), a=k0*n/(2*h), b=2*ke*rho/h, c=kvac*n/(2*h),
                s=np.pi*ka*(m*n/mp)/(2*r2*ell*emin**3)*n/(2*h),
                j=np.pi*ka*cp*cp/(2*4**(2/3))*(m*n/split**2)*n**(1/3)/h,
                d=kb*m*m*n/h)


def loss(x, co, *, t=0.):
    r,a,b,c,s,j,d=(co[k] for k in ('r','a','b','c','s','j','d'))
    R=np.maximum(r-x,0.)
    return (a*x*x+np.minimum(b*x,c*x*x)+s*R*x*x*(1-t)**3+
            j*R*x**(4/3)*(1-t)+d*R*R)


def slab_passage_lower(x, co, slabs):
    # Independent endpoint rectangles for every positive factor.
    # min of the two vector upper envelopes is still an upper envelope.
    lo,hi=np.linspace(x,1.,slabs+1)[:-1],np.linspace(x,1.,slabs+1)[1:]
    r,a,b,c,s,j,d=(co[k] for k in ('r','a','b','c','s','j','d'))
    R=np.maximum(r-lo,0.)
    L=a*hi*hi+np.minimum(b*hi,c*hi*hi)+s*R*hi*hi+j*R*hi**(4/3)+d*R*R
    # Relative padding is much larger than a short sequence of IEEE operations.
    return float(np.sum((hi-lo)/L)*(1-1e-11))


def certified_floor(co, slabs):
    if slab_passage_lower(0.,co,slabs)<=1.:
        return 0.
    root=brentq(lambda x:slab_passage_lower(x,co,slabs)-1.,0.,1.,xtol=5e-13)
    safe=max(root-1e-10,0.)
    assert slab_passage_lower(safe,co,slabs)>=1.
    return safe


def point(bg,T, st, eta):
    _,sv=visible(bg,T)
    H=np.pi*T*T/(np.sqrt(45)*MPL)
    Tcap=np.exp((bg.gs(np.log(T))-np.log(2))/3)*T
    ce,_=energy_rate_coefficient(Tcap)
    n=sv*st['yields'][0]
    roh=sv*sum(st['energy_per_visible_entropy'])
    rho=(roh+ce*Tcap**9/(4*H))*np.exp(eta*Tcap**5/(3*H))
    return n, rho


def run(eos,traj,prior):
    bg=Background(eos)
    normalize=conversion(eos)['Omega_h2_per_GeV_times_yield']*m
    with gzip.open(traj,'rt') as f:
        paths=json.load(f)['arms']
    previous=json.loads(prior.read_text())
    eta=288*COUP2/(np.pi**3*VECTOR_MASS**4)
    ans=[]
    for st in paths[-1]['trajectory']:
        T=st['T_visible']
        n,rho=point(bg,T,st,eta)
        norm=normalize*st['yields'][0]
        vals=[]
        for split in (.002,.004,.0045,.0049):
            co=coefficients(T,n,rho,split)
            coarse=certified_floor(co,128)
            fine=certified_floor(co,2048)
            assert fine>=coarse-2e-10
            vals.append((fine,split,co,coarse))
        fine,split,co,coarse=max(vals, key=lambda x:x[0])
        exact=quad(lambda x:1./loss(x,co),fine,1.,epsabs=2e-10)[0]
        assert exact+1e-11>=slab_passage_lower(fine,co,2048)
        # Independent continuous comparator with the two positive dilution
        # factors kept, not used to assert the certified passage bound.
        ode=solve_ivp(lambda u,x:[-loss(max(float(x[0]),0.),co,t=u)],
                      (0.,1.),[1.],rtol=1e-10,atol=1e-12)
        assert ode.success and max(0.,ode.y[0,-1])>=fine-1e-8
        ans.append(dict(T_visible_GeV=T, n0=n, rho_cap=rho, split_GeV=split,
                        coeffs=co, coarse_retained_lower=coarse,
                        certified_retained_lower=fine,
                        omega_h2_lower=norm*fine,
                        continuous_comparator_with_dilution_diagnostic=max(0.,float(ode.y[0,-1]))))
    late=ans[-1]
    assert late['omega_h2_lower']>previous['anchors'][-1]['omega_h2_lower']
    tol=[]
    for a in [x for x in ans if x['omega_h2_lower']>.12]:
        T=a['T_visible_GeV']; n=a['n0']; rho=a['rho_cap']
        base=next(x['yields'][0] for x in paths[-1]['trajectory'] if x['T_visible']==T)*normalize
        for target in (.12, .5*base, .8*base):
            v=target/base
            def passage(mult):
                return slab_passage_lower(v,coefficients(T,n,rho*mult,a['split_GeV']),2048)
            cap=(brentq(lambda mult:passage(mult)-1.,1.,10.)-1e-9
                 if passage(1.)>1. else None)
            tol.append(dict(T_visible_GeV=T, target_floor=target,
                            sufficient_energy_ceiling_multiplier=cap))
    return dict(scope="Same conditional cap-free isotropic Boltzmann model and fixed charged-lepton "
                      "inputs as the earlier no-cap bound. Only the sampled anchors are tested. "
                      "Energy robustness scans change the cap at fixed n0, not the production history.",
                sampled_anchors=ans,rho_sensitivity=tol)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('eos','trajectory','previous'):
        p.add_argument(k,type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.trajectory,a.previous),indent=2,allow_nan=False))
