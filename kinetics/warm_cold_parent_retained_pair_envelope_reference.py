"""Support and shell-number upper on tagged parents from dark pair gains."""
import argparse
import json

import numpy as np
from pathlib import Path
from scalar_pair_expanding_continuation_reference import moments
from three_species_decay_onsager_reference import M,m,mp,y,width,P
from warm_light_occupation_stage_audit_reference import read, digest
from warm_near_rest_light_inverse_decay_reference import window

def scalar_amp_over_E1E2(smax,n=384):
 se=np.linspace(4*M*M,smax,n+1)
 U=np.linspace(0,1,2*n+1)
 sl,sh=se[:-1,None],se[1:,None]
 ul,uh=U[None,:-1],U[None,1:]
 uc=np.clip(.5,ul,uh)
 lo,hi=sl/2-M*M-mp*mp+m*m,sh/2-M*M-mp*mp+m*m
 pl,ph=sl/4-M*M,sh/4-M*M
 kl,kh=sl/4-mp**2,sh/4-mp**2
 ll=(M+m)*lo-2*M*kh*uh;lh=(M+m)*hi-2*M*kl*ul
 G=(M*M-mp*mp)**2+m*m*(sl-2*M*M-2*mp*mp)+m**4
 D=G+4*pl*kl*(1-uh)
 return (8*y**4*ph*(sh*kh**2*uc*(1-uc)+np.maximum(ll**2,lh**2))/(D**2*(sl/4-mp**2/2))).max()

def light_amp_over_E1E2(smax,eH):
 tt=(M*M-m*m)**2/smax
 return y**4*(1+(M+m)**2/(tt+mp*mp))**2/(m*eH)

def bound(S,n):
    p=np.asarray(S["lab_momentum_edges_GeV"])
    _,_,_,lo,hi=window(p[1])
    Et=m+lo; Emax=np.hypot(p[1],m)+hi
    ph=np.sqrt(Et*Et-M*M)
    Ehard=(Et+M)/2
    mom=moments(p,[S["entropy_variable_scalar"],S["entropy_variable_fermion"]])
    result=[]
    for j,d in enumerate((mp,m)):
        E=np.hypot(p,d)
        if 2*E[-1]<=Et+M:
            result.append(0.);continue
        N=mom[j,0]*(2*np.pi**2)/(1 if j==0 else 4)
        a=N.sum()
        Nt=N[E[1:]>=Ehard].sum()
        pmin=np.sqrt(Ehard*Ehard-d*d)
        A=(scalar_amp_over_E1E2(4*E[-1]**2,n) if j==0 else light_amp_over_E1E2(4*E[-1]**2,Ehard))
        # Four times the heavy spin-averaged matrix element overcounts
        # the entire spin sum at fixed tag; omit identical-initial 1/2.
        C=4*A*a*Nt/pmin/(32*np.pi**3*Et*ph)
        # Fermi rearrangement on all represented light radial number:
        # integral p f dp <= (3*N_L)^(2/3)/2, even for narrow saturated shells.
        N_L=mom[1,0].sum()*np.pi**2/2
        angular_L=(3*N_L)**(2/3)/(2*m)/(2*ph*P/M)
        assert angular_L<.001
        stat=C*Emax/(width*M*(1-angular_L))*(1+2e-10)
        result.append(stat)
    return result

def run(D):
    baseline=read(D/"warm-near-rest-light-inverse-decay-2026-11-26.json.gz")
    vis=read(D/"warm-cold-light-visible-source-margin-2026-11-26.json")
    exact=read(D/"warm-cold-light-parent-dark-pair-source-2026-11-26.json.gz")
    cases=[]
    for q in (16,32):
        tapes=read(D/f"warm-light-occupation-next-q{q}-tape-2026-11-26.json.gz")
        tapes+=[x["checkpoint"] for x in read(D/f"warm-observed-suboperator-q{q}-tape-2026-11-26.json.gz")["blocks"]]
        base=next(x["stages"] for x in baseline["cases"] if x["q"]==q)
        rows=[]
        for j,(c,a) in enumerate(zip(tapes,base,strict=True)):
            S=c["snapshot"]
            assert a["snapshot_sha256"]==c["snapshot_sha256"]==digest(S)
            x=np.array([bound(S,r) for r in (192,384,768)])
            assert np.all(x[-1]<=x[0]*(1+1e-13))
            f0=a["fine"]
            # Use the nonnegative gain derivative at the frozen reference.
            aa,bb=np.asarray(S["entropy_variable_scalar"])
            smax=float(1/np.expm1(min(np.min(aa),np.min(aa+bb))))
            possible=(M/m)**3*(1+smax)*x[-1].sum()/f0["H_over_rest_decay_width"]
            out=dict(T_GeV=S["T_visible_GeV"],interval_orders=[192,384,768],
                     parent_stationary_dark_pair_upper=x.tolist(),
                     upper_light_positive_rate_per_Hubble_at_frozen_baseline=possible)
            if q==32:
                v=vis["rows"][j]["sampled_upper_net_including_injection_per_Hubble"]
                assert abs(S["T_visible_GeV"]/vis["rows"][j]["T_GeV"]-1)<1e-12
                assert v+possible<0.
                out["augmented_frozen_sampled_net_upper_per_Hubble"]=v+possible
                out["upper_added_dark_pair_share_of_prior_frozen_residual"]=possible/-v
            if j in (0,3):
                archived=next(z for z in exact["results"] if (z["q"],z["j"])==(q,j))
                assert np.all(np.asarray(archived["orders"][-1]["divided_by_lower_decay_damping"])<=x[-1])
            rows.append(out)
        cases.append(dict(q=q,stages=rows))
    return dict(scope="Uniform upper on added retained dark pair creation using only "
                      "number, support, vacuum matrix-element envelopes and Pauli <= 1. "
                      "Supplemental signed-rate sum holds the decay-only reference frozen, "
                      "does not prove a signed-rate bound for re-shaped daughter spectra.",
                cases=cases)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--tapes",type=Path,default=Path("results"))
    print(json.dumps(run(p.parse_args().tapes),indent=2,allow_nan=False))
