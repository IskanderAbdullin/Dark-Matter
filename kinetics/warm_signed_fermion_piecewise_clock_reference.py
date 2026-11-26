"""Refine the conditional signed warm comparison with local radiation floors."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from conditional_relic_normalization_reference import conversion
from late_radiation_background_reference import TabulatedRadiation
from three_species_decay_onsager_reference import m
from warm_isolated_recycling_energy_envelope_reference import radiation_bound
from warm_light_occupation_stage_audit_reference import digest,read
from warm_prompt_reentry_shape_audit_reference import get_mom
from warm_signed_fermion_time_comparison_reference import comparison,setup


def weights(eos,T,b_stop,steps,Hph,powers=(3,4,5)):
    R=TabulatedRadiation(eos,T)
    u=R.rg.x
    lo=np.log(R.state(b_stop)[0])
    u=np.r_[lo,u[(u>lo)&(u<R.initial)],R.initial]
    rho,s=R.rg(u),R.sg(u)
    bb=((R.s0-s-3*u)/3)[::-1]
    eta=np.exp(.5*(np.minimum(rho[:-1],rho[1:])-np.log(2.))+
               2/3*(s[-1]-np.maximum(s[:-1],s[1:]))-1e-12)[::-1]
    begin=np.linspace(0.,np.log(b_stop),steps+1)
    assert abs(bb[0])<3e-15 and abs(bb[-1]-begin[-1])<2e-12
    bb[[0,-1]]=begin[[0,-1]]
    q=np.unique(np.r_[begin,bb])
    md=(q[1:]+q[:-1])/2
    j=np.searchsorted(bb,md,"right")-1
    i=np.searchsorted(begin,md,"right")-1
    a=np.empty((steps+1,len(powers)))
    for c,p in enumerate(powers):
        x=p-2
        g=np.exp(-x*q[:-1])*(-np.expm1(-x*np.diff(q)))/(x*Hph*eta[j])
        a[:-1,c]=np.bincount(i,weights=g,minlength=steps)*(1+2e-14)
        a[-1,c]=b_stop**(-x)/(x*Hph)
    return a


def assess(eos,T,n,E,K,Kpi,z,label):
    clk=radiation_bound(eos,T,500.)
    out=[]
    for ns in (350,1000,2800):
        A,B,*r=setup(n,E,K,clk,ns)
        B=B+.5*n*Kpi
        new=weights(eos,T,500.,ns,clk["photon_anchor_GeV"])
        assert np.all(new<=A*(1+2e-10))
        errors=[]
        for cap in (.01,.1,.5):
            uniform=comparison((A,B,*r),cap)[0]*z
            piece=comparison((new,B,*r),cap)[0]*z
            assert piece>=uniform and piece>.12
            errors.append(dict(cap=cap,uniform_lower=uniform,local_floor_lower=piece))
        out.append(dict(time_steps=ns,lower_Omega_h2=errors,
                        integrated_weight_over_uniform=(new.sum(axis=0)/A.sum(axis=0)).tolist()))
    # As a sensitivity control only, change the gross available-energy
    # ceiling while keeping the measured starting N and all loss kernels.
    def diagnostic(f):
        A,B,*r=setup(n,E*f,K,clk,2800)
        return comparison((new,B+.5*n*Kpi,*r),.5)[0]*z
    critical=brentq(lambda f:diagnostic(f)-.120,1.,3.,xtol=1e-10)
    return dict(label=label,T_GeV=T,energy_upper_GeV4=E,
                number_at_anchor_GeV3=n,unconverted_start_proxy=z,controls=out,
                energy_multiplier_first_losing_exclusion_at_cap_half=critical)


def run(eos,D,later=None):
    first=read(D/"warm-isolated-fermion-survival-bound-2026-11-26.json")["cases"][-1]
    L=read(D/"warm-fermion-survival-with-leptons-2026-11-26.json")
    W=read(D/"warm-signed-fermion-with-charged-pions-2026-11-26.json")["cases"][-1]["controls"][-1]
    oldT=first["anchor_T_GeV"]
    base=assess(eos,oldT,first["core_start_fermion_number_GeV3"],W["new_available_energy_upper_GeV4"],
                L["coefficient"],W["pion_pair_sigma_v_cap_with_stated_support"],
                first["conservative_initial_fermion_Omega"],"old q32 anchor")
    assert abs(base["controls"][-1]["lower_Omega_h2"][-1]["uniform_lower"]/
               W["with_charged_pions_lower_Omega_h2"]-1)<1e-12
    A=read(later or D/"warm-q64-history-serial-b-2026-11-26.json.gz")[-1]
    s=A["snapshot"]
    assert A["snapshot_sha256"]==digest(s)
    T=s["T_visible_GeV"]
    rad=TabulatedRadiation(eos,T)
    early=TabulatedRadiation(eos,oldT)
    b=np.exp((early.s0-rad.s0)/3)
    N=get_mom(np.asarray(s["lab_momentum_edges_GeV"]),
              [s["entropy_variable_scalar"],s["entropy_variable_fermion"]])[0][1].sum()
    C=conversion(eos)["Omega_h2_per_GeV_times_yield"]
    V=m*N*C/(2*np.pi**2/45*np.exp(rad.s0))
    if "anchor_T_GeV" not in first: raise ValueError("unknown first energy anchor")
    E=W["new_available_energy_upper_GeV4"]/b**3
    new=assess(eos,T,N,E,L["coefficient"],W["pion_pair_sigma_v_cap_with_stated_support"],V,
               "later q64 checkpoint, conditional")
    assert new["controls"][-1]["lower_Omega_h2"][-1]["local_floor_lower"]>base["controls"][-1]["lower_Omega_h2"][-1]["local_floor_lower"]
    return dict(scope="Both fermion species remain below the prescribed cap and use the specified "
                      "homogeneous isotropic selected visible budget. This is not proof of the cap "
                      "or of the true post-checkpoint history. The prior full available energy "
                      "bound is transported with b^3 rho; no new outside energy is assumed.",
                later_snapshot_sha256=A["snapshot_sha256"],
                later_over_old_expansion_from_radiation_entropy=b,cases=[base,new])


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("eos",type=Path)
    p.add_argument("--tapes",type=Path,default=Path("results"))
    p.add_argument("--later",type=Path);a=p.parse_args()
    print(json.dumps(run(a.eos,a.tapes,a.later),indent=2,allow_nan=False))
