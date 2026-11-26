"""Audit matched lower-temperature q32 and q64 runs and the optional time ladder."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre

from conditional_relic_normalization_reference import conversion
from frw_exterior_restart_reference import import_restart
from joint_visible_heavy_adiabatic_reference import snapshot_functions
from late_radiation_background_reference import MPL,TabulatedRadiation
from three_species_decay_onsager_reference import m,mp
from warm_coupled_serial_restart_reference import sha
from warm_light_occupation_stage_audit_reference import digest,read
from warm_observed_suboperator_continuation_reference import metrics
from warm_prompt_high_parent_ledger_reference import high
from warm_prompt_reentry_shape_audit_reference import get_mom


def measures(s):
    p=np.asarray(s["lab_momentum_edges_GeV"])
    B,_,F,_=snapshot_functions(s)
    N,K,E=get_mom(p,[s["entropy_variable_scalar"],s["entropy_variable_fermion"]])
    z,w=roots_legendre(24)
    t=(p[:-1,None]+p[1:,None])/2+np.diff(p)[:,None]*z/2
    W=np.diff(p)[:,None]*w/2
    P=[]
    for d,g,fn in ((mp,1,B),(m,4,F)):
        u=np.hypot(t,d)
        P.append(g/(6*np.pi**2)*np.sum(W*t**4*fn(u)/u,axis=1))
    return np.array([N,E,P])


def branch(eos,seed,work,steps,q):
    x=read(seed)
    S=x[-1]["snapshot"]
    assert x[-1]["snapshot_sha256"]==digest(S)
    I=measures(S)
    a0=S["scale_factor"]
    ref=I[0].sum()*a0**3
    old=import_restart(S["positive_exterior_restart"]).array()
    n0=old.shape[1]
    ranges=[]; V=[]
    last=digest(S)
    guide_count=0
    for i in range(1,steps+1):
        C=read(work/f"step-{i}.json")
        assert C["index"]==i and C["previous"]==last
        assert C["plan"]["tape_sha256"]==sha(seed) and C["plan"]["EOS_sha256"]==sha(eos)
        s=C["snapshot"]; last=digest(s)
        assert C["snapshot_sha256"]==last
        z=s["scale_factor"]; b=z/a0
        assert abs(b/(1.08**(i/steps))-1.)<7e-14
        x=measures(s)
        ext=import_restart(s["positive_exterior_restart"])
        ar=ext.array()
        check=ar[:,:n0].copy()
        check[[0,5]]*=np.array([b,b**3])[:,None]
        np.testing.assert_allclose(check,old,rtol=1.2e-13,atol=0.)
        number=abs((x[0].sum()+ar[5,n0:].sum())*z**3/ref-1.)
        assert number<5e-13
        temp=TabulatedRadiation(eos,s["T_visible_GeV"])
        Svis=(2*np.pi**2/45)*np.exp(temp.s0)
        omega=conversion(eos)["Omega_h2_per_GeV_times_yield"]*m*x[0,1].sum()/Svis
        native=C["operator_summary"]
        audit=C["audit"]["solver"]
        guide_count+=int(audit["used_history"])
        ranges.append(dict(step=i,T_GeV=s["T_visible_GeV"],snapshot_sha256=last,
                           cumulative_count_residual=float(number),local_energy_defect=native["energy_ledger"],
                           history_iterations=audit.get("newton_iterations"),retained_stable_rest_Omega_h2_proxy=omega,
                           peak=metrics([s["entropy_variable_scalar"],s["entropy_variable_fermion"]])))
        V.append(x*z**3)
    return dict(rows=ranges,accepted_with_history=guide_count,original_kinetic_grade="selected",
                q=q,steps=steps),V,s


def run(D,eos,work,refined):
    A=read(D/"warm-light-occupation-next-q32-tape-2026-11-26.json.gz")[0]["snapshot"]["scale_factor"]
    # Same normalization for all runs: initial q32 count and full
    # dark energy. Pressure and energy are not normalized separately.
    start=read(work/"seed32.json")[-1]["snapshot"]
    init=measures(start)*start["scale_factor"]**3
    den=np.vstack((init[0].sum(axis=1),np.repeat(init[1].sum(),2),np.repeat(init[1].sum(),2)))
    output={}; data={}
    for n in ((16,32) if refined else (16,)):
        R=[]
        for q in (32,64):
            folder=work/(f"q{q}" if n==16 else f"q{q}-fine")
            z,v,s=branch(eos,work/f"seed{q}.json",folder,n,q)
            R.append((z,v,s))
            data[n,q]=v
        ix=R[0][1][-1].shape[-1]
        support=np.asarray([abs(x-y[:,:,:ix]).sum(axis=2)/den for x,y in zip(R[0][1],R[1][1],strict=True)])
        s=R[1][2]; narrow=R[0][2]
        np.testing.assert_allclose(narrow["lab_momentum_edges_GeV"],s["lab_momentum_edges_GeV"][:ix+1],rtol=1e-13)
        assert abs(s["T_visible_GeV"]/narrow["T_visible_GeV"]-1)<1e-13
        p=measures(s)
        E=p[1].sum()
        B=s["scale_factor"]/A
        H=np.sqrt((TabulatedRadiation(eos,s["T_visible_GeV"]).state(1.)[1]+E+
                   import_restart(s["positive_exterior_restart"]).energy(1.))/3)/MPL
        vals=np.array([high(s,1./B,g)[1:3] for g in (8,14)])/(E*H)
        assert np.max(abs(vals[1]/vals[0]-1.))<.03
        output[str(n)]=dict(controls=[z for z,_,_ in R],support_L1_by_step=support.tolist(),
                            q64_outer_number_fraction=(p[0,:,ix:].sum(axis=1)/p[0].sum(axis=1)).tolist(),
                            q64_frozen_bare_high_inverse_and_tagged_energy=vals.tolist())
    if refined:
        output["fine_step_L1_by_support"]={str(q):(abs(data[16,q][-1]-data[32,q][-1]).sum(axis=2)/den).tolist()
                                           for q in (32,64)}
    return dict(scope="All accepted physical target equations use their own grid independently, "
                      "only selected channels. Same expansion anchor and time ladders, "
                      "old known exterior used as gravity, not restored as a phase distribution. "
                      "Conditional stable rest proxy is not a cosmological prediction.",
                metric_order=["number","full_energy","pressure"],results=output)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("eos",type=Path);p.add_argument("work",type=Path)
    p.add_argument("--tapes",type=Path,default=Path("results"));p.add_argument("--refined",action="store_true")
    a=p.parse_args()
    print(json.dumps(run(a.tapes,a.eos,a.work,a.refined),indent=2,allow_nan=False))
