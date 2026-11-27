"""Preserve the exact saved phase while splitting its outer prompt mesh."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre

from joint_prompt_two_moment_evolution_reference import Collision,kin
from joint_resonant_frw_split_reference import pack
from joint_visible_heavy_adiabatic_reference import snapshot_functions
from late_radiation_background_reference import MPL,TabulatedRadiation
from three_species_decay_onsager_reference import M,m,mp
from warm_corrected_late_prompt_cut_aligned_reference import observable_weights
from warm_coupled_serial_restart_reference import atomic
from warm_direct_prompt_stress_reference import DirectPrompt
from warm_light_occupation_stage_audit_reference import digest,read
from warm_prompt_excess_energy_tail_reference import integrate
from warm_prompt_reentry_shape_audit_reference import get_mom


def split(S,minimum,factor):
    if factor==1:
        return S,0.
    old=np.asarray(S["lab_momentum_edges_GeV"])
    t=np.linspace(0.,1.,factor+1)[1:-1]
    selected=old[:-1]>=minimum
    p=np.unique(np.r_[old,(old[:-1,None]+np.diff(old)[:,None]*t)[selected].ravel()])
    i=np.searchsorted(old,(p[:-1]+p[1:])/2,side="right")-1
    np.testing.assert_equal(old[i]<=p[:-1],True)
    np.testing.assert_equal(old[i+1]>=p[1:],True)
    th=[]
    for mass,k in ((mp,"scalar"),(m,"fermion")):
        u=kin(old,mass)
        x=kin(p,mass)
        a,b=np.asarray(S["entropy_variable_"+k])
        width=np.diff(u)[i]
        th.append(np.array([a[i]+b[i]*(x[:-1]-u[i])/width,b[i]*np.diff(x)/width]))
    old_m=get_mom(old,[S["entropy_variable_"+k] for k in ("scalar","fermion")])
    new_m=get_mom(p,th)
    residual=max(float(np.max(abs(np.sum(a,axis=1)/np.sum(b,axis=1)-1.))) for a,b in zip(new_m,old_m))
    assert residual<2e-12
    T=dict(S,lab_momentum_edges_GeV=p.tolist(),entropy_variable_scalar=th[0].tolist(),entropy_variable_fermion=th[1].tolist())
    a,x,b,y=snapshot_functions(S)
    A,X,B,Y=snapshot_functions(T)
    z,_=roots_legendre(6)
    pp=p[:-1,None]+np.diff(p)[:,None]*(1+z)/2
    for d,u,v in ((mp,a,A),(m,b,B)):
        E=np.hypot(d,pp)
        np.testing.assert_allclose(u(E),v(E),rtol=3e-12,atol=1e-245)
    return T,residual


def run(eos,gate,seed,output):
    A=read(seed)
    S=A["snapshot"]
    assert A["hash"]==digest(S)
    C=np.hypot(M,1./read(gate)["branches"][-1]["entropy_expansion_from_original_anchor"])
    rows=[]
    for n in (1,2):
        T,res=split(S,.4299838884938,n)
        p=np.array(T["lab_momentum_edges_GeV"])
        state=[np.asarray(T["entropy_variable_"+k]) for k in ("scalar","fermion")]
        r=get_mom(p,state)[2].sum()
        H=np.sqrt((TabulatedRadiation(eos,S["T_visible_GeV"]).state(1.)[1]+r)/3)/MPL
        note=[]
        for o,a in ((5,4),(6,6)):
            op=DirectPrompt(Collision(p,o,a),H,[],True)
            w=observable_weights(op,C)
            Q=pack(state)
            note.append(dict(radial_and_angular=[o,a],first_native_tail_flux=float(w@op.collision(Q)[0]/r),
                             tail_energy_over_rho=float(w@op.moment(Q)[0]/r)))
        test=integrate(T,C,12,True)/(r*H)
        if n==2 and output is not None:
            R=dict(A,snapshot=T,hash=digest(T))
            atomic(output,R)
        rows.append(dict(factor=n,shells=len(p)-1,accepted_exact_phase_moment_defect=res,
                         independent_angular_flux=float(test[3]),native=note,snapshot_sha256=digest(T)))
    return dict(scope="Same exact initial entropy functions on every partition; rate comparison is frozen, not a time convergence proof.",
                source_sha256=A["hash"],controls=rows)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("eos","gate","archive"):
        p.add_argument(k,type=Path)
    p.add_argument("--split-two",type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.gate,a.archive,a.split_two),indent=2,allow_nan=False))
