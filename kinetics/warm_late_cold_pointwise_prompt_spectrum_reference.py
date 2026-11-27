"""Resolve signs within the cold scalar band of frozen hard reactivations."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre

from decay_scalar_ir_relaxation_reference import limits
from joint_visible_heavy_adiabatic_reference import angle_panels,snapshot_functions
from late_radiation_background_reference import MPL,TabulatedRadiation
from three_species_decay_onsager_reference import Es,M,P,m,mp,width
from warm_corrected_exterior_reactivated_tail_stress_reference import reactivate
from warm_late_cold_pointwise_bose_dose_reference import envelope
from warm_light_occupation_stage_audit_reference import digest,read
from warm_prompt_reentry_shape_audit_reference import get_mom


def boundary_parent(energy,mass,rest):
    p=np.sqrt((energy-mass)*(energy+mass))
    hi=energy*rest+p*P
    return np.r_[M*(rest*rest+p*p)/hi,M*hi/(mass*mass)]


def spectrum(D,q,r):
    f,A,g,B=snapshot_functions(D)
    e,lo,span=limits(q)
    hi=min(lo+span,A[-1]+B[-1])
    if hi<=lo:
        return np.zeros(3)
    bx=np.unique(np.clip(np.r_[lo,hi,e+B,boundary_parent(A,mp,Es),boundary_parent(B,m,M-Es),
                                (A[:,None]+B).ravel()],lo,hi))
    x,w=roots_legendre(r)
    h=np.diff(bx)[:,None]/2
    E=(bx[:-1,None]+h*(1+x)).ravel()
    measure=(h*w).ravel()
    g1=f(e)
    acc=np.zeros(3)
    for j in range(0,len(E),96):
        v=E[j:j+96]
        a,W=angle_panels(v,A,B,angle=r)
        b=v[:,None,None]-a
        fp,fl=f(a),g(b)
        I=(W*fp*fl).sum(axis=(1,2))
        F=(W*(1.+fp)*(1.-fl)).sum(axis=(1,2))
        occ=np.divide(I,I+F,out=np.zeros_like(I),where=I+F>0)
        gamma=g(v-e)
        gain=occ*(1.+g1)*(1.-gamma)
        loss=(1.-occ)*g1*gamma
        mm=measure[j:j+96]
        acc+=np.array([mm@gain,mm@loss,mm@(gain-loss)])
    return 4*width*M*M/(2*P*q*e)*acc


def run(eos,file):
    v=read(file)
    S=v["snapshot"]
    assert v["snapshot_sha256"]==digest(S)
    T=S["T_visible_GeV"]
    rows=[]
    for h in (3.3,6.5):
        D,*_=reactivate(S,h,.02,cap=7.5)
        p=np.asarray(D["lab_momentum_edges_GeV"])
        _,_,E=get_mom(p,[D["entropy_variable_"+k] for k in ("scalar","fermion")])
        H=np.sqrt((TabulatedRadiation(eos,T).state(1.)[1]+E.sum())/3.)/MPL
        f,*_=snapshot_functions(D)
        qs=np.r_[.015,.025,np.linspace(.03,.1,15)]
        results=[]
        for q in qs:
            pair=[spectrum(D,q*T,r)/H for r in (6,10)]
            gross=np.max(pair[-1][:2])
            if gross:
                assert np.max(abs(pair[0]-pair[1]))<4e-4*gross
                assert abs(pair[-1][0]-pair[-1][1]-pair[-1][2])<4e-14*gross
            else:
                assert np.all(pair[0]==0.)
            occupation=float(f(np.hypot(q*T,mp)))
            results.append(dict(p_over_T=float(q),occupation=occupation,
                                 open_prompt_gain_loss_signed_over_H=[a.tolist() for a in pair],
                                 target_fraction_of_occupation_per_H=float(pair[-1][2]/occupation)))
        rows.append(dict(hard_center_GeV=h,seed_sha256=digest(D),
                         cold_reach_floor_over_T=envelope(.1*T,float(np.hypot(p[-1],mp)+np.hypot(p[-1],m)),24)["cold_reach_floor_GeV"]/T,
                         orders=[6,10],samples=results))
    return dict(scope="Frozen initial pointwise open prompt gain and loss with stationary parent. "
                      "Some rates are zero only for this represented support; no line between samples.",
                source_sha256=v["snapshot_sha256"],cases=rows)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("eos",type=Path)
    p.add_argument("archive",type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.archive),indent=2,allow_nan=False))
