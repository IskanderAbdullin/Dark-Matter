"""Probe the initial Newton systems of one matched coupled block, without advancing it."""

import argparse
import json
from pathlib import Path

import numpy as np

from frw_exterior_restart_reference import import_restart
from joint_below_pole_elastic_linear_modes_reference import setup
from joint_prompt_two_moment_evolution_reference import Collision
from joint_resonant_frw_split_reference import pack,redshift
from late_archived_shape_continuation_reference import Late
from late_radiation_background_reference import MPL,TabulatedRadiation
from scalar_pair_expanding_continuation_reference import e_total,moments,rate,reconstruct
from warm_light_occupation_stage_audit_reference import digest,read


def prepare(eos,S,step):
    p=np.asarray(S["lab_momentum_edges_GeV"])
    th=[np.asarray(S["entropy_variable_"+a]) for a in ("scalar","fermion")]
    a=np.sqrt(step)
    th,*_=redshift(p,1.,a,th)
    q=p/a
    U=moments(q,th)
    bg=TabulatedRadiation(eos,S["T_visible_GeV"])
    T,r=bg.state(a)
    ex=import_restart(S["positive_exterior_restart"])
    H=np.sqrt((r+e_total(U)+ex.energy(a)/a**4)/3)/MPL
    h=np.log(step)
    A,_,_=rate(q,th,(3,7,6,2))
    mid=reconstruct(q,U+h*A[:,:,1:-1]/(4*H))
    B,_,_=rate(q,mid,(3,7,6,2))
    th=reconstruct(q,U+h*B[:,:,1:-1]/(2*H))
    v=pack(th)
    seed=Late(Collision(q),H,[],True)
    _,_,thetaT=seed.reference(v,T)
    panels=[]
    setup(max(thetaT,T),seed.k,thetaT,seed.eta,3,2,2,8,32.,.82,cache=panels,cache_only=True)
    model=Late(seed.k,H,panels,True,compact=True,local_jac=True)
    a,b,_=model.reference(v,T)
    return model,v,a,b,h,T,H


def probe(mod,v,ref,R0,h,eq):
    U0=mod.moment(v)[0]
    init=ref if eq else v
    z=init-ref
    U,_,C=mod.moment(init)
    F,jvp=mod.centered(z,ref,R0)
    weights=np.repeat(np.maximum(U0[::2],U[::2]),2)
    assert np.min(weights)>0.
    r=(U-U0-h*F)/weights
    def din(w):
        a=w.reshape(-1,2)
        return (-np.stack((C[:,0]*a[:,0]+C[:,1]*a[:,1],C[:,1]*a[:,0]+C[:,2]*a[:,1]),axis=1).ravel()-h*jvp(w))/weights
    matrix=np.stack([din(d) for d in np.eye(mod.d)],axis=1)
    row=1/np.max(abs(matrix),axis=1)
    M=row[:,None]*matrix
    col=1/np.max(abs(M),axis=0)
    M*=col
    s0=np.linalg.svd(matrix,compute_uv=False)
    s1=np.linalg.svd(M,compute_uv=False)
    da0=np.linalg.solve(matrix,-r)
    da1=col*np.linalg.solve(M,-row*r)
    relative=lambda X: float(np.linalg.norm(matrix@X+r)/max(np.linalg.norm(r),1e-280))
    i=int(np.argmax(abs(r)))
    species=int(i>=2*mod.n)
    cell=(i%(2*mod.n))//2
    return dict(equilibrium_trial=eq,max_residual=float(np.max(abs(r))),species_and_shell=[species,cell],
                shell_momentum_range_GeV=mod.k.c[species].p[cell:cell+2].tolist(),
                min_to_max_number_ratio=float(np.min(weights)/np.max(weights)),
                singular_condition_original=s0[0]/s0[-1],singular_condition_equilibrated=s1[0]/s1[-1],
                ordinary_linear_residual=relative(da0),equilibrated_linear_residual=relative(da1),
                difference_between_updates_over_norm=float(np.linalg.norm(da0-da1)/max(np.linalg.norm(da1),1e-280)))


def run(eos,D,q64):
    C=read(D/"warm-light-occupation-next-q32-tape-2026-11-26.json.gz")
    A=C[1]
    B=read(q64)
    X=A["snapshot"];Y=B["snapshot"]
    assert A["snapshot_sha256"]==digest(X) and B["snapshot_sha256"]==digest(Y)
    assert abs(X["T_visible_GeV"]/Y["T_visible_GeV"]-1)<1e-14
    step=C[2]["snapshot"]["scale_factor"]/X["scale_factor"]
    rows=[]
    for tag,S in (("q32",X),("q64",Y)):
        a,v,ref,R0,h,T,H=prepare(eos,S,step)
        rows.append(dict(tag=tag,log_step=h,T_mid_GeV=T,H_GeV=H,compact_panels=len(a.compiled),
                         trials=[probe(a,v,ref,R0,h,eq) for eq in (False,True)]))
    return dict(scope="Initial frozen nonlinear residual and row-scaled entropy-coordinate Jacobian "
                      "of one implicit suboperator, after the actual redshift and first pair kick. "
                      "No iteration, accepted block or cause of a long runtime inferred.",cases=rows)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("eos",type=Path)
    p.add_argument("q64",type=Path)
    p.add_argument("--tapes",type=Path,default=Path("results"))
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.tapes,a.q64),indent=2,allow_nan=False))
