"""Matched q32/q64 explicit prompt re-entry in a short expanding control."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from joint_prompt_two_moment_evolution_reference import Collision
from joint_resonant_frw_split_reference import pack,redshift
from late_archived_shape_continuation_reference import Late,collide
from late_radiation_background_reference import MPL,TabulatedRadiation
from warm_light_occupation_stage_audit_reference import digest,read
from warm_prompt_reentry_shape_audit_reference import get_mom,pad


def episode(eos,S,q,ncore,n,b):
    th=pad(S,q)
    bg=TabulatedRadiation(eos,S["T_visible_GeV"])
    N0,K0,E0=get_mom(q,th)
    INIT=E0.sum()
    Ns=N0.sum(axis=1)
    w=0.
    rows=[]
    for lo,hi in zip(np.linspace(0.,np.log(b),n+1)[:-1],
                     np.linspace(0.,np.log(b),n+1)[1:],strict=True):
        a,c,d=np.exp([lo,(lo+hi)/2,hi])
        th,s,*_=redshift(q,a,c,th)
        w+=s.sum()
        T,R=bg.state(c)
        nc,kc,ec=get_mom(q/c,th)
        H=np.sqrt((R+ec.sum())/3)/MPL
        model=Late(Collision(q/c),H,[],True,local_jac=True)
        v,res,ok,no=collide(model,pack(th),hi-lo,T,continuation=True,balanced=True)
        th,s,*_=redshift(q,c,d,model.arrays(v))
        w+=s.sum()
        N,K,E=get_mom(q/d,th)
        err=float(np.max(abs(N.sum(axis=1)*d**3/Ns-1)))
        led=float(abs(d**4*E.sum()-INIT-w)/INIT)
        assert err<5e-12 and led<5e-11
        rows.append(dict(scale_factor=d,T_GeV=bg.state(d)[0],q32_outside_number_fraction=(N[:,ncore:].sum(axis=1)/N.sum(axis=1)).tolist(),
                         q32_outside_energy_fraction=(E[:,ncore:].sum(axis=1)/E.sum(axis=1)).tolist(),
                         comoving_number_error=err,energy_work_closure=led,solver=[res,ok,no]))
    return dict(blocks=n,steps=rows,relative_number_per_q32_shared_shell=(N[:,:ncore]*b**3/Ns[:,None]).tolist(),
                energy_per_q32_shared_shell_normalized_to_initial_total=(E[:,:ncore]*b**4/INIT).tolist())


def run(eos,D,b):
    C=read(D/"warm-light-occupation-next-q32-tape-2026-11-26.json.gz")[0]
    S=C["snapshot"]
    assert digest(S)==C["snapshot_sha256"]
    q=np.asarray(S["lab_momentum_edges_GeV"])
    wide=np.r_[q,q[-1]*np.geomspace(1.,2.,6)[1:]]
    args=[(q,4),(q,8),(wide,4),(wide,8)]
    res=[]
    for P,n in args:
        R=episode(eos,S,P,len(q)-1,n,b)
        R["support"]="q32" if len(P)==len(q) else "q64"
        res.append(R)
    def err(a,b,key):
        return np.sum(abs(np.asarray(a[key])-np.asarray(b[key])),axis=1).tolist()
    keys=["relative_number_per_q32_shared_shell","energy_per_q32_shared_shell_normalized_to_initial_total"]
    measures={key:dict(q32_step_change=err(res[0],res[1],key),q64_step_change=err(res[2],res[3],key),
                       added_support_at_fine_time_step=err(res[1],res[3],key)) for key in keys}
    return dict(scope="Matched short expanding trial with only prompt transfer, using real q32 "
                     "retained archive plus explicit vacuum padding. Excluded inherited exterior "
                     "and all other kinetic channels, no long-run conclusion. Radiation clock uses "
                     "conserved visible entropy; dark energy only changes the Hubble normalization.",
                snapshot_sha256=C["snapshot_sha256"],EOS_sha256=hashlib.sha256(eos.read_bytes()).hexdigest(),
                final_relative_b=b,res=res,relative_L1_changes=measures)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("eos",type=Path)
    p.add_argument("--tapes",type=Path,default=Path("results"))
    p.add_argument("--b",type=float,default=1.0001)
    ns=p.parse_args()
    print(json.dumps(run(ns.eos,ns.tapes,ns.b),indent=2,allow_nan=False))
