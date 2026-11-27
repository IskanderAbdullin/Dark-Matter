"""Resolve prompt channel turnover rates per momentum shell."""
import argparse
import json
from pathlib import Path

import numpy as np

from joint_prompt_two_moment_evolution_reference import Collision
from late_radiation_background_reference import MPL,TabulatedRadiation
from warm_coupled_serial_restart_reference import sha
from warm_light_occupation_stage_audit_reference import digest,read
from warm_prompt_reentry_shape_audit_reference import get_mom


def run(eos,path):
    r=read(path)
    S=r["snapshot"]
    assert digest(S)==r["hash"] and r["plan"]["eos_sha256"]==sha(eos)
    p=np.asarray(S["lab_momentum_edges_GeV"])
    th=[np.asarray(S["entropy_variable_"+name]) for name in ("scalar","fermion")]
    T=S["T_visible_GeV"]
    N,_,E=get_mom(p,th)
    H=np.sqrt((TabulatedRadiation(eos,T).state(1.)[1]+E.sum())/(3.*MPL**2))
    k=Collision(p)
    # Avoid compiling the local derivative; this is an ordered-rate fold.
    F=k.c[0].f(th[0],k.ph,k.I)
    G=k.c[1].f(th[1],k.ch,k.J)
    a=k.pre*F*G
    b=k.pre*(1.+F)*(1.-G)
    A=np.bincount(k.y,weights=a,minlength=k.n)
    B=np.bincount(k.y,weights=b,minlength=k.n)
    h=np.divide(A,A+B,out=np.zeros_like(A),where=A+B>0.)[k.y]
    plus=k.fac/H*h*b
    minus=k.fac/H*(1.-h)*a
    assert abs(plus.sum()/minus.sum()-1.)<3e-14
    cases=[]
    for j,I in enumerate((k.I,k.J)):
        vals=np.stack([np.bincount(I,weights=z,minlength=len(p)-1) for z in (plus,minus,plus-minus)],axis=1)
        overlap=(p[1:]>=.795)&(p[:-1]<=.805)
        ratio=np.divide(vals,N[j,:,None],out=np.zeros_like(vals),where=N[j,:,None]>0.)
        core=[]
        for i in np.flatnonzero(overlap):
            core.append(dict(shell=int(i),edges=p[i:i+2].tolist(),density=float(N[j,i]),
                             gain_loss_and_signed_per_particle_per_Hubble_time=ratio[i].tolist()))
        cases.append(core)
    return dict(input_hash=r["hash"],hard_GeV=r["plan"]["hard_GeV"],time=r["time"],
                columns=["scalar","light"],scope="Frozen positive rates, not stability eigenvalues or a timestep positivity criterion.",
                shells_near_point805=cases)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("eos",type=Path);p.add_argument("paths",nargs="+",type=Path)
    a=p.parse_args()
    print(json.dumps([run(a.eos,x) for x in a.paths],indent=2,allow_nan=False))
