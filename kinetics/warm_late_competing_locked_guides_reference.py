"""Save all native prompt guides from one archived cut-aligned interval."""

import argparse
from pathlib import Path

import numpy as np

from joint_prompt_two_moment_evolution_reference import Collision
from joint_resonant_centered_stiff_reference import advance
from joint_resonant_frw_split_reference import pack
from late_radiation_background_reference import MPL, TabulatedRadiation
from warm_compact_absolute_prompt_reference import LocalAbsolutePrompt
from warm_corrected_late_prompt_cut_aligned_reference import observable_weights, predictor
from warm_coupled_serial_restart_reference import atomic, sha
from warm_light_occupation_stage_audit_reference import digest, read
from warm_prompt_reentry_shape_audit_reference import get_mom


def run(eos,start,target,out,n):
    a,b=read(start),read(target)
    assert a["plan"]==b["plan"] and a["plan"]["eos_sha256"]==sha(eos)
    assert all(x["hash"]==digest(x["snapshot"]) for x in (a,b))
    s=a["snapshot"]
    p=np.asarray(s["lab_momentum_edges_GeV"])
    assert np.array_equal(p,b["snapshot"]["lab_momentum_edges_GeV"])
    T=s["T_visible_GeV"]
    keys=("entropy_variable_scalar","entropy_variable_fermion")
    th=[s[k] for k in keys]
    rho=get_mom(p,th)[2].sum()
    H=np.sqrt((TabulatedRadiation(eos,T).state(1.)[1]+rho)/3.)/MPL
    op=LocalAbsolutePrompt(Collision(p),H,[],True)
    W=observable_weights(op,a["plan"]["normalized_cut_GeV"])
    u=pack(th)
    v=pack([b["snapshot"][k] for k in keys])
    out.mkdir(exist_ok=True,parents=True)
    atomic(out/"step-0.json",a)
    h=(b["time"]-a["time"])/n
    for j in range(1,n+1):
        u,meta=advance(op,u,h,T,initial_theta=predictor(op,u,h),vacuum_scaled=True,balanced_merit=True)
        t=b["time"] if j==n else a["time"]+j*h
        S=dict(s,**{k:r.tolist() for k,r in zip(keys,op.arrays(u),strict=True)})
        atom=dict(snapshot=S,hash=digest(S),plan=a["plan"],time=t,
                  signed_native_flux=float(W@op.collision(u)[0]/rho),solver=meta,
                  source_start_hash=a["hash"],source_end_hash=b["hash"])
        if j==n:
            N=op.moment(u)[0].reshape(2,op.n,2)
            N1=op.moment(v)[0].reshape(2,op.n,2)
            err=float(np.max((abs(N-N1).sum(axis=1)/N.sum(axis=1))))
            assert err<2e-12
            atom["archived_terminal_relative_moment_L1"]=err
        atomic(out/f"step-{j}.json",atom)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("eos","start","target","output"):
        p.add_argument(k,type=Path)
    p.add_argument("--steps",type=int,required=True)
    a=p.parse_args()
    assert a.steps>0
    run(a.eos,a.start,a.target,a.output,a.steps)
