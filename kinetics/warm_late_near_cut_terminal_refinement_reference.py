"""Refine only the terminal interval from an unchanged accepted saved phase."""

import argparse
import json
from pathlib import Path

import numpy as np

from joint_prompt_two_moment_evolution_reference import Collision
from joint_resonant_centered_stiff_reference import advance
from joint_resonant_frw_split_reference import pack
from warm_compact_absolute_prompt_reference import LocalAbsolutePrompt
from warm_corrected_late_prompt_cut_aligned_reference import observable_weights, predictor
from warm_coupled_serial_restart_reference import atomic
from warm_light_occupation_stage_audit_reference import digest, read
from warm_prompt_excess_energy_tail_reference import integrate
from warm_prompt_reentry_shape_audit_reference import get_mom


def run(source, output, total_steps):
    ref=read(source/"summary.json")
    a,b=[read(source/f"step-{j}.json") for j in (5,6)]
    assert all(z["hash"]==digest(z["snapshot"]) and z["plan"]==ref["plan"] for z in (a,b))
    assert total_steps>=3
    S=a["snapshot"]
    H=ref["H_GeV"]
    cut=ref["plan"]["normalized_cut_GeV"]
    keys=("entropy_variable_scalar", "entropy_variable_fermion")
    p=np.asarray(S["lab_momentum_edges_GeV"])
    r=get_mom(p,[S[k] for k in keys])
    N0=r[0].sum(axis=1)
    rho0=r[2].sum()
    op=LocalAbsolutePrompt(Collision(p),H,[],True,stream_local=True)
    W=observable_weights(op,cut)
    u=pack([S[k] for k in keys])
    h=(b["time"]-a["time"])/total_steps
    rows=[]
    for j in range(1,total_steps+1):
        start=W@op.moment(u)[0]/rho0
        u,info=advance(op,u,h,S["T_visible_GeV"],initial_theta=predictor(op,u,h),
                       vacuum_scaled=True,balanced_merit=True)
        f=W@op.collision(u)[0]/rho0
        w=W@op.moment(u)[0]/rho0
        assert abs(w-start-h*f)<4e-10*max(w,start,1e-70)
        n,_,e=get_mom(p,op.arrays(u))
        assert max(*abs(n.sum(axis=1)/N0-1),abs(e.sum()/rho0-1))<8e-12
        Q=dict(S,**{k:x.tolist() for k,x in zip(keys,op.arrays(u),strict=True)})
        angular=integrate(Q,cut,8,True)/(rho0*H)
        assert abs(angular[3]-f)<3e-5*max(angular[1:3])
        rows.append(dict(index=j,Hdt=a["time"]+j*h,W_over_rho_from_saved_phase=float(w),
            signed_current=float(f), angular_signed_current=float(angular[3]),
            W_balance_error=float(w-start-h*f),solver=info,hash=digest(Q)))
    if total_steps==3:
        assert rows[-1]["hash"]==b["hash"]
    atomic(output,dict(scope="Refined last interval only; same accepted coarse prefix and strict entropy guard.",
        anchor_hash=a["hash"],coarse_terminal_hash=b["hash"],coarse_terminal_W=b["tail_energy"],
        coarse_terminal_current=b["signed_native_flux"],interval_substeps=total_steps,rows=rows,
        terminal_snapshot=Q))


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source",type=Path)
    parser.add_argument("output",type=Path)
    parser.add_argument("--substeps",type=int,default=6)
    args=parser.parse_args()
    run(args.source,args.output,args.substeps)
