"""Replay the complete fixed six-substep clock with explicit entropy certification."""

import argparse
import json
from pathlib import Path

import numpy as np

from joint_prompt_two_moment_evolution_reference import Collision
from joint_resonant_centered_stiff_reference import advance
from joint_resonant_frw_split_reference import pack
from warm_compact_absolute_prompt_reference import LocalAbsolutePrompt
from warm_corrected_late_prompt_cut_aligned_reference import excess, observable_weights, predictor
from warm_coupled_serial_restart_reference import atomic
from warm_late_near_cut_early_entropy_audit_reference import entropy_line, stable_total
from warm_light_occupation_stage_audit_reference import digest, read
from warm_prompt_excess_energy_tail_reference import integrate
from warm_prompt_reentry_shape_audit_reference import get_mom


def run(source, out):
    a=read(source/"step-0.json")
    coarse=read(source/"summary.json")
    S=a["snapshot"]
    assert a["hash"]==digest(S)==coarse["plan"]["seed_sha256"]
    p=np.asarray(S["lab_momentum_edges_GeV"])
    keys=("entropy_variable_scalar","entropy_variable_fermion")
    sh=[S[k] for k in keys]
    N0,_,E0=get_mom(p,sh)
    rho=E0.sum()
    H=coarse["H_GeV"]
    cut=coarse["plan"]["normalized_cut_GeV"]
    op=LocalAbsolutePrompt(Collision(p),H,[],True,stream_local=True)
    W=observable_weights(op,cut)
    out.mkdir(exist_ok=True,parents=True)
    u=pack(sh)
    records=[]
    for a,b in zip(coarse["rows"],coarse["rows"][1:]):
        h=(b["time"]-a["time"])/6
        trace=[]
        for j in range(6):
            before=stable_total(op,u)
            base=W@op.moment(u)[0]/rho
            def guard(v,s,s0):
                stable=float(stable_total(op,v)-before)
                tested=[]
                tol=float(16*np.finfo(np.longdouble).eps*abs(before))
                if stable<=tol:
                    raise ValueError(f"Stabilized entropy failed: {stable=:g}, {tol=:g}")
                for n in (8,16,32,64):
                    z, gross, _=entropy_line(op,u,v,n)
                    tested.append([n,z])
                    tolerance=max(2e-5*stable,tol)
                    if len(tested)>1 and min(z,tested[-2][1])>tol and max(abs(z-stable),abs(z-tested[-2][1]))<tolerance:
                        trace.append(dict(legacy=float(s-s0),stable=stable,line=tested,
                                          total_subtraction_resolution=tol,gross=float(gross)))
                        return True
                raise ValueError(f"Entropy path not certified: {tested!r}, {stable=:g}")
            u,info=advance(op,u,h,S["T_visible_GeV"],initial_theta=predictor(op,u,h),
                           vacuum_scaled=True,balanced_merit=True,entropy_guard=guard)
            current=W@op.collision(u)[0]/rho
            w=W@op.moment(u)[0]/rho
            assert abs(w-base-h*current)<4e-10*max(w,base,1e-70)
            trace[-1].update(solver=info,W_over_rho=float(w),current=float(current))
        Q=dict(S,**{k:x.tolist() for k,x in zip(keys,op.arrays(u),strict=True)})
        assert abs(sum(excess(Q,cut))/rho/w-1.)<4e-11
        n,_,e=get_mom(p,op.arrays(u))
        nd=float(max(abs(n.sum(axis=1)/N0.sum(axis=1)-1.)))
        ed=float(abs(e.sum()/rho-1.))
        assert max(nd,ed)<8e-12
        angular=integrate(Q,cut,8,True)/(rho*H)
        assert abs(angular[3]-current)<3e-5*max(angular[1:3])
        fit=dict(time=b["time"],W=float(w),current=float(current),angular=float(angular[3]),
            discrete_N_E_relative_defects=[nd,ed],snapshot_sha256=digest(Q),substeps=trace,
            coarse_W=b["tail_energy"],coarse_current=b["signed_native_flux"])
        records.append(fit)
        atomic(out/f"step-{len(records)}.json",dict(snapshot=Q,row=fit))
        print(b["time"],w,current,flush=True)
    atomic(out/"summary.json",dict(scope="Whole fixed prompt-only clock; independent entropy certification. "
           "No relaxed numerical root, occupation floor, or temporal subdivision fallback.",
           seed_sha256=coarse["plan"]["seed_sha256"],rows=records))


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source",type=Path)
    parser.add_argument("out",type=Path)
    args=parser.parse_args()
    run(args.source,args.out)
