"""Portable wider continuation using its own two accepted shapes as nonlinear seeds."""

import argparse
import json
import sys
from pathlib import Path

import numpy as np

import coupled_warm_pair_continuation_reference as kernel
from frw_exterior_restart_reference import import_restart
from joint_resonant_centered_stiff_reference import advance
from joint_resonant_frw_split_reference import pack
from warm_coupled_serial_restart_reference import atomic,sha
from warm_light_occupation_stage_audit_reference import digest,read
from warm_prompt_reentry_shape_audit_reference import get_mom


def run(eos,tape,b,n,out):
    seed=read(tape)
    P,S=[v["snapshot"] for v in seed[-2:]]
    assert all(digest(v["snapshot"])==v["snapshot_sha256"] for v in seed[-2:])
    assert n>=1 and b>1.
    out.mkdir(parents=True,exist_ok=True)
    plan=dict(tape_sha256=sha(tape),EOS_sha256=sha(eos),wrapper_sha256=sha(Path(__file__)),
              kernel_sha256=sha(Path(kernel.__file__)),step_expansion=b,steps=n)
    a0=S["scale_factor"]
    old=import_restart(S["positive_exterior_restart"]).array()
    nold=old.shape[1]
    def mom(s):
        return get_mom(np.asarray(s["lab_momentum_edges_GeV"]),
                       [s["entropy_variable_scalar"],s["entropy_variable_fermion"]])
    initial=mom(S)[0].sum()*a0**3
    keys=("entropy_variable_scalar","entropy_variable_fermion")
    row=[]
    reset=kernel.collide
    for j in range(1,n+1):
        path=out/f"step-{j}.json"
        if path.exists():
            Z=read(path)
            assert Z["plan"]==plan and Z["index"]==j and Z["previous"]==digest(S)
            assert Z["snapshot_sha256"]==digest(Z["snapshot"])
            P,S=S,Z["snapshot"]
            row.append(Z["audit"])
            continue
        ratio=np.log(b)/np.log(S["scale_factor"]/P["scale_factor"])
        prev=[np.asarray(S[k]) for k in keys]
        p=[u+ratio*(u-np.asarray(P[k])) for u,k in zip(prev,keys,strict=True)]
        if np.any(p[0][0]+np.minimum(0.,p[0][1])<=0.):
            p=prev
        stats={}
        hist=[]

        def solve(model,u,h,T):
            def obs(i,r,k):hist.append((i,r,k))
            try:
                N,payload=advance(model,u,h,T,reference_guess=False,vacuum_scaled=True,
                                  initial_theta=pack(p),balanced_merit=True,observer=obs)
            except (RuntimeError,FloatingPointError,RuntimeWarning,ValueError,AssertionError) as e:
                stats.update(used_history=False,history_trial_error=str(e),observations=hist)
                return reset(model,u,h,T,continuation=True,balanced=True)
            stats.update(used_history=True,observations=hist,**payload)
            return N,payload["residual"],1,0

        ck={}
        kernel.collide=solve
        try:
            r=kernel.run(eos,S,1,b,checkpoint=ck,track_exterior=True,local_prompt_jacobian=True,
                         scheduled_single_block=True)
        finally:
            kernel.collide=reset
        now=ck["scale_factor"]
        C=import_restart(ck["positive_exterior_restart"]).array()
        D=C[:,:nold].copy()
        D[[0,5]]*=np.array([now/a0,(now/a0)**3])[:,None]
        np.testing.assert_allclose(D,old,rtol=9e-14,atol=0.)
        q=abs((mom(ck)[0].sum()+C[5,nold:].sum())*now**3/initial-1.)
        assert q<2e-11
        a=dict(step=j,T_GeV=ck["T_visible_GeV"],scale_factor=now,solver=stats,
               local_energy_defect=r["energy_ledger"],cumulative_number_residual=float(q))
        atomic(path,dict(plan=plan,index=j,previous=digest(S),snapshot=ck,
                         snapshot_sha256=digest(ck),audit=a,operator_summary=r))
        P,S=S,ck
        row.append(a)
        print(json.dumps(dict(step=j,T_GeV=S["T_visible_GeV"],history=stats["used_history"])),file=sys.stderr,flush=True)
    return dict(scope="Accepted q64 equations use no narrower future distribution, one bounded "
                      "history-based numerical initialization and ordinary wider fallback.",
                plan=plan,steps=row,final_snapshot_sha256=digest(S))


if __name__=="__main__":
    q=argparse.ArgumentParser(description=__doc__)
    for a in ("eos","tape","out"):q.add_argument(a,type=Path)
    q.add_argument("--step-expansion",type=float,required=True)
    q.add_argument("--blocks",type=int,required=True)
    a=q.parse_args()
    print(json.dumps(run(a.eos,a.tape,a.step_expansion,a.blocks,a.out),indent=2,allow_nan=False))
