"""Test a history-only q64 Newton initialization against a narrow-guided control."""

import argparse
import json
from pathlib import Path

import numpy as np

import coupled_warm_pair_continuation_reference as kernel
from joint_resonant_centered_stiff_reference import advance
from joint_resonant_frw_split_reference import pack
from warm_coupled_serial_restart_reference import atomic,bytes_of
from warm_light_occupation_stage_audit_reference import digest,read
from warm_prompt_reentry_shape_audit_reference import get_mom
from warm_q64_guided_serial_reference import run as guided_run


def single(eos,tape,narrow,work):
    W=read(tape)
    assert len(W)>=2
    P,S=[x["snapshot"] for x in W[-2:]]
    assert all(digest(x["snapshot"])==x["snapshot_sha256"] for x in W[-2:])
    g=read(narrow)
    target=g["snapshot"]
    assert g["snapshot_sha256"]==digest(target)
    b=target["scale_factor"]/S["scale_factor"]
    old=S["scale_factor"]/P["scale_factor"]
    ratio=np.log(b)/np.log(old)
    assert .5<ratio<1.5
    work.mkdir(parents=True,exist_ok=True)
    key=("entropy_variable_scalar","entropy_variable_fermion")
    pred=[np.asarray(S[k])+ratio*(np.asarray(S[k])-np.asarray(P[k])) for k in key]
    assert np.min(pred[0][0]+np.minimum(0,pred[0][1]))>0.
    obs=[]
    stats={}
    old_solver=kernel.collide

    def solver(model,u,h,T):
        assert len(obs)==0 and model.prompt
        def observe(i,r,j):obs.append([int(i),float(r),int(j)])
        try:
            v,d=advance(model,u,h,T,reference_guess=False,vacuum_scaled=True,
                        initial_theta=pack(pred),balanced_merit=True,observer=observe)
            stats.update(d,history_only_accepted=True)
            return v,d["residual"],1,0
        except (RuntimeError,FloatingPointError,RuntimeWarning,ValueError,AssertionError) as e:
            stats.update(history_only_accepted=False,fallback_error=str(e))
            return old_solver(model,u,h,T,continuation=True,balanced=True)

    ck={}
    kernel.collide=solver
    try:
        v=kernel.run(eos,S,1,b,checkpoint=ck,track_exterior=True,local_prompt_jacobian=True,scheduled_single_block=True)
    finally:
        kernel.collide=old_solver
    assert stats["history_only_accepted"]
    assert abs(ck["T_visible_GeV"]/target["T_visible_GeV"]-1)<2e-13
    atomic(work/"history-endpoint.json",dict(snapshot=ck,snapshot_sha256=digest(ck),
               parent_snapshot_sha256=digest(S),nonlinear_seed_snapshot_sha256=digest(P),
               solver=stats,observations=obs,kernel=v))
    # Use the actual independent q32 end to check the same wider equation,
    # but do not make that shape available to the history-only solver.
    orig=work/"wide-source.json"
    guides=work/"narrow-actual.json"
    orig.write_bytes(bytes_of(W[-1]))
    guides.write_bytes(bytes_of(dict(blocks=[dict(checkpoint=g)])))
    result=guided_run(eos,orig,guides,work/"narrow-guided")
    other=read(work/"narrow-guided"/"step-1.json")["snapshot"]
    n=[]
    for s in (ck,other):
        p=np.asarray(s["lab_momentum_edges_GeV"])
        N,K,E=get_mom(p,[s[k] for k in key])
        n.append(np.asarray([N,E]))
    d=np.sum(abs(n[0]-n[1]),axis=2)
    d/=np.vstack((np.sum(n[1][0],axis=1),np.repeat(n[1][1].sum(),2)))
    assert np.max(d)<5e-10
    return dict(scope="The q32 guide is unavailable to the history predictor's solve; it sets "
                      "only the matching interval and provides a separate wider control. "
                      "All accepted targets and collision operators remain q64. Valid for this tested block only.",
                final_T_GeV=ck["T_visible_GeV"],history_seed=stats,observations=obs,
                narrow_seed=result["steps"][0]["trial"],
                retained_number_and_energy_L1_between_solutions=d.tolist(),
                history_sha256=digest(ck),narrow_seed_control_sha256=digest(other))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for a in ("eos","tape","q32_end","work"):p.add_argument(a,type=Path)
    x=p.parse_args()
    print(json.dumps(single(x.eos,x.tape,x.q32_end,x.work),indent=2,allow_nan=False))
