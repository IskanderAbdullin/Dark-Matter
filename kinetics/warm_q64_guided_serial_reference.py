"""Continue a physical q64 state using q32 endpoint shapes only as Newton guesses."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

import coupled_warm_pair_continuation_reference as kernel
from joint_resonant_centered_stiff_reference import advance
from joint_resonant_frw_split_reference import pack
from warm_coupled_serial_restart_reference import atomic,bytes_of,sha
from warm_light_occupation_stage_audit_reference import digest,read
from warm_prompt_reentry_shape_audit_reference import get_mom


def run(eos,start,guides,out):
    A=read(start)
    S=A["snapshot"]
    assert A["snapshot_sha256"]==digest(S)
    tape=read(guides)
    actions=[v["checkpoint"] for v in tape["blocks"]]
    out.mkdir(parents=True,exist_ok=True)
    plan=dict(source_sha256=sha(start),guides_sha256=sha(guides),eos_sha256=sha(eos),
              kernel_sha256=sha(Path(kernel.__file__)),solver_wrapper_sha256=sha(Path(__file__)),
              guide_snapshots=[x["snapshot_sha256"] for x in actions])
    old_collider=kernel.collide
    prev=S
    results=[]
    for i,g in enumerate(actions,1):
        guide=g["snapshot"]
        assert g["snapshot_sha256"]==digest(guide)
        path=out/f"step-{i}.json"
        if path.exists():
            d=read(path)
            assert d["plan"]==plan and d["step"]==i and d["parent_snapshot_sha256"]==digest(prev)
            assert digest(d["snapshot"])==d["snapshot_sha256"]
            prev=d["snapshot"]
            results.append(d["audit"])
            continue
        b=guide["scale_factor"]/prev["scale_factor"]
        before=np.asarray(prev["lab_momentum_edges_GeV"])
        narrow=np.asarray(guide["lab_momentum_edges_GeV"])
        np.testing.assert_allclose(before[:len(narrow)]/b,narrow,rtol=8e-14,atol=0.)
        shape=[np.asarray(guide["entropy_variable_"+j]) for j in ("scalar","fermion")]
        trials=[]

        def solver(model,u,h,T):
            assert model.prompt
            trial=[a.copy() for a in model.arrays(u)]
            for a,s in zip(trial,shape,strict=True):
                a[:,:s.shape[1]]=s
            obs=[]
            def observe(i,r,j):
                obs.append((int(i),float(r),int(j)))
            try:
                v,data=advance(model,u,h,T,reference_guess=False,vacuum_scaled=True,
                               initial_theta=pack(trial),balanced_merit=True,observer=observe)
            except (RuntimeError,FloatingPointError,RuntimeWarning,ValueError,AssertionError) as e:
                trials.append(dict(used_q32_guess=False,fallback_error=str(e),observations=obs))
                return old_collider(model,u,h,T,continuation=True,balanced=True)
            trials.append(dict(used_q32_guess=True,observations=obs,**data))
            return v,data["residual"],1,0

        ck={}
        kernel.collide=solver
        try:
            r=kernel.run(eos,prev,1,b,checkpoint=ck,track_exterior=True,
                         local_prompt_jacobian=True,scheduled_single_block=True)
        finally:
            kernel.collide=old_collider
        assert len(trials)==1
        assert abs(ck["T_visible_GeV"]/guide["T_visible_GeV"]-1)<1e-13
        def mom(s):
            return get_mom(np.asarray(s["lab_momentum_edges_GeV"]),
                           [s["entropy_variable_scalar"],s["entropy_variable_fermion"]])
        N,K,E=mom(ck); aN,aK,aE=mom(guide)
        d=len(narrow)-1
        v=dict(T_GeV=ck["T_visible_GeV"],scale_factor=ck["scale_factor"],trial=trials[0],
               common_grid_count_L1=(abs(N[:,:d]-aN).sum(axis=1)/aN.sum(axis=1)).tolist(),
               common_grid_full_energy_L1_over_total=(abs(E[:,:d]-aE).sum(axis=1)/aE.sum()).tolist(),
               final_added_region_count_fraction=(N[:,d:].sum(axis=1)/N.sum(axis=1)).tolist(),
               full_energy_and_count_defects=[r["energy_ledger"],r["count_ledger"]])
        state=dict(plan=plan,step=i,parent_snapshot_sha256=digest(prev),snapshot=ck,
                   snapshot_sha256=digest(ck),audit=v,operator_summary=r)
        atomic(path,state)
        results.append(v)
        prev=ck
        print(json.dumps({"finished_step":i,"T_GeV":ck["T_visible_GeV"],
                          "guided":trials[0]["used_q32_guess"]}),file=__import__("sys").stderr,flush=True)
    return dict(scope="Archived q32 end profiles are solver seeds only; the accepted trajectory "
                      "uses the full original q64 operators and physical moment targets. "
                      "Old exterior cohorts remain collisionless. Not a full cosmology.",
                plan=plan,steps=results,final_snapshot_sha256=digest(prev))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for n in ("eos","start","guides","out"):p.add_argument(n,type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.start,a.guides,a.out),indent=2,allow_nan=False))
