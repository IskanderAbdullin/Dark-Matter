"""Replace the full pair quadrature on a late corrected q64 interval."""

import argparse
import json
from pathlib import Path

import numpy as np

import coupled_warm_pair_continuation_reference as kernel
from joint_resonant_centered_stiff_reference import advance
from joint_resonant_frw_split_reference import pack
from joint_visible_heavy_adiabatic_reference import snapshot_functions
from warm_coupled_serial_restart_reference import atomic, sha
from warm_history_lower_matched_audit_reference import measures
from warm_light_occupation_stage_audit_reference import digest, read
from warm_pair_common_fermion_chart_reference import rate as common


def replay(eos, start, endpoint, output, fine):
    A,B=map(read,(start,endpoint))
    S,G=A["snapshot"],B["snapshot"]
    assert A["snapshot_sha256"]==digest(S)
    assert B["snapshot_sha256"]==digest(G)
    assert A["snapshot_sha256"]==B["parent_snapshot_sha256"] and A["step"]+1==B["step"]
    assert A["plan"]==B["plan"] and B["plan"]["eos_sha256"]==sha(eos)
    b=G["scale_factor"]/S["scale_factor"]
    assert 1.<b<1.01
    p0=np.asarray(S["lab_momentum_edges_GeV"])
    np.testing.assert_allclose(p0/b,G["lab_momentum_edges_GeV"],rtol=8e-14,atol=0.)
    target_theta=pack([np.asarray(G["entropy_variable_"+k]) for k in ("scalar","fermion")])
    src=dict(parent_snapshot_sha256=A["snapshot_sha256"],preconditioner_sha256=B["snapshot_sha256"],
             eos_sha256=sha(eos),whole_chart_sha256=sha(Path(common.__code__.co_filename)),fine=fine,
             stage=B["step"], relative_expansion=b)
    if output.exists():
        R=read(output)
        assert src==R["source"] and R["snapshot_sha256"]==digest(R["snapshot"])
        return {k:v for k,v in R.items() if k!="snapshot"}

    observations, charts, ck=[],[],{}
    solver_stats={}
    old_rate,old_collide=kernel.rate,kernel.collide

    def pair(p,theta,order):
        target=dict(lab_momentum_edges_GeV=p,entropy_variable_scalar=theta[0],entropy_variable_fermion=theta[1])
        B,_,F,_=snapshot_functions(target)
        use=(4,10,10,4) if fine else order
        q,g,I,d,hi=common(p,F,B,use)
        charts.append(dict(order=list(use),endpoint_pair_defect=d,max_outside_momentum_GeV=hi,
                           inverse_exterior_light_count_fraction=float(g[0,1,0,-1]/(2*I[0,0]))))
        return q,d,I[0,0]/I[1,0]

    def obs(i,r,j):
        observations.append([int(i),float(r),int(j)])

    def collide(model,u,h,T):
        try:
            V,data=advance(model,u,h,T,reference_guess=False,vacuum_scaled=True,
                           initial_theta=target_theta,balanced_merit=True,observer=obs)
        except (ValueError,AssertionError,RuntimeWarning,FloatingPointError,RuntimeError) as e:
            solver_stats.update(used_nine_cell_guess=False,fallback_reason=str(e))
            return old_collide(model,u,h,T,continuation=True,balanced=True)
        solver_stats.update(used_nine_cell_guess=True,**data)
        return V,data["residual"],1,0

    try:
        kernel.rate,kernel.collide=pair,collide
        R=kernel.run(eos,S,1,b,checkpoint=ck,track_exterior=True,
                     local_prompt_jacobian=True,scheduled_single_block=True)
    finally:
        kernel.rate,kernel.collide=old_rate,old_collide
    assert len(charts)>0
    assert abs(ck["T_visible_GeV"]/G["T_visible_GeV"]-1.)<4e-14
    assert max(abs(R[k]) for k in ("energy_ledger","count_ledger"))<5e-13
    def cold(s):
        return measures(s)[0,1,0]/(4*s["lab_momentum_edges_GeV"][1]**3/(6*np.pi**2))
    X,Y=measures(ck),measures(G)
    Z=dict(source=src,snapshot=ck,snapshot_sha256=digest(ck),T_GeV=ck["T_visible_GeV"],
           integrated_ledgers=[R["energy_ledger"],R["count_ledger"]],
           cold_slopes_nine_cell_and_common=[(cold(s)-cold(S))/np.log(b) for s in (G,ck)],
           moment_L1_against_nine_cell=(abs(X-Y).sum(axis=-1)/Y.sum(axis=-1)).tolist(),
           solver=solver_stats,observations=observations,charts=charts)
    output.parent.mkdir(parents=True,exist_ok=True)
    atomic(output,Z)
    return {k:v for k,v in Z.items() if k!="snapshot"}


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("eos","start","endpoint","output"):
        p.add_argument(k,type=Path)
    p.add_argument("--fine",action="store_true")
    a=p.parse_args()
    print(json.dumps(replay(a.eos,a.start,a.endpoint,a.output,a.fine),indent=2,allow_nan=False))
