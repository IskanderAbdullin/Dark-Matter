"""Check common-pair order at an explicitly locked prompt/elastic time grid."""

import argparse
import json
from pathlib import Path

import numpy as np

import coupled_warm_pair_continuation_reference as kernel
from joint_resonant_centered_stiff_reference import advance
from joint_resonant_frw_split_reference import pack,redshift
from joint_visible_heavy_adiabatic_reference import snapshot_functions
from warm_coupled_serial_restart_reference import atomic,sha
from warm_history_lower_matched_audit_reference import measures
from warm_light_occupation_stage_audit_reference import digest,read
from warm_pair_common_fermion_chart_reference import rate


def run(eos,parent,reference,guide,output,order,n):
    A,B,C=map(read,(parent,reference,guide))
    S,U,G=[v["snapshot"] for v in (A,B,C)]
    assert A["snapshot_sha256"]==B["parent_snapshot_sha256"]==C["source"]["parent_snapshot_sha256"]
    assert B["snapshot_sha256"]==C["source"]["preconditioner_sha256"]
    assert C["source"]["fine"]==(order==(4,10,10,4))
    assert all(r["snapshot_sha256"]==digest(r["snapshot"]) for r in (A,B,C))
    b=G["scale_factor"]/S["scale_factor"]
    q=np.asarray(S["lab_momentum_edges_GeV"])
    np.testing.assert_allclose(q/b,G["lab_momentum_edges_GeV"],rtol=3e-14)
    p=b**.5
    end=[np.asarray(G["entropy_variable_"+k]) for k in ("scalar","fermion")]
    # Reverse passive remap only to make a seed. The output did
    # include the second pair kick, so this is not yet the root.
    near=pack(redshift(q,b,p,end)[0])
    src=dict(parent_sha256=A["snapshot_sha256"],nine_cell_sha256=B["snapshot_sha256"],
             guide_sha256=C["snapshot_sha256"],eos_sha256=sha(eos),selected_substeps=n,
             order=list(order),chart_sha256=sha(Path(rate.__code__.co_filename)),
             trajectory_wrapper_sha256=sha(Path(kernel.__file__)))
    if output.exists():
        d=read(output)
        assert d["source"]==src and digest(d["snapshot"])==d["snapshot_sha256"]
        return {k:v for k,v in d.items() if k!="snapshot"}
    status=[]
    old_r,old_s=kernel.rate,kernel.collide

    def whole(p,z,unused):
        tmp=dict(lab_momentum_edges_GeV=p,entropy_variable_scalar=z[0],entropy_variable_fermion=z[1])
        B,_,F,_=snapshot_functions(tmp)
        v,_,I,e,_=rate(p,F,B,order)
        return v,e,I[0,0]/I[1,0]

    def collide(model,u,h,T):
        start=u.copy()
        for j in range(1,n+1):
            x=j/n
            guesses=[(1-x)*start+x*near,u.copy(),near]
            faults=[]
            for phase,g in enumerate(guesses):
                try:
                    end,res=advance(model,u,h/n,T,initial_theta=g,vacuum_scaled=True,balanced_merit=True)
                    break
                except (ValueError,AssertionError,RuntimeError,RuntimeWarning,FloatingPointError) as exc:
                    faults.append(str(exc))
            else:
                raise RuntimeError(f"Cannot accept locked root {j}/{n}: {faults}")
            status.append(dict(index=j,accepted_seed=phase,rejected_seed_reasons=faults,
                               step_H=h/n,**res))
            u=end
        return u,max(s["residual"] for s in status),n,0

    X={}
    try:
        kernel.rate,kernel.collide=whole,collide
        result=kernel.run(eos,S,1,b,checkpoint=X,track_exterior=True,local_prompt_jacobian=True,
                          scheduled_single_block=True)
    finally:
        kernel.rate,kernel.collide=old_r,old_s
    assert result["acceptances"]==n and len(status)==n
    assert max(abs(result[k]) for k in ("energy_ledger","count_ledger"))<3e-13
    D=measures(X);C0=measures(U)
    def occ(s,z):
        return z[0,1,0]/(4*s["lab_momentum_edges_GeV"][1]**3/(6*np.pi**2))
    out=dict(source=src,snapshot=X,snapshot_sha256=digest(X),logged_stiff_substeps=status,
             cold_slope_per_ln_a=(occ(X,D)-occ(S,measures(S)))/np.log(b),
             moment_L1_against_nine_cell=(abs(D-C0).sum(axis=-1)/C0.sum(axis=-1)).tolist(),
             energy_and_count_ledger=[result["energy_ledger"],result["count_ledger"]])
    output.parent.mkdir(parents=True,exist_ok=True)
    atomic(output,out)
    return {k:v for k,v in out.items() if k!="snapshot"}


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("eos","parent","nine","guide","output"):
        p.add_argument(k,type=Path)
    p.add_argument("--fine",action="store_true")
    p.add_argument("--substeps",default=2,type=int)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.parent,a.nine,a.guide,a.output,
                         (4,10,10,4) if a.fine else (3,7,6,2),a.substeps),indent=2,allow_nan=False))
