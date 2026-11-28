"""Compare the whole common-coordinate pair kernel on one actual q32 block."""

import argparse
import json
from pathlib import Path

import numpy as np

import coupled_warm_pair_continuation_reference as kernel
from joint_visible_heavy_adiabatic_reference import snapshot_functions
from warm_coupled_serial_restart_reference import atomic, sha
from warm_history_lower_matched_audit_reference import measures
from warm_light_occupation_stage_audit_reference import digest, read
from warm_pair_common_fermion_chart_reference import rate
from warm_q64_history_serial_reference import run


def block(eos, prev, start, original, nine, folder, chosen=None):
    A,S,X,Y=[read(x) for x in (prev,start,original,nine)]
    for row in (A,S,X,Y):
        assert digest(row["snapshot"]) == row["snapshot_sha256"]
    assert A["snapshot_sha256"]==S["previous"]
    assert S["snapshot_sha256"]==X["previous"]==Y["previous"]
    b=X["plan"]["step_expansion"]
    assert X["plan"]["EOS_sha256"]==sha(eos)
    folder.mkdir(parents=True,exist_ok=True)
    seed=folder/"seed.json"
    manifest=folder/"common_chart_manifest.json"
    plan=dict(quadrature_source_sha256=sha(Path(rate.__code__.co_filename)),
              previous_sha256=A["snapshot_sha256"],start_sha256=S["snapshot_sha256"],
              control_native_sha256=X["snapshot_sha256"],control_nine_sha256=Y["snapshot_sha256"],
              eos_sha256=sha(eos),outer_expansion=b)
    if chosen is not None:
        plan["fixed_pair_order"]=list(chosen)
    for path, value in ((seed,[A,S]),(manifest,plan)):
        if path.exists():
            assert read(path)==value
        else:
            assert not list(folder.glob("step-*.json"))
            atomic(path,value)
    old=kernel.rate
    audit=[]

    def common(p, th, order):
        snap=dict(lab_momentum_edges_GeV=p,entropy_variable_scalar=th[0],entropy_variable_fermion=th[1])
        B,_,F,_=snapshot_functions(snap)
        order = order if chosen is None else chosen
        net,pro,r,err,plim=rate(p,F,B,order)
        audit.append(dict(order=order, count_and_energy_defect=err, maximum_exterior_p=plim,
                          inverse_fermions_outside_fraction=float(pro[0,1,0,-1]/(2*r[0,0]))))
        return net,err,r[0,0]/r[1,0]

    try:
        kernel.rate=common
        Q=run(eos,seed,b,1,folder)
    finally:
        kernel.rate=old
    R=read(folder/"step-1.json")
    candidates=[v["snapshot"] for v in (X,Y,R)]
    assert Q["final_snapshot_sha256"]==digest(candidates[-1])==R["snapshot_sha256"]
    def occ(t):return measures(t)[0,1,0]/(4*t["lab_momentum_edges_GeV"][1]**3/(6*np.pi**2))
    u=occ(S["snapshot"])
    slopes=[(occ(s)-u)/np.log(b) for s in candidates]
    lhs=[measures(s) for s in candidates]
    L=[(np.abs(lhs[2]-j).sum(axis=-1)/j.sum(axis=-1)).tolist() for j in lhs[:2]]
    return dict(scope="Selected one-block comparison from identical measured inputs, not a converged full trajectory.",
                T_GeV=candidates[-1]["T_visible_GeV"], native_lifted_common_slopes=slopes,
                shell_L1_against_native_and_nine_cell=L, controls=audit,
                final_snapshot_sha256=R["snapshot_sha256"],energy_defect=R["operator_summary"]["energy_ledger"],
                count_defect=R["operator_summary"]["count_ledger"],manifest=plan)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for x in ("eos","prev","start","original","nine","folder"):
        p.add_argument(x,type=Path)
    p.add_argument("--fine",action="store_true")
    a=p.parse_args()
    print(json.dumps(block(a.eos,a.prev,a.start,a.original,a.nine,a.folder,(4,10,10,4) if a.fine else None),indent=2,allow_nan=False))
