"""Subdivide matched q64 blocks with interpolation used solely for Newton seeds."""

import argparse
import json
from pathlib import Path

import numpy as np

from late_radiation_background_reference import TabulatedRadiation
from warm_coupled_serial_restart_reference import bytes_of
from warm_light_occupation_stage_audit_reference import digest,read
from warm_prompt_reentry_shape_audit_reference import get_mom
from warm_q64_guided_serial_reference import run as guided_run


def make_guides(eos,anchor,observed,num,path):
    c0=read(anchor)
    left=(c0[-1] if isinstance(c0,list) else c0)["snapshot"]
    targets=[x["checkpoint"] for x in read(observed)["blocks"]]
    blocks=[]
    for c in targets:
        right=c["snapshot"]
        assert digest(right)==c["snapshot_sha256"]
        T=TabulatedRadiation(eos,left["T_visible_GeV"])
        old=left["scale_factor"]; later=right["scale_factor"]
        ratio=later/old
        for i in range(1,num+1):
            if i==num:
                blocks.append(dict(checkpoint=c,kind="archived end"))
                continue
            x=i/num
            b=ratio**x
            S=dict(lab_momentum_edges_GeV=(np.asarray(left["lab_momentum_edges_GeV"])/b).tolist(),
                   scale_factor=old*b,T_visible_GeV=T.state(b)[0])
            for kind in ("scalar","fermion"):
                k="entropy_variable_"+kind
                a,bb=np.asarray(left[k]),np.asarray(right[k])
                S[k]=((1-x)*a+x*bb).tolist()
            blocks.append(dict(checkpoint=dict(snapshot=S,snapshot_sha256=digest(S)),kind="interpolated guess only"))
        left=right
    b=bytes_of(dict(blocks=blocks))
    if path.exists():
        assert path.read_bytes()==b
    else:
        path.write_bytes(b)
    return blocks


def end_moment(s):
    p=np.asarray(s["lab_momentum_edges_GeV"])
    N,K,E=get_mom(p,[s["entropy_variable_scalar"],s["entropy_variable_fermion"]])
    return np.asarray([N,E])*s["scale_factor"]**3


def run(eos,A,q32_start,guides,base,work,splits):
    work.mkdir(parents=True,exist_ok=True)
    g=work/"constructed-seeds.json"
    make_guides(eos,q32_start,guides,splits,g)
    result=guided_run(eos,A,g,work/"run")
    a=read(base)[-1]["snapshot"]
    b=read(work/"run"/f"step-{len(result['steps'])}.json")["snapshot"]
    assert abs(a["scale_factor"]/b["scale_factor"]-1)<4e-14
    np.testing.assert_allclose(a["lab_momentum_edges_GeV"],b["lab_momentum_edges_GeV"],rtol=8e-14,atol=0.)
    ref=end_moment(a); end=end_moment(b); normal=np.vstack((ref[0].sum(axis=1),np.repeat(ref[1].sum(),2)))
    return dict(scope="Same q64 initial moments/operators and final physical time. Interpolated "
                      "q32 entropy coordinates are numerical seeds only, not physical data. "
                      "Compares selected-operator finite-grid trajectories.",
                subdivisions_per_archived_block=splits,
                final_L1_over_original_q64=[(abs(ref-end).sum(axis=2)/normal).tolist()],
                guide_result=result)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for a in ("eos","q64","q32","guides","base","work"):p.add_argument(a,type=Path)
    p.add_argument("--splits",type=int,choices=(2,4),default=2)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.q64,a.q32,a.guides,a.base,a.work,a.splits),indent=2,allow_nan=False))
