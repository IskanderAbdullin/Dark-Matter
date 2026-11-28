"""Byte-exact value and tangent parity for stored versus streamed angular panels."""

import argparse
import json
from pathlib import Path

import numpy as np

from compact_prompt_jacobian_reference import compile_prompt
from joint_resonant_frw_split_reference import pack
from warm_compact_absolute_prompt_reference import LocalAbsolutePrompt
from warm_comoving_parent_chart_reference import ComovingCollision,common_rule
from warm_late_finite_parent_prompt_reference import FiniteParent
from warm_light_occupation_stage_audit_reference import digest,read


def audit(path):
    x=read(path)
    S=x["snapshot"]
    assert x["hash"]==digest(S)
    q=np.asarray(S["lab_momentum_edges_GeV"])
    Q,w=common_rule(q,(1.,),2)
    k=ComovingCollision(q,1.,Q,w)
    u=pack([S["entropy_variable_"+name] for name in ("scalar","fermion")])
    d=np.sin(.71*np.arange(len(u))+.41)
    rows=[]
    for finite in (False,True):
        model=FiniteParent(k,2.6e-23) if finite else LocalAbsolutePrompt(k,2.6e-23,[],True,stream_local=True)
        if finite:
            model.old=model.equilibrium(u)
            model.dt=1.5e-8
        a,j=model.collision(u)
        jj=j(d)
        model.parents=compile_prompt(model)
        b,i=model.collision(u)
        ii=i(d)
        assert a.tobytes()==b.tobytes() and jj.tobytes()==ii.tobytes()
        rows.append(dict(finite=finite,value_and_derivative_identical=True,parent_panels=len(model.parents)))
    return dict(input_hash=x["hash"],angular_nodes=len(k.y),rows=rows)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("state",type=Path)
    a=p.parse_args()
    print(json.dumps(audit(a.state),indent=2,allow_nan=False))
