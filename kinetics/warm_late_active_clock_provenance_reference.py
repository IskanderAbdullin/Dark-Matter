"""Verify selected active checkpoints against a fixed three-root clock."""

import argparse
import json
from pathlib import Path

import numpy as np

from joint_prompt_two_moment_evolution_reference import Collision
from joint_resonant_centered_stiff_reference import advance
from joint_resonant_frw_split_reference import pack
from late_radiation_background_reference import MPL,TabulatedRadiation
from warm_compact_absolute_prompt_reference import LocalAbsolutePrompt
from warm_corrected_late_prompt_cut_aligned_reference import predictor
from warm_light_occupation_stage_audit_reference import digest,read
from warm_prompt_reentry_shape_audit_reference import get_mom


def run(eos,D,indices):
    initial=read(D/"step-0.json")
    S=initial["snapshot"]
    assert initial["hash"]==digest(S)
    assert initial["plan"]["substeps"]==3 and initial["plan"]["artificial_empty_shell_filling"]==0.
    p=np.array(S["lab_momentum_edges_GeV"])
    keys=("entropy_variable_scalar","entropy_variable_fermion")
    rho=get_mom(p,[S[k] for k in keys])[2].sum()
    H=np.sqrt((TabulatedRadiation(eos,S["T_visible_GeV"]).state(1.)[1]+rho)/3)/MPL
    op=LocalAbsolutePrompt(Collision(p),H,[],True)
    out=[]
    for i in indices:
        A=read(D/f"step-{i-1}.json")
        B=read(D/f"step-{i}.json")
        assert A["plan"]==B["plan"]==initial["plan"]
        for row in (A,B):
            assert row["hash"]==digest(row["snapshot"])
            assert row["snapshot"]["lab_momentum_edges_GeV"]==S["lab_momentum_edges_GeV"]
        a=pack([A["snapshot"][k] for k in keys])
        b=pack([B["snapshot"][k] for k in keys])
        dt=(B["time"]-A["time"])/3
        info=[]
        for _ in range(3):
            guess=predictor(op,a,dt)
            a,res=advance(op,a,dt,S["T_visible_GeV"],initial_theta=guess,
                          vacuum_scaled=True,balanced_merit=True)
            info.append(res)
        x=op.moment(a)[0].reshape(2,op.n,2)
        y=op.moment(b)[0].reshape(2,op.n,2)
        scale=y[:,0:op.n,0].sum(axis=-1)
        discrepancy=(abs(x-y).sum(axis=1)/scale[:,None]).tolist()
        out.append(dict(interval=i,begin=A["time"],end=B["time"],source_sha256=A["hash"],
                        end_sha256=B["hash"],three_microstep_size_H=dt,residuals=[s["residual"] for s in info],
                        archived_equivalent_normalized_moment_l1=discrepancy))
    return dict(scope="Independent locked three-microstep replay of selected active intervals, "
                      "not a new time grid or an assumption on untested adaptive paths.",
                active_seed_sha256=initial["hash"],hard_GeV=initial["plan"]["hard_GeV"],controls=out)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("eos",type=Path);p.add_argument("checkpoints",type=Path)
    p.add_argument("intervals",type=int,nargs="+")
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.checkpoints,a.intervals),indent=2,allow_nan=False))
