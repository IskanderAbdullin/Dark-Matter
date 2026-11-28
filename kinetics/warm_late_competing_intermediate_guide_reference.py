"""Generate an intermediate prompt-only Newton guide on a late frozen clock."""

import argparse
from pathlib import Path

import numpy as np

from joint_prompt_two_moment_evolution_reference import Collision
from joint_resonant_centered_stiff_reference import advance
from joint_resonant_frw_split_reference import pack
from late_radiation_background_reference import MPL, TabulatedRadiation
from warm_compact_absolute_prompt_reference import LocalAbsolutePrompt
from warm_corrected_late_prompt_cut_aligned_reference import predictor
from warm_coupled_serial_restart_reference import atomic, sha
from warm_light_occupation_stage_audit_reference import digest, read
from warm_prompt_reentry_shape_audit_reference import get_mom


def guide(eos, start, end):
    a,b=read(start),read(end)
    assert a["plan"]==b["plan"] and sha(eos)==a["plan"]["eos_sha256"]
    assert all(t["hash"]==digest(t["snapshot"]) for t in (a,b))
    s=a["snapshot"]
    keys=("entropy_variable_scalar","entropy_variable_fermion")
    p=np.asarray(s["lab_momentum_edges_GeV"])
    T=s["T_visible_GeV"]
    q=pack([s[k] for k in keys])
    H=np.sqrt((TabulatedRadiation(eos,T).state(1.)[1]+get_mom(p,[s[k] for k in keys])[2].sum())/3.)/MPL
    op=LocalAbsolutePrompt(Collision(p),H,[],True)
    dt=(b["time"]-a["time"])/6.
    stats=[]
    for _ in range(3):
        guess=predictor(op,q,dt)
        q,stat=advance(op,q,dt,T,initial_theta=guess,vacuum_scaled=True,balanced_merit=True)
        stats.append(stat)
    mid=dict(s,**{k:v.tolist() for k,v in zip(keys,op.arrays(q),strict=True)})
    return dict(snapshot=mid,hash=digest(mid),plan=a["plan"],time=(a["time"]+b["time"])/2,
                input_hash=a["hash"],end_source_hash=b["hash"],microstep_audits=stats)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("eos","start","end","output"):
        p.add_argument(k,type=Path)
    a=p.parse_args()
    atomic(a.output,guide(a.eos,a.start,a.end))
