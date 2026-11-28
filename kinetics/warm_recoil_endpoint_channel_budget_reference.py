"""Signed native first-cell reaction budget at the post-recoil checkpoints."""

import argparse
import json
from pathlib import Path

import numpy as np

from frw_exterior_restart_reference import import_restart
from joint_below_pole_elastic_linear_modes_reference import setup
from joint_prompt_two_moment_evolution_reference import Collision
from joint_resonant_frw_split_reference import pack
from late_archived_shape_continuation_reference import Late
from late_radiation_background_reference import MPL,TabulatedRadiation
from scalar_pair_expanding_continuation_reference import rate
from warm_light_occupation_stage_audit_reference import digest,read
from warm_prompt_reentry_shape_audit_reference import get_mom


def at(eos,S,settings):
    p=np.asarray(S["lab_momentum_edges_GeV"])
    th=[np.asarray(S["entropy_variable_"+b]) for b in ("scalar","fermion")]
    bath=TabulatedRadiation(eos,S["T_visible_GeV"])
    grav=get_mom(p,th)[2].sum()+import_restart(S["positive_exterior_restart"]).energy(1.)
    T,rho=bath.state(1.)
    H=np.sqrt((rho+grav)/3)/MPL
    c=Collision(p); base=Late(c,H,[],True)
    u=pack(th)
    ref,_,td=base.reference(u,T)
    prompt=c.rate(th)[1,0,0]/(H*c.c[1].vol[0])
    out=[]
    for e,pair in settings:
        allp=[]
        setup(max(td,T),c,td,base.eta,*e,32.,.82,cache=allp,cache_only=True)
        model=Late(c,H,allp,True,compact=True)
        q,_,_=model.reference(u,T)
        elastic=model.elastic(u-q)[0].reshape(2,c.c[1].p.size-1,2)[1,0,0]/c.c[1].vol[0]
        B,err,_=rate(p,th,pair)
        assert err<2e-11
        f=B[1,0,1]/(H*c.c[1].vol[0])
        out.append(dict(elastic_order=e,pair_order=pair,first_cold_shell_rates_per_log_a=[prompt,elastic,f],
                        sum_per_log_a=prompt+elastic+f))
    return out


def run(D,eos):
    a=read(D/"warm-recoil-crossing-matched-2026-11-26.json.gz")
    b=a["frames"][-1]
    rows=[]
    for q,z in zip((32,64),b["controls"][:2],strict=True):
        A=read(D/f"warm-recoil-crossing-q{q}-terminal-2026-11-26.json.gz")[-1]
        S=A["snapshot"]
        assert digest(S)==A["snapshot_sha256"]==b["checkpoints"][int(q==64)]
        x=at(eos,S,[((3,2,2,8),(3,7,6,2)),((6,4,4,12),(4,10,10,4))])
        assert all(y["sum_per_log_a"]<0 for y in x)
        rows.append(dict(q=q,T_GeV=S["T_visible_GeV"],controls=x,
                         accepted_full_step_rate_per_log_a=z["mean_occupation_change_per_log_a_over_preceding_interval"]))
    return dict(scope="Frozen native operator budget for the first retained shell at 8.133 MeV. "
                      "Exact component support and observed finite-step mean listed separately; "
                      "no assignment of the measured finite-step difference to untested channels.",
                labels=["prompt","selected_below_pole_elastic","signed_scalar_pair"],cases=rows)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("eos",type=Path)
    p.add_argument("--tapes",type=Path,default=Path("results"));a=p.parse_args()
    print(json.dumps(run(a.tapes,a.eos),indent=2,allow_nan=False))
