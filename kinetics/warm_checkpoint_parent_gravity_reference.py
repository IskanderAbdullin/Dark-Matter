"""Measure the active-minus-passive parent contribution to the local Hubble rate."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from late_radiation_background_reference import TabulatedRadiation
from three_species_decay_onsager_reference import M
from warm_comoving_durable_restart_reference import restore
from warm_comoving_parent_chart_reference import ComovingCollision,common_rule
from warm_coupled_serial_restart_reference import sha
from warm_finite_parent_comoving_expansion_reference import parent_equilibrium
from warm_light_occupation_stage_audit_reference import digest,read
from warm_prompt_reentry_shape_audit_reference import get_mom


def h(x):
    return hashlib.sha256(np.asarray(x,dtype="<f8").tobytes()).hexdigest()


def run(eos,start,last,d,order,sub):
    a,b=read(start),read(last)
    S=a["snapshot"]
    assert a["hash"]==digest(S) and b["hash"]==digest(b["snapshot"])
    assert a["plan"]["eos_sha256"]==sha(eos) and a["plan"]==b["plan"]
    p=np.asarray(S["lab_momentum_edges_GeV"])
    scale=np.exp(np.linspace(0.,b["time"]-a["time"],2*sub+1))
    Q,w=common_rule(p,scale,order)
    path=d/"restart-state.json"
    ps=read(path)
    k=ps["payload"]["accepted_blocks"]
    rows=[read(d/f"step-{i}.json")["row"] for i in range(1,k+1)]
    assert all(r["qsha256"]==h(Q) and r["wsha256"]==h(w) for r in rows)
    archive,_,T,ex,active=restore(path,ps["payload"]["key"],
                    np.r_[1.,[r["c"] for r in rows]],a["hash"],Q,h)
    # If a writer advances during this audit, require the caller to retry on
    # the new exact accepted image; do not silently compare mixed epochs.
    assert archive["accepted_blocks"]==k
    keys=("entropy_variable_scalar","entropy_variable_fermion")
    theta=[np.asarray(S[k]) for k in keys]
    theta_out=[np.asarray(T[k]) for k in keys]
    kernel=ComovingCollision(p,1.,Q,w)
    passive=parent_equilibrium(kernel,theta)
    c=ex.a
    energy_per_q=w*np.hypot(M,Q/c)/c**3
    pa=float(energy_per_q@passive)
    ac=float(energy_per_q@active)
    delta=float(energy_per_q@(active-passive))
    daughter=get_mom(np.asarray(T["lab_momentum_edges_GeV"]),theta_out)[2].sum()
    radiation=TabulatedRadiation(eos,S["T_visible_GeV"]).state(c)[1]
    rest=daughter+radiation+ex.energy(c)/c**4
    x=delta/(rest+pa)
    frac=x/(1.+np.sqrt(1.+x))
    assert abs((ac-pa)-delta)<2e-14*max(ac,pa,1e-240)
    return dict(accepted_steps=k,absolute_H_clock_time=rows[-1]["time_J"],
                relative_scale=c,input_sha256=a["hash"],active_parent_rho=ac,passive_parent_rho=pa,
                physical_delta_parent_rho=delta,active_minus_passive_H_over_passive_H=float(frac),
                daughter_rho=float(daughter),radiation_rho=float(radiation),checkpoint_sha256=ps["sha256"],
                scope="Change from replacing passive parent density in the same instantaneous Friedmann denominator; not a rerun with active gravity or an integrated error bound.")


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("eos","start","last","folder"):
        p.add_argument(name,type=Path)
    p.add_argument("--order",type=int,default=2)
    p.add_argument("--common-subdivisions",type=int,default=2)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.start,a.last,a.folder,a.order,a.common_subdivisions),indent=2,allow_nan=False))
