"""Replay the same first three exact-vacuum microsteps with the local Jacobian."""

import argparse
import json
from pathlib import Path
import time

import numpy as np

from joint_prompt_two_moment_evolution_reference import Collision
from joint_resonant_centered_stiff_reference import advance
from joint_resonant_frw_split_reference import pack
from late_radiation_background_reference import MPL,TabulatedRadiation
from warm_compact_absolute_prompt_reference import LocalAbsolutePrompt
from warm_corrected_late_prompt_cut_aligned_reference import predictor
from warm_light_occupation_stage_audit_reference import digest,read
from warm_prompt_reentry_shape_audit_reference import get_mom


def run(eos,seed,endpoint):
    C,D=map(read,(seed,endpoint))
    S,X=C["snapshot"],D["snapshot"]
    assert C["plan"]==D["plan"] and C["time"]==0. and D["time"]==1e-11 and D["plan"]["substeps"]==3
    assert C["plan"]["artificial_empty_shell_filling"]==0.
    assert C["hash"]==digest(S) and D["hash"]==digest(X)
    p=np.asarray(S["lab_momentum_edges_GeV"])
    th=[np.asarray(S["entropy_variable_"+k]) for k in ("scalar","fermion")]
    other=[np.asarray(X["entropy_variable_"+k]) for k in ("scalar","fermion")]
    N,K,E=get_mom(p,th)
    H=np.sqrt((TabulatedRadiation(eos,S["T_visible_GeV"]).state(1.)[1]+E.sum())/3)/MPL
    model=LocalAbsolutePrompt(Collision(p),H,[],True)
    u=pack(th)
    t=time.perf_counter()
    trace=[]
    n=3
    for _ in range(n):
        h=D["time"]/n
        u,z=advance(model,u,h,S["T_visible_GeV"],initial_theta=predictor(model,u,h),
                    vacuum_scaled=True,balanced_merit=True)
        trace.append(z)
    r=get_mom(p,model.arrays(u))
    g=get_mom(p,other)
    # In addition to bulk closure, check the tiny activated shell ratios.
    bulk=[float(np.sum(abs(a-b))/np.sum(b)) for a,b in zip(r,g,strict=True)]
    local=[float(np.max(abs(a-b)/np.maximum(np.maximum(a,b),1e-200))) for a,b in zip(r,g,strict=True)]
    count=float(np.max(abs(r[0].sum(axis=1)/N.sum(axis=1)-1.)))
    energy=float(abs(r[2].sum()/E.sum()-1.))
    assert max(bulk)<2e-12 and max(local)<2e-8 and max(count,energy)<5e-13
    return dict(scope="Same exact-vacuum old moments and timestep. Local Jacobian changes no rate.",
                input_state_sha256=C["hash"],reference_end_sha256=D["hash"],
                maximum_shell_relative_errors_number_kinetic_total=local,normalized_global_errors=bulk,
                initial_moment_ledgers=[count,energy],solver=trace,
                local_control_wall_time_seconds=time.perf_counter()-t)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("eos","seed","endpoint"):
        p.add_argument(k,type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.seed,a.endpoint),indent=2,allow_nan=False))
