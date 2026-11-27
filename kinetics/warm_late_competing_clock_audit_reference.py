"""Separate the coupled local perturbation from the turning clock's time error."""

import argparse
import json
from pathlib import Path

import numpy as np

from warm_corrected_late_prompt_cut_aligned_reference import excess
from warm_light_occupation_stage_audit_reference import digest, read
from warm_prompt_reentry_shape_audit_reference import get_mom


def audit(raw, D):
    a=read(raw)
    A=a["snapshot"]
    assert a["hash"]==digest(A)
    C=a["plan"]["normalized_cut_GeV"]
    p=np.asarray(A["lab_momentum_edges_GeV"])
    rho=get_mom(p,[A["entropy_variable_"+k] for k in ("scalar","fermion")])[2].sum()
    initial=sum(excess(A,C))/rho
    values=[]
    rows=[]
    for mode in ("prompt","all"):
        c,m,f=[read(D/f"compare-{mode}-{tag}.json") for tag in ("coarse","mid","fine")]
        assert all(r["fixed_energy_cut_GeV"]==C and r["mode"]==mode for r in (c,m,f))
        assert all(r["hash"]==digest(r["snapshot"]) for r in (c,m,f))
        assert c["input_hash"]==m["input_hash"]==a["hash"] and f["input_hash"]==m["hash"]
        assert c["raw_target_hash"]==f["raw_target_hash"]
        assert c["relative_log_expansion"]==m["relative_log_expansion"]+f["relative_log_expansion"]
        end=[r["snapshot"] for r in (c,f)]
        assert all(np.allclose(end[0][k],end[1][k],rtol=4e-15,atol=0.) for k in (
            "lab_momentum_edges_GeV","T_visible_GeV","scale_factor"))
        W=[sum(excess(z,C))/rho for z in end]
        values.append(W)
        rows.append(dict(mode=mode,first_hash=a["hash"],coarse_mid_fine_hashes=[r["hash"] for r in (c,m,f)],
                         W_coarse_and_fine_over_initial_rho=W,relative_time_resolution_change=W[1]/W[0]-1.,
                         exact_local_step_residuals=[r["solver"]["residual"] for r in (c,m,f)]))
    values=np.array(values)
    delta=values[1]/values[0]-1.
    assert np.max(abs(delta))<2e-7
    return dict(initial_W_over_rho=initial,rows=rows,coupled_relative_delta_coarse_and_fine=delta.tolist())


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("raw",type=Path)
    p.add_argument("tapes",type=Path)
    a=p.parse_args()
    print(json.dumps(audit(a.raw,a.tapes),indent=2,allow_nan=False))
