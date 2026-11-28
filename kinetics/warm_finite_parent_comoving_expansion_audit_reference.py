"""Verify matched-clock tails and the decoded daughter endpoints."""

import argparse
import json
from pathlib import Path

import numpy as np

from warm_corrected_late_prompt_cut_aligned_reference import excess
from warm_light_occupation_stage_audit_reference import digest,read
from warm_prompt_reentry_shape_audit_reference import get_mom


def audit(raw, root, tag, directions):
    A=read(raw)
    S=A["snapshot"]
    assert digest(S)==A["hash"]
    q=np.asarray(S["lab_momentum_edges_GeV"])
    _,_,e=get_mom(q,[S["entropy_variable_"+k] for k in ("scalar","fermion")])
    rho=e.sum()
    C=A["plan"]["normalized_cut_GeV"]
    init=sum(excess(S,C))/rho
    runs={}
    for direction in directions:
        name=tag+direction
        data=read(root/name/"summary.json")
        assert data["mode"]==direction and data["input_sha256"]==A["hash"]
        assert abs(data["W_initial_over_rho"]/init-1.)<3e-14
        rows=[]
        for j,part in enumerate(data["rows"],1):
            raw=read(root/name/f"step-{j}.json")
            assert raw["row"]==part
            V=raw["daughter_snapshot"]
            assert part["daughter_hash"]==digest(V)
            assert np.allclose(q/part["c"],V["lab_momentum_edges_GeV"],rtol=2e-15)
            assert abs(sum(excess(V,C))/rho/part["W_over_initial_daughter_rho"]-1.)<3e-14
            assert part["solver"]["residual"]<8e-12
            assert max(abs(part["relative_energy_defect"]),max(abs(x) for x in part["relative_N_defects"]))<3e-10
            rows.append(part)
        runs[direction]=rows
    sha={(r["qsha256"],r["wsha256"]) for rows in runs.values() for r in rows}
    assert len(sha)==1
    rel={}
    for key,a in runs.items():
        b=a[-1]
        assert all(r["time_J"]==a[j]["time_J"] for j,r in enumerate(runs[directions[0]]))
        rel[key]=dict(W_over_initial_rho=b["W_over_initial_daughter_rho"],
                      W_change_from_anchor=b["W_over_initial_daughter_rho"]-init,
                      expansion_factor=b["c"],accepted_blocks=len(a))
    for phase,tail in (("frozen",""),("expand",""),("expand","-all")):
        fs,ff=(rel.get(phase+"-"+v+tail) for v in ("stationary","finite"))
        if fs and ff:
            ff["finite_relative_shift"]=ff["W_over_initial_rho"]/fs["W_over_initial_rho"]-1.
            assert abs(ff["finite_relative_shift"])<5e-5
    return dict(id=A["hash"],initial=init,endpoint=rel,comoving_chart=list(sha)[0])


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("raw",type=Path)
    p.add_argument("root",type=Path)
    p.add_argument("tag")
    p.add_argument("directions",nargs="+")
    a=p.parse_args()
    print(json.dumps(audit(a.raw,a.root,a.tag,a.directions),indent=2,allow_nan=False))
