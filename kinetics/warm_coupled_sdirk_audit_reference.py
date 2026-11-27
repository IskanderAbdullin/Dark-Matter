"""Validate accepted phases and compare the collision-clock refinements."""
import argparse
import json
from pathlib import Path

import numpy as np

from warm_light_occupation_stage_audit_reference import digest,read


def run(path,be1,be2):
    choices=[be1,be2,path/"068-abs1",path/"068-abs2"]
    rows=[]
    ref=None
    for name,step in zip(choices,(1,2,1,2),strict=True):
        a=read(name/"summary.json")
        r=a["rows"]
        assert len(r)==step and a["mode"]=="expand-finite-all"
        q=(a["input_sha256"],a["cut_GeV"],r[-1]["c"],r[-1]["qsha256"],r[-1]["wsha256"],a["W_initial_over_rho"])
        if ref is None: ref=q
        assert q==ref
        peak_stage_defect=0.
        for j,b in enumerate(r,1):
            c=read(name/f"step-{j}.json")
            assert c["row"]==b and digest(c["daughter_snapshot"])==b["daughter_hash"]
            assert b["solver"]["residual"]<8e-12
            assert max(abs(b["relative_energy_defect"]),*map(abs,b["relative_N_defects"]))<3e-10
            if name.name.startswith("068-abs"):
                assert b["collision_clock_method"]=="sdirk2"
                ss=b["solver"]
                assert len(ss["max_relative_invariant_defect_at_each_stage"])==2
                peak_stage_defect=max(peak_stage_defect,*ss["max_relative_invariant_defect_at_each_stage"])
                assert 0.<=ss["first_stage_parent_min_max"][0]<=ss["first_stage_parent_min_max"][1]<=1.
        rows.append(dict(name=name.name,W=r[-1]["W_over_initial_daughter_rho"],rows=r,
                         max_collision_stage_invariant_defect=peak_stage_defect))
    assert all(r["W"]<ref[-1] for r in rows)
    ww=np.array([r["W"] for r in rows]).reshape(2,2)
    return dict(scope="Same fixed chart and local FRW/pair/elastic split; only the implicit collision clock differs.",
                rows=rows,one_to_two_relative_change_by_method=(ww[:,1]/ww[:,0]-1.).tolist(),
                W_differences_absolute_by_method=abs(ww[:,1]-ww[:,0]).tolist())


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("sdirk","be1","be2"):
        p.add_argument(name,type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.sdirk,a.be1,a.be2),indent=2,allow_nan=False))
