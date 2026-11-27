"""Separate the local finite-parent correction from temporal resolution."""
import argparse
import json
from pathlib import Path

import numpy as np

from warm_light_occupation_stage_audit_reference import read


def audit(first, refined, mode0, mode1):
    a,b=map(read,(first,refined))
    for key in ("id","initial","comoving_chart"):
        assert a[key]==b[key]
    x0,x1=[a["endpoint"][v] for v in (mode0,mode1)]
    y0,y1=[b["endpoint"][v] for v in (mode0,mode1)]
    assert x0["accepted_blocks"]==x1["accepted_blocks"]==1
    assert y0["accepted_blocks"]==y1["accepted_blocks"]==2
    assert len({t["expansion_factor"] for t in (x0,x1,y0,y1)})==1
    w=np.array([[x0["W_over_initial_rho"],y0["W_over_initial_rho"]],
                [x1["W_over_initial_rho"],y1["W_over_initial_rho"]]])
    assert np.all(w>0.)
    return dict(input_hash=a["id"],initial=a["initial"],comoving_chart=a["comoving_chart"],
                modes=[mode0,mode1],W_stationary_and_finite_by_clock=w.tolist(),
                finite_parent_shift_at_fixed_clock=(w[1]/w[0]-1.).tolist(),
                time_refinement_shift_at_fixed_physics=(w[:,1]/w[:,0]-1.).tolist(),
                decline_from_anchor=(w<a["initial"]).tolist())


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("coarse",type=Path)
    p.add_argument("fine",type=Path)
    p.add_argument("--all",action="store_true")
    a=p.parse_args()
    ext="-all" if a.all else ""
    print(json.dumps(audit(a.coarse,a.fine,"expand-stationary"+ext,"expand-finite"+ext),indent=2,allow_nan=False))
