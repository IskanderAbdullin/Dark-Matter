"""Check split count slopes against pair rates along a collision-only trajectory."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def run(eos,seed,leading,standard,small,smaller):
    d=json.loads(leading.read_text())
    A,B=d
    assert A["nodes"]==22 and B["nodes"]==43
    ref=np.asarray(B["composite_simpson_count_delta"])
    check=max(np.max(abs(np.asarray(z[k])/ref-1.)) for z in d
              for k in ("composite_trapezoid_count_delta","composite_simpson_count_delta"))
    assert check<2e-5 and np.max([z["collision_only_count_and_energy_error"] for z in d])<1e-9
    s=ref[0]/ref[1]
    assert ref[0]>0. and ref[1]<0.
    paths=(standard,small,smaller)
    out=[]
    for eps,p in zip((1.,.1,.03),paths,strict=True):
        q=json.loads(p.read_text())
        for name,items in q.items():
            assert items[0]["outer_subintervals"]==1 and items[-1]["outer_subintervals"]==4
            vals=[]
            for r in (items[0],items[-1]):
                v=(np.asarray(r["count_relative_to_initial"])-1.)/eps
                assert abs(v[0]/v[1]/s-1.)<1e-7
                err=v/ref-1.
                vals.append(dict(substeps=r["outer_subintervals"],normalized_difference_of_counts=v.tolist(),
                                 relative_to_collision_only_derivative_error=err.tolist()))
            out.append(dict(pair_strength_scale=eps,schedule=name,values=vals))
    for n in (0,1):
        x,y=out[-2:]
        assert abs(y["values"][n]["relative_to_collision_only_derivative_error"][0])<abs(
               x["values"][n]["relative_to_collision_only_derivative_error"][0])/10.
    return dict(leading_count_delta=ref.tolist(),quadrature_relative_spread=float(check),
                initial_and_final_net_rate_per_species=[B["net_rate_per_initial_species_per_log_b"][j] for j in (0,-1)],
                initial_and_final_inverse_over_forward=[B["inverse_over_forward"][j] for j in (0,-1)],
                comparisons=out,
                input_sha256={p.name:digest(p) for p in (eos,seed,leading,*paths)},
                scope="Derivative at zero artificial pair strength on a fixed warm collision geometry. The collision trajectory alone carries no pair feedback; finite-strength comparisons use actual split kicks. This does not pick a production solver or a full expanding solution.")


if __name__=="__main__":
    p=argparse.ArgumentParser()
    for key in ("eos","seed","leading","standard","small","smaller"):
        p.add_argument(key,type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.seed,a.leading,a.standard,a.small,a.smaller),indent=2,allow_nan=False))
