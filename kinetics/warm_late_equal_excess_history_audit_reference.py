"""Verify the live descendants of the two archived equal-W phases."""
import argparse
import json
from pathlib import Path

import numpy as np

from warm_corrected_late_prompt_cut_aligned_reference import excess
from warm_light_occupation_stage_audit_reference import digest, read
from warm_prompt_reentry_shape_audit_reference import get_mom


def run(equal, dirs):
    z=read(equal)
    paths=[]
    times=None
    for x,p in zip(z["cases"],dirs,strict=True):
        r=read(p/"summary.json")
        assert r["plan"]["seed_sha256"]==x["snapshot_sha256"]
        assert r["plan"]["normalized_cut_GeV"]==z["cut_GeV"]
        s0=read(p/"step-0.json")["snapshot"]
        q=np.asarray(s0["lab_momentum_edges_GeV"])
        n,_,e=get_mom(q,[s0["entropy_variable_"+k] for k in ("scalar","fermion")])
        rho=float(e.sum())
        np.testing.assert_allclose(n.sum(axis=1),x["species_N"],rtol=6e-15)
        np.testing.assert_allclose(e.sum(axis=1),x["species_E"],rtol=6e-15)
        W=[]
        J=[]
        here=[]
        for j,a in enumerate(r["rows"]):
            b=read(p/f"step-{j}.json")
            s=b.pop("snapshot")
            assert a==b and digest(s)==a["hash"]
            assert max(a["discrete_N_and_E_errors"])<8e-12
            assert abs(sum(excess(s,z["cut_GeV"]))/rho/a["tail_energy"]-1)<5e-11
            W.append(a["tail_energy"])
            J.append(a["signed_native_flux"])
            here.append(a["time"])
        assert abs(W[0]/z["W_over_rho_target"]-1)<3e-11
        assert np.sign(J[0])==np.sign(x["angular"][-1]["signed"])
        if times is None:
            times=here
        else:
            assert times==here
        paths.append(dict(initial_sha256=x["snapshot_sha256"],Hdt=times,
                          W_over_rho=W,current=J,W_terminal_over_initial=W[-1]/W[0],
                          sampled_brackets=[times[i:i+2] for i in range(len(J)-1) if J[i]*J[i+1]<0],
                          max_relative_rise=float(max(W)/W[0]-1.)))
    return dict(cases=paths,terminal_W_ratio=paths[0]["W_over_rho"][-1]/paths[1]["W_over_rho"][-1],
                scope="Only sampled fixed-background closed prompt trajectories, not fitted reversal times.")


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("equal",type=Path)
    p.add_argument("folders",nargs=2,type=Path)
    args=p.parse_args()
    print(json.dumps(run(args.equal,args.folders),indent=2,allow_nan=False))
