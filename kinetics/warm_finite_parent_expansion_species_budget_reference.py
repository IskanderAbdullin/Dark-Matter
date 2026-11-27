"""Resolve number, full energy and threshold-excess changes by daughter species."""
import argparse
import json
from pathlib import Path

import numpy as np

from warm_corrected_late_prompt_cut_aligned_reference import excess
from warm_light_occupation_stage_audit_reference import digest,read
from warm_prompt_reentry_shape_audit_reference import get_mom


def measure(s,C):
    q=np.array(s["lab_momentum_edges_GeV"])
    n,k,e=get_mom(q,[s["entropy_variable_"+j] for j in ("scalar","fermion")])
    return n.sum(axis=1),k.sum(axis=1),e.sum(axis=1),np.array(excess(s,C))


def run(anchor,*cases):
    r=read(anchor)
    s=r["snapshot"]
    assert r["hash"]==digest(s)
    C=r["plan"]["normalized_cut_GeV"]
    N0,K0,E0,W0=measure(s,C)
    rho=E0.sum()
    rows=[]
    for p in cases:
        t=read(p)
        x=t["daughter_snapshot"]
        z=t["row"]
        assert digest(x)==z["daughter_hash"] and z["initial_daughter_hash"]==r["hash"]
        a=z["c"]
        N,K,E,W=measure(x,C)
        assert abs(W.sum()/rho/z["W_over_initial_daughter_rho"]-1.)<3e-14
        rows.append(dict(case=p.parent.name,daughter_sha256=z["daughter_hash"],
            count_change_comoving_relative_by_species=(a**3*N/N0-1.).tolist(),
            daughter_a4_kinetic_energy_change_over_initial_rho=((a**4*K-K0)/rho).tolist(),
            daughter_a4_full_energy_change_over_initial_rho=((a**4*E-E0)/rho).tolist(),
            physical_excess_change_over_initial_rho=((W-W0)/rho).tolist()))
    return dict(scope="Observed per-species net changes; projected FRW trace and reaction exchange both contribute.",
                initial_hash=r["hash"],hard_GeV=r["plan"]["hard_GeV"],rows=rows,
                columns=["scalar","light_fermion"])


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("anchor",type=Path)
    p.add_argument("cases",type=Path,nargs="+")
    a=p.parse_args()
    print(json.dumps(run(a.anchor,*a.cases),indent=2,allow_nan=False))
