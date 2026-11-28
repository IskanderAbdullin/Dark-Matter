"""Compose any number of contiguous warm count archives."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from three_species_decay_onsager_reference import m,mp
from warm_coupled_count_ledger_reference import endpoint,identities,read,row


def run(budget,start_shape,stages):
    B=read(budget)
    b=(B["cases"] if isinstance(B,dict) else B)[-1]
    assert b["enabled"] and b["signed_visible_portal"]
    s = read(start_shape)["snapshots"][-1]
    a0=s["scale_factor"]
    mass=np.array([mp,m])
    n=np.asarray(s["represented_species_count_a3"])
    norm=np.asarray(b["retained_fixed_Omega_h2"])/(mass*n)
    assert abs(norm[0]/norm[1]-1.)<3e-12
    ex=np.asarray(b["exterior_fixed_Omega_h2"])/(mass*norm)
    before=float((n+ex).sum())
    result=[]
    for f,g in stages:
        r=read(f)["rows"][-1]
        v=read(g)["snapshots"][-1]
        last=float(n[1]+ex[1])
        n,ex=endpoint(r,n,ex)
        identities(r,v,s,n)
        assert r["steps"]==v["steps"] and r["with_pair"]
        s=v
        res=row(r,n,ex,norm,mass)
        result.append(dict(relative_expansion_from_12_MeV=s["scale_factor"]/a0,
             fractional_inclusive_stable_loss_from_previous=1.-(n[1]+ex[1])/last,
             inclusive_number_ledger_from_start=float((n+ex).sum()/before-1.),
             predecessor_occupation_audit_only=True,
             two_file_sha256=[hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in (f,g)],
             **res))
    return result


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("budget")
    p.add_argument("start_shape")
    p.add_argument("--stage",action="append",nargs=2,required=True,metavar=("ROWS","SHAPE"))
    a=p.parse_args()
    print(json.dumps(dict(stages=run(a.budget,a.start_shape,a.stage),
        scope="Prior exterior counts included exactly once; prior exterior energy and absorbing/stimulating tail occupancies are not evolved by resumed finite-support operators."),
        indent=2,allow_nan=False))
