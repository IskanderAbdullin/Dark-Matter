"""Report actual above-cut daughter energy separately from convex excess."""

import argparse
import json
from pathlib import Path

import numpy as np

from three_species_decay_onsager_reference import M,m,mp
from warm_light_occupation_stage_audit_reference import digest,read
from warm_prompt_reentry_shape_audit_reference import get_mom


def run(prior,*paths):
    C=np.hypot(M,1./read(prior)["branches"][-1]["entropy_expansion_from_original_anchor"])
    out=[]
    for D in paths:
        rows=[]
        for p in sorted(D.glob("step-*.json"),key=lambda f:int(f.stem.split("-")[-1])):
            a=read(p)
            S=a["snapshot"]
            assert a["hash"]==digest(S)
            pp=np.asarray(S["lab_momentum_edges_GeV"])
            q=[np.asarray(S["entropy_variable_"+k]) for k in ("scalar","fermion")]
            N,K,E=get_mom(pp,q)
            n,e=[],[]
            for i,d in enumerate((mp,m)):
                edges=np.hypot(d,pp)
                assert np.min(abs(edges-C))<4e-16 and not np.any((edges[:-1]<C-4e-16)&(edges[1:]>C+4e-16))
                take=edges[1:]>C+4e-16
                n.append(N[i,take].sum())
                e.append(E[i,take].sum())
            e=np.asarray(e);n=np.asarray(n)
            rho=E.sum()
            excess=(e-C*n)/rho
            residue=np.sum(abs(excess-a["excess_over_initial_energy"]))
            assert residue<1e-10*max(excess.sum(),1e-60)
            rows.append(dict(time_H=a["time"],snapshot_sha256=a["hash"],shell_aligned_excess_over_retained_rho=excess.tolist(),
                             total_energy_above_cut_over_retained_rho=(e/rho).tolist(),
                             count_above_cut_over_same_species_total=(n/N.sum(axis=1)).tolist(),
                             excess_closure_difference=float(residue)))
        out.append(dict(hard_GeV=a["plan"]["hard_GeV"],reported_steps=len(rows),rows=rows))
    return dict(scope="Both metrics are state moments for the aligned finite, closed prompt-only calculation. "
                      "No inference drawn for unrepresented particles or a cosmological resource.",
                cut_GeV=C,cases=out)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("prior",type=Path)
    p.add_argument("directories",type=Path,nargs="+")
    a=p.parse_args()
    print(json.dumps(run(a.prior,*a.directories),indent=2,allow_nan=False))
