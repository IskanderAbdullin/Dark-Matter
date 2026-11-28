"""Separate the two daughter contributions to the signed prompt cut current."""

import argparse
import json
from pathlib import Path

import numpy as np

from late_radiation_background_reference import MPL,TabulatedRadiation
from three_species_decay_onsager_reference import M
from warm_light_occupation_stage_audit_reference import digest,read
from warm_prompt_excess_energy_tail_reference import integrate
from warm_prompt_reentry_shape_audit_reference import get_mom


def run(eos,prior,paths):
    C=np.hypot(M,1./read(prior)["branches"][-1]["entropy_expansion_from_original_anchor"])
    out=[]
    for p in paths:
        A=read(p)
        S=A["snapshot"]
        assert A["hash"]==digest(S)
        q=np.asarray(S["lab_momentum_edges_GeV"])
        v=[S["entropy_variable_"+k] for k in ("scalar","fermion")]
        rho=get_mom(q,v)[2].sum()
        H=np.sqrt((TabulatedRadiation(eos,S["T_visible_GeV"]).state(1.)[1]+rho)/3)/MPL
        r=[]
        for order in (6,10,14):
            diag=np.empty((2,2))
            a=integrate(S,C,order,True,split_result=diag)/(rho*H)
            diag/=rho*H
            signed=diag[:,1]-diag[:,0]
            assert abs(signed.sum()-a[3])<=6e-14*max(diag.sum(),1e-245)
            r.append(dict(order=order,per_species_consumed_emitted=diag.tolist(),
                          per_species_signed_change=signed.tolist(),total_signed=float(a[3])))
        out.append(dict(snapshot_sha256=A["hash"],time=A["time"],hard_GeV=A["plan"]["hard_GeV"],rate_controls=r))
    return dict(scope="Closed stationary prompt diagnostic on already evolved phases. "
                      "The species ledgers are convex excess, not net particle production.",
                energy_cut_GeV=C,species=["scalar","light"],cases=out)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("eos",type=Path);p.add_argument("prior",type=Path);p.add_argument("paths",type=Path,nargs="+")
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.prior,a.paths),indent=2,allow_nan=False))
