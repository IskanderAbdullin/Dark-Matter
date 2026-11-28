"""Bound the ballistic mean path of the direct spectra without within-shell placement."""
import argparse
import json
from pathlib import Path

import numpy as np

from archived_collisionless_streaming_reference import Background
from three_species_decay_onsager_reference import mp
from warm_light_occupation_stage_audit_reference import read


def run(eos,spec):
    bg=Background(eos)
    r=read(spec)
    q=np.array([x if x is not None else np.inf for x in r["q_over_T0_edges"]])
    assert np.isinf(q[-1]) and q[0]==0 and np.all(np.diff(q)>0)
    # Start every spectrum at the same late point, after its represented
    # production. This is an apples-to-apples ballistic diagnostic.
    ai=bg.scale(.006)
    x=q[1:-1]*bg.T0
    d=bg.streams(x,mp,ai,bg.aeq,7)
    cf=bg.streams(x,mp,ai,bg.aeq,12)
    assert np.max(abs(d/cf-1))<1e-8
    finite=np.r_[0.,cf]
    # H >= H0 sqrt(omega_m) a^(-3/2) and v <= c, including for
    # the handful of unbounded q beyond the last finite edge.
    massless_upper=2*bg.D0*(np.sqrt(bg.aeq)-np.sqrt(ai))/np.sqrt(bg.om)
    assert finite[-1]<massless_upper
    rows=[]
    for item in r["rows"]:
        weights=np.asarray(item["shell_lepton_then_pion_yields"])
        n=weights.sum(axis=1)
        assert np.all(weights>=0.) and np.all(n>0)
        low=weights@finite
        high=weights@np.r_[finite[1:],massless_upper]
        assert np.all(low<high)
        mix=np.array([low.sum(),high.sum()])/n.sum()
        per=np.stack((low/n,high/n),axis=-1)
        rows.append(dict(T_reheat_GeV=item["T_reheat_GeV"],mean_path_Mpc_mixture_lower_upper=mix.tolist(),
                         mean_path_Mpc_per_species_lower_upper=per.tolist(),
                         pion_mixture_fraction=n[1]/n.sum()))
    return dict(scope="Number-weighted ballistic path after the common 6 MeV end of the represented "
                      "cutoff, through equality, in one prescribed background. Monotone momentum "
                      "endpoints, exact for the shell bounds up to positive quadrature error. "
                      "Not a transfer function or an observation-derived constraint.",
                scale_factor_start=float(ai),scale_factor_equality=float(bg.aeq),massless_tail_upper=massless_upper,rows=rows)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("eos",type=Path)
    p.add_argument("spec",type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.spec),indent=2,allow_nan=False))
