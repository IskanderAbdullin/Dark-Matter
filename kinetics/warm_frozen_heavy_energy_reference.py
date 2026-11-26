"""Integrate the decay-only stationary heavy gas on archived retained daughters."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre

from warm_archived_bose_pair_weighting_reference import snapshot
from joint_visible_heavy_adiabatic_reference import snapshot_functions,parent_branch
from warm_isolated_recycling_energy_envelope_reference import moments
from three_species_decay_onsager_reference import M


def run(paths):
    out=[]
    for p in paths:
        S=snapshot(p)
        f,ef,g,eg=snapshot_functions(S)
        upper=ef[-1]+eg[-1]
        cut=np.sqrt(max(upper**2-M**2,0.))
        values=[]
        for n in (4,8,14):
            edges=np.linspace(0,cut,70)
            x,w=roots_legendre(n)
            a,b=(edges[:-1,None]+edges[1:,None])/2,np.diff(edges)[:,None]/2
            if cut:
                P=(a+b*x).ravel()
                E=np.hypot(P,M)
                F=parent_branch(E,f,g,ef,eg,angle=15)[:,0]
                R=2/np.pi**2*np.sum((b*w).ravel()*P**2*E*F)
            else:
                R=0.
            values.append(float(R/moments(S).sum(axis=0)[1]))
        assert abs(values[-1]-values[-2]) < 1e-6*max(values[-1],1e-20)
        out.append(dict(name=Path(p).name,T=S["T_visible_GeV"],fractions=values,
                        maximum_heavy_P=float(cut)))
    return dict(scope="Decay/inverse-decay stationary heavy state on frozen retained "
                      "daughter shapes, vacuum outside their mesh, no additional sources. "
                      "Not a time-dependent upper on the heavy energy fraction.",
                radial_panel_orders=[4,8,14],cases=out)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("paths",nargs="+")
    a=p.parse_args()
    print(json.dumps(run(a.paths),indent=2,allow_nan=False))
