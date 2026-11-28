"""Quantify fixed-kernel scalar Bose multiplier changes under cold regridding."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_endpoint_scalar_spectral_moment_budget import load_card
from early_pole_endpoint_scalar_shape_transport import CUT


THRESHOLD=(.01,.05,.1,.20)
POINTS=(2.6e-6,3.5e-6,4.0e-6,4.7e-6,8.2e-6)


def pair(a,b):
    (d,u,_),(e,v,_)=a,b
    if d.ds.ed[0]!=e.ds.ed[0] or d.ds.ed[-1]!=e.ds.ed[-1]:
        raise ArithmeticError('Need one fixed scalar support')
    lo=d.ds.ed[0]
    z=np.union1d(np.union1d(d.ds.ed,e.ds.ed),np.array([CUT]))
    edge=z[(z>=lo)&(z<=CUT)]
    l,h=edge[:-1],edge[1:]
    q=(l+h)/2
    def f(state,grid,q):
        return state[np.searchsorted(grid,q,side='right')-1]
    a,b=f(u,d.ds.ed,q),f(v,e.ds.ed,q)
    V=(h-l)*(h*h+h*l+l*l)/(6*np.pi**2)
    X,Y=V*a,V*b
    r=(1+b)/(1+a)-1
    weights=[fsum(X),fsum(Y)]
    o=[]
    for k in THRESHOLD:
        region=abs(r)>=k
        o.append(dict(absolute_relative_bose_factor_threshold=k,
                      fraction_of_old_and_new_cold_counts=[fsum(w[region])/t for w,t in zip((X,Y),weights,strict=True)]))
    fixed=[float(z) for z in ((1+f(v,e.ds.ed,np.array(POINTS)))/(1+f(u,d.ds.ed,np.array(POINTS)))-1)]
    # Reconstruct the known integrals of f, without inferring a rate.
    dot=fsum(abs(X-Y))/weights[0]
    return dict(absolute_unscaled_shape_change_over_old_cold_count=dot,
                fractions_above_change_thresholds=o,
                pointwise_Bose_factor_relative_changes=fixed)


def run(paths):
    levels=(2,4,8,16)
    data=[load_card(p,n) for p,n in zip(paths,levels,strict=True)]
    return dict(cold_upper=CUT,point_q=POINTS,
                cases=[dict(scalar_levels=[a,b],**pair(x,y)) for a,b,x,y in
                    zip(levels[:-1],levels[1:],data[:-1],data[1:],strict=True)],
                scope='Frozen kernel would multiply each forward event by (1+f). Count-weighted support, not rate-weighted or an integrated rate.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('state',type=Path,nargs=4)
    a=p.parse_args()
    print(json.dumps(run(a.state),indent=2,allow_nan=False))
