"""Locate the two non-equivalent local FD mode-clock crossings."""
import argparse
import hashlib
import json
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from fd_portal_differential_reciprocity_reference import PositivePortal
from fd_portal_expansion_forced_modes import parcel
from sm_dark_entropy_exchange_reference import StandardModel


def run(path):
    eos=StandardModel(path)
    p,q=[PositivePortal(o) for o in ((13,36,26),(17,44,32))]
    low,high=[lru_cache(None)(lambda T,j=j,n=n:parcel(eos,j,n,T)['operators'])
              for j,n in ((p,112),(q,168))]
    rows=[]
    for i,tag in enumerate(('timelike','combined')):
        for k,label,interval in ((0,'slow',(.55,.64135)),(1,'fast',(.30,.45))):
            ans=[]
            for O in (low,high):
                def objective(t):
                    return O(float(t))[i]['positive_mode_rates_over_H'][k]-1.
                ts=np.linspace(*interval,17)
                pre=np.array([objective(t) for t in ts])
                if pre[0]>=0. or pre[-1]<=0. or np.any(np.diff(pre)<=0.):
                    raise ArithmeticError('Expected one sign reversal in the sampled bracket')
                j=np.flatnonzero(pre[:-1]*pre[1:]<0)
                if len(j)!=1:
                    raise ArithmeticError('Clock crossing cannot be identified by the declared samples')
                a,b=ts[j[0]:j[0]+2]
                t=brentq(objective,a,b,xtol=4e-13)
                o=O(t)[i]
                if abs(objective(t))>2e-10:
                    raise ArithmeticError('Numerical mode clock root failed to close')
                surrogate=np.array([objective(t+h)-objective(t-h) for h in (5e-6,1e-5)])/np.array([1e-5,2e-5])
                if min(surrogate)<=0. or abs(surrogate[0]/surrogate[1]-1.)>2e-4:
                    raise ArithmeticError('Independent local slope was not resolved')
                ans.append(dict(crossing_T=t,sampled_bracket=[a,b],rate_crossing_slope=surrogate[0],
                                both_rates_over_H=o['positive_mode_rates_over_H'],
                                algebraic_eta_and_beta=o['algebraic_eta_and_beta']))
            d=abs(ans[0]['crossing_T']-ans[1]['crossing_T'])
            if d>3e-6:
                raise ArithmeticError('Two independent collision resolutions disagree on the local crossing')
            rows.append(dict(operator=tag,mode=label,maximum_cross_mesh_shift=d,**ans[1]))
    return dict(rows=rows,table_input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                scope='Clock equality for local linearized vacuum FD modes, not a decoupling temperature or trajectory determination.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('eos',type=Path)
    print(json.dumps(run(p.parse_args().eos),indent=2,allow_nan=False))
