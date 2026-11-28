"""Bracket accumulated entropy equality for two distinct direct FD histories."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from fd_portal_combined_entropy_split_expansion import CombinedExpansion
from fd_portal_differential_reciprocity_reference import PositivePortal
from fd_portal_tabulated_entropy_split_expansion import Stage
from sm_dark_entropy_exchange_reference import StandardModel


def solve(stage,old,upper,method,atol):
    T=.90
    S=stage.entropy(0.,[T,0.,np.log(T)])
    x0,Ts,eta,ld,_,P,_,_=old
    u0=[Ts/T,eta,ld-np.log(T),P]
    def f(x,z):
        H,k,s=stage.rhs(x,[T*z[0],z[1],z[2]+np.log(T)])
        return np.r_[k/[T,1.,1.],s/S]/H
    u=solve_ivp(f,(x0,upper),u0,method=method,dense_output=True,rtol=atol*40,atol=atol)
    if not u.success:
        raise ArithmeticError('No independently solved continuous interval')
    return u.sol,f


def run(path,base,follow):
    y,z=[x.read_bytes() for x in (base,follow)]
    d,u=json.loads(y),json.loads(z)
    H=hashlib.sha256(path.read_bytes()).hexdigest()
    if H!=d['table_input_sha256'] or H!=u['table_input_sha256'] or hashlib.sha256(y).hexdigest()!=u['archive_sha256']:
        raise ArithmeticError('Exact joined entropy histories and tabulation are required')
    D=u['combined_minus_timelike']
    if not D[0][5]<0.<D[1][5]:
        raise ArithmeticError('Expected an observed numerical sign difference')
    starts=[v['controlled_full_path'][-1] for v in d['rows']]
    if [x[0] for x in starts]!=[.25,.25] or [a['points'][1][0] for a in u['rows']]!=[.30,.30]:
        raise ArithmeticError('Wrong exact archived comparison interval')
    eos=StandardModel(path)
    p,q=[PositivePortal(o) for o in ((13,36,26),(17,44,32))]
    out=[]
    for port,n,method,tol in ((p,112,'BDF',1e-13),(q,168,'Radau',2e-14)):
        a,b=[solve(k(eos,port,n),v,.30,method,tol) for k,v in zip((Stage,CombinedExpansion),starts,strict=True)]
        def diff(x):
            return b[0](x)[3]-a[0](x)[3]
        grid=np.linspace(.25,.30,41)
        f=np.array([diff(x) for x in grid])
        ix=np.flatnonzero(f[1:]*f[:-1]<0.)
        if f[0]>=0. or f[-1]<=0. or len(ix)!=1 or abs(f[-1]-D[1][5])>3e-12:
            raise ArithmeticError('Interpolated interval did not reproduce the archived sign comparison')
        j=ix[0]
        bound=grid[j:j+2]
        c=brentq(diff,*bound,xtol=2e-13)
        s,t=[o[0](c) for o in (a,b)]
        derivative=b[1](c,t)[3]-a[1](c,s)[3]
        if derivative<=0. or abs(diff(c))>4e-17:
            raise ArithmeticError('Continuous accumulated entropy crossing not resolved')
        out.append(dict(method=method,log_a_crossing=c,positive_local_difference_derivative=derivative,
                        sampled_bracket=bound.tolist(),interval_difference_endpoints=f[[0,-1]].tolist(),
                        visible_temperatures=[.90*s[0],.90*t[0]]))
    shift=abs(out[0]['log_a_crossing']-out[1]['log_a_crossing'])
    if shift>4e-5:
        raise ArithmeticError('Independent histories did not resolve the local cumulative crossing')
    return dict(rows=out,maximum_cross_mesh_log_a_shift=shift,table_input_sha256=H,
                base_and_follow_sha256=[hashlib.sha256(v).hexdigest() for v in (y,z)],
                scope='Only the difference of two cumulative entropy histories crosses zero. No reaction threshold, entropy violation, or cosmological prediction.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('eos','base','follow'):
        p.add_argument(k,type=Path)
    t=p.parse_args()
    print(json.dumps(run(t.eos,t.base,t.follow),indent=2,allow_nan=False))
