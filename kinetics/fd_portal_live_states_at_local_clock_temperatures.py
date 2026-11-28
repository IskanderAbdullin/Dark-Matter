"""Measure actual numerical states where the visible temperature crosses local FD clock roots."""
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


def measure(s,lo,hi,T,method,tol):
    u=np.array([.90,0.,np.log(.90)])
    S0,E0=s.entropy(0.,u),s.energy(0.,u)
    N0=s.bath(u[1:])['n']
    def flow(x,v):
        H,z,q=s.rhs(x,v[:3])
        P=s.eos.state(v[0])[1]+s.bath(v[1:3])['pressure']
        return np.r_[z/H,q/(H*S0),3*np.exp(3*x)*P/E0]

    a=solve_ivp(flow,(lo[0],hi[0]),np.r_[lo[1:4],lo[5:7]],method=method,
                rtol=tol,atol=tol*.02,dense_output=True)
    if not a.success or np.max(abs(a.y[:,-1]-np.r_[hi[1:4],hi[5:7]]))>5e-8:
        raise ArithmeticError('Dense comparison did not reproduce the physical bracket endpoint')
    xx=brentq(lambda x:a.sol(x)[0]-T,lo[0],hi[0],xtol=3e-13)
    b=a.sol(xx)
    rate=flow(xx,b)
    n=s.bath(b[1:3])['n']
    norm=np.exp(3*xx)*n/N0
    defect=max(abs(s.entropy(xx,b[:3])/S0-1.-b[3]),abs(s.energy(xx,b[:3])/E0-1.+b[4]))
    if (abs(b[0]-T)>3e-12 or rate[0]>=0. or defect>2e-7 or b[3]<=0.):
        raise ArithmeticError('Actual local-clock crossing or physical budget did not close')
    return dict(log_a=xx,actual_eta=b[1],actual_beta=1.-b[0]/np.exp(b[2]),comoving_count_over_initial=norm,
                combined_total_ledger_defect=defect,visible_temperature_derivative=rate[0])


def run(path,prior,final,clocks):
    a,b,c=[p.read_bytes() for p in (prior,final,clocks)]
    d,e,k=[json.loads(z) for z in (a,b,c)]
    h=hashlib.sha256(path.read_bytes()).hexdigest()
    if (any(p['table_input_sha256']!=h for p in (d,e,k))
            or e['restart_and_clock_sha256']!=[hashlib.sha256(z).hexdigest() for z in (a,c)]):
        raise ArithmeticError('Actual path join and local clock input must be exact')
    p,q=[PositivePortal(o) for o in ((13,36,26),(17,44,32))]
    eos=StandardModel(path)
    out=[]
    for left,right,K in zip(d['rows'],e['rows'],(Stage,CombinedExpansion),strict=True):
        tag=left['operator']
        if right['operator']!=tag:
            raise ArithmeticError('Separate physical trajectories cannot be interchanged')
        u,v=[np.asarray(z['points']) for z in (left,right)]
        if np.max(abs(u[-1]-v[0]))>2e-15:
            raise ArithmeticError('Wrong continuation state')
        r=np.r_[u,v[1:]]
        if np.any(np.diff(r[:,1])>=0.):
            raise ArithmeticError('No monotone measured temperature bracket')
        for x in k['rows']:
            if tag!=x['operator']:
                continue
            T=x['crossing_T']
            hits=np.flatnonzero((r[:-1,1]>=T)&(r[1:,1]<=T))
            if len(hits)!=1:
                raise ArithmeticError('Local clock temperature not crossed in the actual data')
            m=hits[0]
            a1,b1=[measure(K(eos,z,n),r[m],r[m+1],T,method,tol)
                  for z,n,method,tol in ((p,112,'BDF',5e-9),(q,168,'Radau',5e-10))]
            delta=max(abs(a1[f]-b1[f]) for f in ('log_a','actual_eta','actual_beta','comoving_count_over_initial'))
            if delta>3e-7:
                raise ArithmeticError('Independent dense paths disagree on the sampled clock state')
            out.append(dict(operator=tag,mode=x['mode'],local_crossing_T=T,bracket_log_a=r[m:m+2,0].tolist(),
                            independent_maximum_difference=delta,instantaneous_algebraic_eta_and_beta=x['algebraic_eta_and_beta'],**b1))
    return dict(rows=out,table_input_sha256=h,prior_continuation_and_clocks_sha256=[hashlib.sha256(z).hexdigest() for z in (a,b,c)],
                scope='Just a label on a smooth, finite restricted FD trajectory; not a decoupling time or an equilibrium transition.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('eos','prior','final','clocks'):
        p.add_argument(k,type=Path)
    t=p.parse_args()
    print(json.dumps(run(t.eos,t.prior,t.final,t.clocks),indent=2,allow_nan=False))
