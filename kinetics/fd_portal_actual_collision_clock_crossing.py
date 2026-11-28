"""Locate local fixed-energy collision-rate/H crossings on the exact direct path."""
import argparse
import hashlib
import json
from pathlib import Path
from functools import lru_cache

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from fd_portal_combined_entropy_split_expansion import CombinedExpansion
from fd_portal_differential_reciprocity_reference import PositivePortal
from fd_portal_live_fixed_energy_collision_spectrum import eigen,fixed_energy,jac
from fd_portal_tabulated_entropy_split_expansion import Stage
from sm_dark_entropy_exchange_reference import StandardModel


def dense(s,lo,hi,method,tol):
    orig=np.array([.90,0.,np.log(.90)])
    S0,E0=s.entropy(0.,orig),s.energy(0.,orig)
    N0=s.bath(orig[1:])['n']
    def rhs(x,u):
        H,F,G=s.rhs(x,u[:3])
        T,eta,l=u[:3]
        p=s.eos.state(T)[1]+s.bath(u[1:3])['pressure']
        q,_=s.collision(T,np.exp(l),eta)
        return np.r_[F/H,G/(S0*H),3*np.exp(3*x)*p/E0,-np.exp(3*x)*q[0]/(N0*H)]
    a,b=np.r_[lo[1:4],lo[5:8]],np.r_[hi[1:4],hi[5:8]]
    r=solve_ivp(rhs,(lo[0],hi[0]),a,method=method,rtol=tol,atol=.022*tol,max_step=.002,dense_output=True)
    if not r.success or max(abs(r.y[:,-1]-b))>6e-8:
        raise ArithmeticError('Dense solution did not recover the archived right endpoint')
    def defects(x,u):
        n=s.bath(u[1:3])['n']
        return np.array([s.entropy(x,u[:3])/S0-1.-u[3],s.energy(x,u[:3])/E0-1.+u[4],np.exp(3*x)*n/N0-1.-u[5]])
    return r.sol,defects


def locate(s,mode,lo,hi,method,tol,h):
    i={'slow':0,'fast':1}[mode]
    f,L=dense(s,lo,hi,method,tol)
    @lru_cache
    def rates(x):
        v,_,_=fixed_energy(s,f(x)[:3])
        z=eigen(jac(v,h))
        if min(z.real)<=0. or max(abs(z.imag))>1e-9:
            raise ArithmeticError('Selected local spectrum is not two positive real rates')
        return z.real

    def g(x):
        return rates(float(x))[i]-1.
    a,b=float(lo[0]),float(hi[0])
    probe=np.linspace(a,b,5)
    z=[g(p) for p in probe]
    if not z[0]>0.>z[-1] or np.any(np.diff(z)>=0.):
        raise ArithmeticError(('Missing ordered crossing in archived bracket',mode,[a,b],z))
    x=brentq(g,a,b,xtol=4e-13)
    u=f(x)
    defect=L(x,u)
    if max(abs(defect))>2e-7 or abs(g(x))>2e-10:
        raise ArithmeticError('Root or original direct ledgers failed')
    return dict(x=x,state=u[:3].tolist(),rates_over_H=rates(x).tolist(),
                bracket_sample_rate_minus_one=z,absolute_ledger_defects=abs(defect).tolist())


def run(path,early,late,reference):
    raw=[p.read_bytes() for p in (early,late,reference)]
    hashes=[hashlib.sha256(r).hexdigest() for r in raw]
    a,b,c=map(json.loads,raw)
    h=hashlib.sha256(path.read_bytes()).hexdigest()
    if any(t['table_input_sha256']!=h for t in (a,b,c)) or c['prior_continuation_and_clocks_sha256'][:2]!=hashes[:2]:
        raise ArithmeticError('Need exact direct prefix and original temperature reference labels')
    E=StandardModel(path)
    P,Q=[PositivePortal(m) for m in ((13,36,26),(17,44,32))]
    rows=[]
    for kind,K,early_record,late_record in zip(('timelike','combined'),(Stage,CombinedExpansion),a['rows'],b['rows'],strict=True):
        if early_record['operator']!=kind or late_record['operator']!=kind:
            raise ArithmeticError('Operator mismatch')
        for mode,(u,v) in (('slow',(.50,.55)),('fast',(.85,.91))):
            # The older equilibrium reference for the combined fast case lies above .91.
            # The smaller instantaneous live rate is expected to cross earlier; fail explicitly otherwise.
            Pts=np.asarray((early_record if mode=='slow' else late_record)['points'])
            lo,hi=[next(w for w in Pts if w[0]==q) for q in (u,v)]
            ref=next(r for r in c['rows'] if (r['operator'],r['mode'])==(kind,mode))
            A,B=[locate(K(E,k,m),mode,lo,hi,alg,tol,eps) for k,m,alg,tol,eps in (
                    (P,112,'BDF',1e-10,6e-4),(Q,168,'Radau',4e-12,3e-4))]
            err=max(abs(A['x']-B['x']),max(abs(np.asarray(A['state'])-B['state'])))
            if err>4e-7:
                raise ArithmeticError('Independent numerical histories disagree on local crossing')
            rows.append(dict(operator=kind,mode=mode,bracket=[u,v],reference_log_a=ref['log_a'],
                reference_visible_T=ref['local_crossing_T'],live_minus_reference_log_a=B['x']-ref['log_a'],
                live_minus_reference_visible_T=B['state'][0]-ref['local_crossing_T'],independent_root_defect=err,**B))
    return dict(rows=rows,table_input_sha256=h,input_sha256=hashes,
                scope='Local unity crossing of collision-only fixed-total-energy Jacobian rate divided by H along a restricted direct history; not a freeze-out threshold or integrated expanding-system stability test.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('path','early','late','reference'):
        p.add_argument(key,type=Path)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
