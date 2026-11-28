"""Resolve actual bath-temperature equality, without assuming chemical equilibrium."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from fd_portal_combined_entropy_split_expansion import CombinedExpansion
from fd_portal_combined_tabulated_bath_entropy import elastic_heat
from fd_portal_differential_reciprocity_reference import PositivePortal
from fd_portal_tabulated_entropy_split_expansion import Stage
from sm_dark_entropy_exchange_reference import StandardModel


def follow(s,a,b,method,tol):
    u0=np.array([.90,0.,np.log(.90)])
    S0,E0=s.entropy(0.,u0),s.energy(0.,u0)
    def rhs(x,u):
        H,y,z=s.rhs(x,u[:3])
        P=s.eos.state(u[0])[1]+s.bath(u[1:3])['pressure']
        return np.r_[y/H,z/(H*S0),3*np.exp(3*x)*P/E0]
    r=solve_ivp(rhs,(a[0],b[0]),np.r_[a[1:4],a[5:7]],method=method,rtol=tol,atol=tol*.02,dense_output=True)
    if not r.success or np.max(abs(r.y[:,-1]-np.r_[b[1:4],b[5:7]]))>8e-8:
        raise ArithmeticError('Dense path did not reproduce its exact numerical bracket')
    def delta(x):
        u=r.sol(x)
        return 1.-u[0]/np.exp(u[2])
    x=brentq(delta,a[0],b[0],xtol=2e-13)
    z=r.sol(x)
    defects=[abs(s.entropy(x,z[:3])/S0-1.-z[3]),abs(s.energy(x,z[:3])/E0-1.+z[4])]
    f,ep=s.collision(z[0],z[0],z[1])
    if (max(defects)>2e-7 or abs(delta(x))>2e-12 or z[1]>-.01
            or f[0]>=0. or ep<=0. or elastic_heat(z[0],z[0],z[1],(18,29.))!=0.):
        raise ArithmeticError('Actual equality with nonzero chemical operator failed')
    return dict(x=x,T_SM=z[0],eta_dark=z[1],timelike_net_number_and_heat=f.tolist(),
                total_local_entropy_rate=ep,absolute_ledger_defects=defects)


def run(path,archive):
    u=archive.read_bytes(); d=json.loads(u); h=hashlib.sha256(path.read_bytes()).hexdigest()
    if d['table_input_sha256']!=h or [x['operator'] for x in d['rows']]!=['timelike','combined']:
        raise ArithmeticError('Exact direct two-channel path and EOS must match')
    e=StandardModel(path)
    p,q=[PositivePortal(o) for o in ((13,36,26),(17,44,32))]
    out=[]
    for old,J in zip(d['rows'],(Stage,CombinedExpansion),strict=True):
        a=np.array(old['points'])
        i=np.flatnonzero((a[:-1,4]>0.)&(a[1:,4]<0.))
        if len(i)!=1:
            raise ArithmeticError('Need an independently observed actual thermal ordering reversal')
        br=a[i[0]:i[0]+2]
        b,c=[follow(J(e,u,n),*br,m,t) for u,n,m,t in ((p,112,'BDF',7e-10),(q,168,'Radau',3e-11))]
        k=max(abs(b[f]-c[f]) for f in ('x','T_SM','eta_dark'))
        if k>3e-7:
            raise ArithmeticError('Independent dense phase-space models could not resolve equality')
        out.append(dict(operator=old['operator'],physical_bracket=br[:,0].tolist(),
                        maximum_independent_crossing_difference=k,**c))
    return dict(rows=out,input_sha256=hashlib.sha256(u).hexdigest(),table_input_sha256=h,
                scope='Actual equality of the two numerical temperatures only, not chemical equilibrium, decoupling, or a relic abundance prediction.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('eos',type=Path)
    p.add_argument('archive',type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.archive),indent=2,allow_nan=False))
