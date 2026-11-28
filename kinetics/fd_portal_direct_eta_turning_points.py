"""Locate a dark reduced chemical-potential minimum on the direct FD history."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from fd_portal_combined_entropy_split_expansion import CombinedExpansion
from fd_portal_differential_reciprocity_reference import PositivePortal
from fd_portal_tabulated_entropy_split_expansion import Stage, hubble
from sm_dark_entropy_exchange_reference import StandardModel


def decomposition(s,u):
    T,eta,l=u[:3]
    td=np.exp(l)
    d=s.bath(u[1:3])
    n,r,p,B,g,C=[d[k] for k in ('n','rho','pressure','dn_deta','chemical_energy','cv_fixed_number')]
    rho,pressure,_,cv=s.eos.state(T)
    H=hubble(rho,r)
    q,sigma=s.collision(T,td,eta)
    M=np.array([[B,B*g/td],[B*g,td*C+B*g*g/td]])
    vals=np.linalg.solve(M,np.array([[-3*n,-q[0]/H,0.],[-3*(r+p),0.,-q[1]/H]]))
    return vals,np.array([(-3*(rho+pressure)+q[1]/H)/cv,*vals.sum(axis=1)]),q,sigma,H


def locate(s,lo,hi,method,tol):
    orig=np.array([.90,0.,np.log(.90)])
    S0,E0=s.entropy(0.,orig),s.energy(0.,orig)
    N0=s.bath(orig[1:])['n']
    def rhs(x,u):
        _,D,q,sigma,H=decomposition(s,u)
        p=s.eos.state(u[0])[1]+s.bath(u[1:3])['pressure']
        a=np.exp(3*x)
        return np.r_[D,a*sigma/(S0*H),3*a*p/E0,-a*q[0]/(N0*H)]
    a,b=np.r_[lo[1:4],lo[5:8]],np.r_[hi[1:4],hi[5:8]]
    r=solve_ivp(rhs,(lo[0],hi[0]),a,rtol=tol,atol=.02*tol,method=method,dense_output=True,max_step=.001)
    if not r.success or max(abs(r.y[:,-1]-b))>2e-8:
        raise ArithmeticError('Independent dense integration fails exact archive terminal comparison')
    def f(x):
        return rhs(x,r.sol(x))[1]
    span=[f(lo[0]),f(hi[0])]
    if not span[0]<0.<span[1]:
        raise ArithmeticError('Expected negative-to-positive eta derivative')
    x=brentq(f,lo[0],hi[0],xtol=2e-14)
    u=r.sol(x)
    parts,deriv,q,sigma,H=decomposition(s,u)
    report=rhs(x,u)
    L=[abs(s.entropy(x,u[:3])/S0-1.-u[3]),abs(s.energy(x,u[:3])/E0-1.+u[4]),
       abs(np.exp(3*x)*s.bath(u[1:3])['n']/N0-1.-u[5])]
    if max(L)>2e-7 or abs(deriv[1])>2e-12 or not (q[0]<0. and sigma>0.):
        raise ArithmeticError('Turning point does not satisfy exact ledger and nonzero source gates')
    return dict(x=x,T_SM=u[0],eta=u[1],T_dark=np.exp(u[2]),rate_eta_bracket=span,
                eta_derivative_pieces=parts[0].tolist(),normalized_count_source_per_log_a=report[-1],
                normalized_entropy_source_per_log_a=report[-3],ledger_defects=L)


def run(eos_path,earlier,later):
    a,b=[p.read_bytes() for p in (earlier,later)]
    A,B=json.loads(a),json.loads(b)
    h=hashlib.sha256(eos_path.read_bytes()).hexdigest()
    if any(d['table_input_sha256']!=h for d in (A,B)) or hashlib.sha256(a).hexdigest()!=B['previous_segment_sha256']:
        raise ArithmeticError('The two direct archives do not form the same continuation')
    eos=StandardModel(eos_path)
    p,q=[PositivePortal(o) for o in ((13,36,26),(17,44,32))]
    out=[]
    for v,w,K in zip(A['rows'],B['rows'],(Stage,CombinedExpansion),strict=True):
        if v['operator']!=w['operator'] or not np.allclose(v['points'][-1],w['points'][0],rtol=0.,atol=2e-15):
            raise ArithmeticError('Direct path does not preserve exact .96 reference state')
        s1,s2=K(eos,p,112),K(eos,q,168)
        points=np.array(v['points'][-4:]+w['points'][1:])
        rates=np.array([decomposition(s2,u[1:4])[1][1] for u in points])
        ix=np.flatnonzero((rates[:-1]<0.)&(rates[1:]>0.))
        if len(ix)!=1:
            raise ArithmeticError('No unique on-path derivative sign bracket in the sampled late range')
        lo,hi=points[ix[0]:ix[0]+2]
        c,d=locate(s1,lo,hi,'BDF',1e-11),locate(s2,lo,hi,'Radau',3e-12)
        defect=max(abs(c[k]-d[k]) for k in ('x','T_SM','T_dark','eta'))
        if defect>4e-7:
            raise ArithmeticError('Turning point unresolved across independent methods and quadratures')
        out.append(dict(operator=v['operator'],archived_bracket=[lo[0],hi[0]],maximum_method_defect=defect,**d))
    return dict(rows=out,eta_piece_fields=['adiabatic_expansion','number_exchange','energy_exchange'],
                input_sha256=[hashlib.sha256(z).hexdigest() for z in (a,b)],table_input_sha256=h,
                scope='Minimum of mu_dark/T_dark on restricted direct histories. Not a zero net source or a physical freeze-out label.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('eos','earlier','later'):
        p.add_argument(k,type=Path)
    z=p.parse_args()
    print(json.dumps(run(z.eos,z.earlier,z.later),indent=2,allow_nan=False))
