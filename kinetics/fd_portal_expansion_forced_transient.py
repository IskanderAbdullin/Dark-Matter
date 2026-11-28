"""Compare local lag forecasts with an actual short moving-equilibrium history."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq

from fd_portal_differential_reciprocity_reference import PositivePortal,integrate
from fd_portal_tabulated_entropy_split_expansion import Stage,hubble
from sm_dark_entropy_exchange_reference import StandardModel


def study(eos,p,n,T):
    s=Stage(eos,p,n)
    a=np.array([T,0.,np.log(T)])
    S0=s.entropy(0.,a)
    grid=np.array([0.,.004,.012,.025,.05])

    def td(T):
        return s.bath([0.,np.log(T)])

    def equilibrium(x):
        return brentq(lambda v:(eos.state(v)[2]+td(v)['entropy'])/S0-np.exp(-3*x),
                      T*.92,T*(1.+2e-14),xtol=4e-15)

    def coefficients(t):
        d=td(t)
        es,ps,_,cs=eos.state(t)
        N,ed,pd,B,E,C=(d[k] for k in ('n','rho','pressure','dn_deta','chemical_energy','cv_fixed_number'))
        A=np.array([[B,B*E/t],[B*E,t*C+B*E*E/t]])
        k,j=np.linalg.solve(A,-3*np.array([N,ed+pd]))
        f=np.array([k,3*(es+ps)/(cs*t*t)+j/t])
        H0=np.diag([1.,1./t])@np.linalg.inv(A)+np.diag([0.,1/(cs*t*t)])
        _,J,_,_=integrate(p,t,0.,0.)
        Q=np.diag([2.,t])
        L=H0@Q@J@Q/hubble(es,ed)
        return np.r_[f,L.ravel()]

    def physical(x,v):
        H,b,u=s.rhs(x,v[:3])
        return np.r_[b/H,u/(S0*H)]

    sol=solve_ivp(physical,(0.,grid[-1]),np.r_[a,0.],t_eval=grid,rtol=2e-11,atol=3e-14,method='DOP853')
    if not sol.success:
        raise ArithmeticError('Actual short differential trajectory failed')
    y=sol.y.T
    control=np.max([abs(s.entropy(x,u[:3])/S0-1.-u[3]) for x,u in zip(grid,y,strict=True)])
    if control>2e-10:
        raise ArithmeticError('Actual trajectory lost comoving entropy closure')

    approximate=[]
    for nodes in (25,49):
        x=np.linspace(0.,grid[-1],nodes)
        C=CubicSpline(x,np.array([coefficients(equilibrium(u)) for u in x]))
        def f(t,u):
            d=C(t)
            return d[:2]-d[2:].reshape(2,2)@u
        z=solve_ivp(f,(0.,grid[-1]),np.zeros(2),t_eval=grid,rtol=2e-11,atol=3e-14,method='DOP853')
        if not z.success:
            raise ArithmeticError('Linear modal forcing trajectory failed')
        approximate.append(z.y.T)

    if np.max(abs(approximate[0]-approximate[1]))>2e-7:
        raise ArithmeticError('Movement of local coefficients was not spatially resolved')
    rows=[]
    for x,v,w in zip(grid,y,approximate[1],strict=True):
        te=equilibrium(x)
        f=coefficients(te)
        stationary=np.linalg.solve(f[2:].reshape(2,2),f[:2]) * [1.,te]
        u=np.array([v[1],1.-v[0]/np.exp(v[2])])
        o=w*[1.,te]
        rows.append(dict(log_a=x,common_entropy_temperature=te,nonlinear_eta_and_beta=u.tolist(),
                         linear_transient_eta_and_beta=o.tolist(),instantaneous_algebraic_eta_and_beta=stationary.tolist(),
                         normalized_comoving_entropy=v[3]))
    return dict(initial_T=T,continuous_entropy_defect=control,coefficient_chart_discrepancy=np.max(abs(approximate[0]-approximate[1])),
                rows=rows)


def run(path):
    eos=StandardModel(path)
    p,q=[PositivePortal(o) for o in ((13,36,26),(17,44,32))]
    rows=[]
    for t in (.55,.6882487043):
        a,b=[study(eos,z,n,t) for z,n in ((p,112),(q,168))]
        difference=max(np.max(abs(np.asarray(x['nonlinear_eta_and_beta'])-y['nonlinear_eta_and_beta']))
                       for x,y in zip(a['rows'],b['rows'],strict=True))
        if difference>2e-7:
            raise ArithmeticError('Two nonlinear phase-space controls disagree')
        rows.append(dict(maximum_two_mesh_live_force_difference=difference,**b))
    return dict(rows=rows,table_input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                scope='Short homogeneous, vacuum timelike FD histories started at exact equilibrium. Not a full expansion history or a decoupling threshold.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('eos',type=Path)
    print(json.dumps(run(p.parse_args().eos),indent=2,allow_nan=False))
