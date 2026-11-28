"""Local equilibrium slopes forced by adiabatic expansion of tabulated FD baths."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from boost_axis_elastic_heat_reference import integrate as elastic, leptons
from fd_portal_differential_reciprocity_reference import PositivePortal, integrate
from fd_portal_reservoir_exponential_filter_reference import symmetric_root
from fd_portal_tabulated_entropy_split_expansion import Stage, hubble
from sm_dark_entropy_exchange_reference import StandardModel


def parcel(eos,p,n,T):
    s=Stage(eos,p,n)
    a=np.array([T,0.,np.log(T)])
    d=s.bath(a[1:])
    rs,ps,_,cs=eos.state(T)
    rn,pn,N,B,E,C=[d[k] for k in ('rho','pressure','n','dn_deta','chemical_energy','cv_fixed_number')]
    H=hubble(rs,rn)
    A=np.array([[B,B*E/T],[B*E,T*C+B*E*E/T]])
    Tprime=-3*(rs+ps)/cs
    eta,lp=np.linalg.solve(A,-3*np.array([N,rn+pn]))
    force=np.array([eta,-Tprime/T**2+lp/T])
    def g(dx):
        t,u,v=s.drift_chart(a,dx)
        return np.array([u,1./t-np.exp(-v)])
    eps=2e-6
    check=(g(eps)-g(-eps))/(2*eps)
    norm=np.array([1.,T])
    drift=np.max(abs((force-check)*norm))
    if drift>3e-8:
        raise ArithmeticError('Independent adiabatic drift does not reproduce expansion forcing')

    X=np.diag([1.,1./T])@np.linalg.inv(A)+np.diag([0.,1/(cs*T*T)])
    R,C=symmetric_root(X)
    _,J,_,_=integrate(p,T,0.,0.)
    Q=np.diag([2.,T])
    M=Q@J@Q
    order,cut=(18,29.) if n==112 else (25,36.)
    el=T*sum(elastic(T,T,0.,l,n=order,cut=cut,linear=True) for l in leptons)
    rows=[]
    for tag,vis in (('timelike',M),('combined',M+np.diag([0.,el]))):
        k,V=np.linalg.eigh(R@vis@R)
        b=V.T@C@force
        modes=(norm[:,None]*(R@V))*((H/k)*b)[None,:]
        lag=modes.sum(axis=1)
        residual=norm*(X@vis@(lag/norm)/H-force)
        if np.max(abs(residual))>2e-11:
            raise ArithmeticError('Frozen forcing failed its relaxation balance')
        rows.append(dict(operator=tag,positive_mode_rates_over_H=(k/H).tolist(),
                         forcing_projection=b.tolist(),modal_response_eta_and_beta=modes.T.tolist(),
                         algebraic_eta_and_beta=lag.tolist(),
                         small_force_candidate=bool(np.max(abs(lag))<.02),response_equation_defect=np.max(abs(residual))))
    return dict(temperature=T,adiabatic_eta_and_beta_slopes=(force*norm).tolist(),
                independently_measured_drift_defect=drift,operators=rows)


def run(path):
    eos=StandardModel(path)
    p,q=[PositivePortal(o) for o in ((13,36,26),(17,44,32))]
    rows=[]
    for T in (.05,.10,.15,.20,.30,.45,.55,.64135,.6882487043,.8,.95):
        a,b=[parcel(eos,x,n,T) for x,n in ((p,112),(q,168))]
        diff=max(np.max(abs(np.asarray(u['algebraic_eta_and_beta'])-v['algebraic_eta_and_beta'])/
                        np.maximum(1e-15,abs(np.asarray(v['algebraic_eta_and_beta']))))
                 for u,v in zip(a['operators'],b['operators'],strict=True))
        if diff>4e-5:
            raise ArithmeticError('Independent phase-space comparison did not resolve forced response')
        rows.append(dict(max_relative_mesh_lag_difference=diff,**b))
    return dict(rows=rows,table_input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                scope='Algebraic local response at common equilibrium, not a solved expanding lag, freeze-out, or relic abundance. Large-force rows diagnose limited validity only.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('eos',type=Path)
    print(json.dumps(run(p.parse_args().eos),indent=2,allow_nan=False))
