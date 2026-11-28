"""Resolve the physical near-equilibrium stiff modes of two local FD steppers."""
import json

import numpy as np

from fd_portal_differential_reciprocity_reference import PositivePortal, integrate
from fd_portal_discrete_gradient_reservoir_reference import step as symmetric
from fd_portal_implicit_reservoir_reference import dark, step as backward


def study(p,T,nodes):
    zero=np.zeros(2)
    _,_,A=dark(T,zero,nodes)
    _,J,_,_=integrate(p,T,*zero)
    B=np.linalg.solve(A,np.diag([2.,T])@J)
    eig,V=np.linalg.eig(B)
    idx=np.argsort(eig)
    eig,V=eig[idx],V[:,idx]
    if not (np.isrealobj(eig) and np.all(eig>0.)):
        raise ArithmeticError('Physical equilibrium linearization is not a relaxation')
    W=np.linalg.inv(V)
    if np.linalg.cond(V)>80.:
        raise ArithmeticError('Ill-conditioned spectral projection')
    output=[]
    for kapp in (0.25,4.,20.,80.):
        dt=kapp/eig[-1]
        e=[]
        for solver in ('discrete_gradient','backward_euler'):
            agreement=[]
            cross=[]
            for h in (5e-4,2.5e-4):
                pkt=[]
                for sign in (1,-1):
                    x=sign*h*V[:,-1]
                    if solver=='discrete_gradient':
                        z,*_=symmetric(p,T,x,dt,5,nodes)
                    else:
                        z=backward(p,T,x,nodes,dt)['final_forces']
                    pkt.append(z)
                modal=W@(np.asarray(pkt[0])-pkt[1])/(2*h)
                agreement.append(modal[-1])
                cross.append(abs(modal[0]))
            expected=(1.-.5*kapp)/(1.+.5*kapp) if solver=='discrete_gradient' else 1./(1.+kapp)
            if abs(agreement[-1]-expected)>2e-6 or cross[-1]>2e-6 or abs(agreement[-1]-expected)>abs(agreement[0]-expected)/2.5+1e-10:
                raise ArithmeticError('Nonlinear small-amplitude responses do not resolve spectral prediction')
            e.append(dict(method=solver,expected_factor=expected,central_amplitude_factors=agreement,cross_mode_amplitudes=cross))
        output.append(dict(fast_mode_dt_times_rate=kapp,damping=e))
    return dict(T_SM=T,positive_rate_eigenvalues=eig.tolist(),force_basis_condition=np.linalg.cond(V),steps=output)


def run():
    p,q=[PositivePortal(o) for o in ((13,36,26),(17,44,32))]
    rows=[]
    for T in (.05,.3,.6882487043):
        a,b=[study(r,T,n) for r,n in ((p,112),(q,168))]
        drift=np.max(abs(np.asarray(a['positive_rate_eigenvalues'])/b['positive_rate_eigenvalues']-1.))
        if drift>1e-5:
            raise ArithmeticError('Independent phase-space mesh failed equilibrium rate check')
        rows.append(dict(relative_mesh_rate_discrepancy=drift,**b))
    return dict(rows=rows,scope='Local two-bath FD near-equilibrium numerical damping, not a global timestep or expanding-history criterion.')


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
