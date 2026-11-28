"""Entropy-discrete FD collision step with positive near-equilibrium rate filtering."""
import json

import numpy as np
from scipy.optimize import root

from fd_portal_differential_reciprocity_reference import PositivePortal, integrate
from fd_portal_discrete_gradient_reservoir_reference import mean
from fd_portal_implicit_reservoir_reference import dark


def symmetric_root(A):
    if not np.allclose(A,A.T,rtol=2e-14,atol=0.):
        raise ArithmeticError('Susceptibility is not symmetric in conjugate units')
    z,V=np.linalg.eigh((A+A.T)/2)
    if min(z)<=0:
        raise ArithmeticError('Need a positive equilibrium susceptibility')
    return (V*np.sqrt(z))@V.T, (V/np.sqrt(z))@V.T


def filtered(p,T,a,dt,nodes):
    if dt<=0:
        raise ValueError('Positive timestep required')
    S=np.array([2.,T])
    y,L,_=dark(T,a,nodes)
    _,_,B=dark(T,np.zeros(2),nodes)
    C,R=symmetric_root(B/S[:,None])
    def packet(b):
        z,M,_=dark(T,b,nodes)
        g=mean(T,a,b,y,z,5,nodes)
        f,_,_,_,H=integrate(p,T,*g,return_secant=True)
        k,V=np.linalg.eigh(R@H@R)
        if min(k)<=0.:
            raise ArithmeticError('Need positive FD second-moment mobility')
        mu=2/dt*np.tanh(dt*k/2)
        v=C@(V@(mu*(V.T@(C@g))))
        return (z-y+dt*S*v)/y,M,float(g@v),float(f[2]),g,z
    b=root(lambda v:packet(v)[0],a,tol=1e-11).x
    r,M,s,p0,g,z=packet(b)
    cap=max(abs(L),abs(M),abs(dt*s),1e-200)
    entropy_error=abs(M-L+dt*s)/cap
    if np.max(abs(r))>5e-13 or entropy_error>5e-13 or np.max(abs(g-mean(T,a,b,y,z,12,nodes)))>4e-12 or s<=0.:
        raise ArithmeticError('Filtered discrete entropy or moment balance failed')
    return b,entropy_error,s/p0


def sample(p,T,nodes):
    S=np.array([2.,T])
    _,_,A=dark(T,np.zeros(2),nodes)
    _,J,_,_=integrate(p,T,0.,0.)
    C,R=symmetric_root(A/S[:,None])
    k,V=np.linalg.eigh(R@J@R)
    basis=R@V
    norm=np.linalg.norm(basis,axis=0)
    U=basis/norm
    W=np.linalg.inv(U)
    e=[]
    for i in (0,1):
        for z in (.25,4.,20.,80.):
            dt=z/k[-1]
            rows=[]
            for h in (5e-4,2.5e-4):
                a,b=[filtered(p,T,q*h*U[:,i],dt,nodes) for q in (1.,-1.)]
                d=W@(a[0]-b[0])/(2*h)
                rows.append(dict(central_factors=d.tolist(),maximum_entropy_error=max(a[1],b[1]),filtered_over_physical_entropy=[a[2],b[2]]))
            target=np.exp(-dt*k[i])
            if abs(rows[-1]['central_factors'][i]-target)>2e-6 or abs(rows[-1]['central_factors'][1-i])>2e-6:
                raise ArithmeticError('Did not match exact small-amplitude linear damping')
            e.append(dict(initial_mode=i,fast_dt_times_rate=z,exact_damping=target,samples=rows))
    return dict(T_SM=T,relaxation_rates=k.tolist(),rows=e)


def run():
    p,q=[PositivePortal(o) for o in ((13,36,26),(17,44,32))]
    rows=[]
    for T in (.05,.3,.6882487043):
        a,b=[sample(s,T,n) for s,n in ((p,112),(q,168))]
        if np.max(abs(np.asarray(a['relaxation_rates'])/b['relaxation_rates']-1.))>1e-5:
            raise ArithmeticError('Physical eigenvalues not independently resolved')
        error=max(np.max(abs(np.asarray(x['samples'][-1]['central_factors'])-y['samples'][-1]['central_factors']))
                  for x,y in zip(a['rows'],b['rows'],strict=True))
        if error>2e-6:
            raise ArithmeticError('Two meshes did not reproduce the damping')
        rows.append(dict(independent_mesh_max_absolute_amplitude_difference=error,**b))
    return dict(rows=rows,scope='Finite-step modified mobility preserves a discrete entropy identity, but its finite-step entropy is not the unfiltered physical rate. Tests are local, not an expanding history.')


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
