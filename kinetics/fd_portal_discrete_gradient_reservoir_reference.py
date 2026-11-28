"""Symmetric discrete-gradient two-moment FD update at fixed visible temperature."""
import json
from functools import lru_cache

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root
from scipy.special import roots_legendre

from fd_portal_differential_reciprocity_reference import PositivePortal,integrate
from fd_portal_implicit_reservoir_reference import dark


@lru_cache(None)
def rule(n):
    u,v=roots_legendre(n)
    return (u+1)/2,v/2


def recover(T,y,guess,dn):
    u=root(lambda a:(dark(T,a,dn)[0]-y)/y,guess,
           jac=lambda a:dark(T,a,dn)[2]/y[:,None],tol=2e-12).x
    state,L,_=dark(T,u,dn)
    if min(y)<=0. or np.max(abs(state/y-1.))>5e-13:
        raise ArithmeticError('No resolved positive thermodynamic inverse')
    return u


def mean(T,a,b,y,z,n,dn):
    t,w=rule(n)
    u=np.array([recover(T,y+s*(z-y),a+s*(b-a),dn) for s in t])
    return w@u


def step(p,T,a,dt,n,dn):
    y,L,_=dark(T,a,dn)
    def packet(b):
        z,M,_=dark(T,b,dn)
        g=mean(T,a,b,y,z,n,dn)
        F=integrate(p,T,*g)[0]
        return (z-y+dt*F[:2]*[2.,T])/y,z,M,F,g
    b=root(lambda v:packet(v)[0],a,tol=5e-12).x
    r,z,M,f,g=packet(b)
    other=mean(T,a,b,y,z,12,dn)
    cap=max(abs(M),abs(L),abs(dt*f[2]),1e-200)
    err=abs(M-L+dt*f[2])/cap
    pairing=abs((g/[2.,T])@(z-y)-(M-L))/cap
    if np.max(abs(r))>5e-13 or np.max(abs(other-g))>4e-12 or max(err,pairing)>5e-13 or f[2]<=0.:
        raise ArithmeticError('Independent path, rates and discrete entropy identity failed')
    return b,max(err,pairing),dt*f[2],np.max(abs(r))


def trajectory(p,T,old,tau,steps,dn):
    b=old.copy()
    v=[]
    for _ in range(steps):
        b,e,f,r=step(p,T,b,tau/steps,5,dn)
        v.append((e,f,r))
    y0,L,_=dark(T,old,dn)
    y1,M,_=dark(T,b,dn)
    z=np.asarray(v)
    if abs(M-L+z[:,1].sum())>2e-12*max(abs(M),abs(L)):
        raise ArithmeticError('Telescoped discrete entropy identity failed')
    return b,y1,np.max(z[:,0]),np.max(z[:,2])


def run():
    p,q=[PositivePortal(o) for o in ((13,36,26),(17,44,32))]
    rows=[]
    for T,old in ((.05,np.array([.6,.16])),(.3,np.array([-1.5,.16])),(.6882487043,np.array([1.8,-.16]))):
        y,_,_=dark(T,old,168)
        f=integrate(q,T,*old)[0]
        tau=.12/max(abs(f[:2]*[2.,T]/y))
        def rhs(t,u):
            J=dark(T,u,168)[2]
            f=integrate(q,T,*u)[0]
            return np.linalg.solve(J,-tau*f[:2]*[2.,T])
        fit=solve_ivp(rhs,(0.,1.),old,method='DOP853',rtol=3e-12,atol=4e-14)
        if not fit.success: raise ArithmeticError('Missing independent high-accuracy reference')
        expected=dark(T,fit.y[:,-1],168)[0]
        calc=[]
        for n in (2,4,8):
            z,m,e,r=trajectory(q,T,old,tau,n,168)
            calc.append(dict(steps=n,force=z.tolist(),error_over_initial_state_scale=np.max(abs(m-expected)/y),
                             max_entropy_error=e,max_rate_residual=r))
        ratio=np.array([calc[j]['error_over_initial_state_scale']/calc[j+1]['error_over_initial_state_scale'] for j in (0,1)])
        if np.min(ratio)<3.3 or np.max(ratio)>4.8:
            raise ArithmeticError('Second-order path control not resolved')
        x,_,_,_=trajectory(p,T,old,tau,4,112)
        if np.max(abs(x-calc[1]['force']))>2e-7:
            raise ArithmeticError('Independent phase-space mesh discrepancy')
        rows.append(dict(T_SM=T,initial_forces=old.tolist(),total_local_time=tau,
                         independent_DOP853_force=fit.y[:,-1].tolist(),levels=calc,error_halving_ratios=ratio.tolist(),
                         independent_mesh_four_step_max_force_difference=np.max(abs(x-calc[1]['force']))))
    return dict(rows=rows,scope='Symmetric local collision update on thermodynamic moments for a fixed reservoir, not a cosmological trajectory or global convergence/positivity guarantee.')


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
