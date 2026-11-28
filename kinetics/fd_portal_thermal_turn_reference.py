"""Short conservative time integration near two dark-temperature nulls."""
import json

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from fd_number_heat_state_reference import evaluate, state
from fd_portal_dark_temperature_null_reference import Control


def one(C,T,r,bracket):
    f=lambda e:C.flow(T,r,float(e))['G']
    eta=brentq(f,*bracket,xtol=1e-12)
    D=T*r
    d=C.flow(T,r,eta)
    rho0,n0=d['rho'],d['n']
    U=min(n0/(2*abs(d['F'])),rho0/abs(d['Q']))
    k=evaluate(D,eta,nodes=C.nodes)
    G_eta=(f(eta+1e-4)-f(eta-1e-4))/(2e-4)
    eta_prime=-2*d['F']*U/k['dn_deta']
    curvature=-G_eta*U*eta_prime/d['cv']

    def reconstruct(n,e):
        Tnew=brentq(lambda x:state(x,n,nodes=C.nodes)['rho']-e,.85*D,1.15*D,xtol=1e-16)
        y=state(Tnew,n,nodes=C.nodes)
        return Tnew,y

    def rhs(_t,y):
        tp,v=reconstruct(y[0]*n0,y[1]*rho0)
        m=C.flow(T,tp/T,v['eta'])
        return U*np.array([-2*m['F']/n0,-m['Q']/rho0,m['sigma']*T/rho0])

    checks=[]
    for h in (.002,.001):
        xs=[]
        for sign in (-1.,1.):
            t=sign*h
            x=solve_ivp(rhs,(0.,t),[1.,1.,0.],rtol=4e-12,atol=1e-14,max_step=h/4)
            if not x.success:
                raise RuntimeError(x.message)
            n,e,si=x.y[:,-1]
            new,sv=reconstruct(n*n0,e*rho0)
            sr=((sv['entropy']-k['entropy'])*T/rho0+1-e)
            if abs(si-sr)>1e-9*max(abs(si),1e-8)+2e-13 or si*sign<=0:
                raise ValueError("Thermostat and kinetic entropy ledgers did not close.")
            xs.append((new,float(sr),float(si)))
        curv=(xs[0][0]+xs[1][0]-2*D)/(h*h)
        if curv*curvature<=0 or abs(curv/curvature-1)>8e-4:
            raise ValueError("Conservative integration did not reproduce turning curvature.")
        checks.append(dict(tau_step=h,temperature_pair=[v[0] for v in xs],
                           entropy_ledger_pair=[v[1:] for v in xs],curvature=curv))
    return dict(eta=eta,T_D=D,normalization_time=U,initial_F=d['F'],initial_Q=d['Q'],
                eta_prime=eta_prime,analytic_curvature=curvature,type='minimum' if curvature>0 else 'maximum',
                finite_conservative_checks=checks)


if __name__=='__main__':
    C=Control(True)
    T,r=.05,.80
    d=[one(C,T,r,b) for b in ((1.9,2.5),(3.45,4.1))]
    if [a['type'] for a in d]!=['minimum','maximum']:
        raise ValueError("Unexpected turn directions.")
    print(json.dumps(dict(T_SM=T,T_D_ratio=r,scope="Short perturbative thermostat comparison, no FRW or evolved visible temperature.",
                          turns=d),indent=2,allow_nan=False))
