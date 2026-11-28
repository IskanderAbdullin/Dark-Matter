"""Total-entropy discrete positive filter for closed tabulated FD collision modes."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from fd_portal_differential_reciprocity_reference import PositivePortal, integrate
from fd_portal_reservoir_exponential_filter_reference import symmetric_root
from fd_portal_tabulated_bath_discrete_entropy import Tabulated
from sm_dark_entropy_exchange_reference import StandardModel


class Filtered(Tabulated):
    def __init__(self,T,u,dn,eos):
        super().__init__(T,u,dn,eos)
        def fixed(x):
            return eos.state(x)[0]+self.dark([0.,np.log(x)])[0][1]-self.E
        self.equilibrium_T=brentq(fixed,self.lo,self.hi,xtol=4e-15)
        x=self.equilibrium_T
        A=self.dark([0.,np.log(x)])[2]
        H=np.diag([1.,1./x])@np.linalg.inv(A)+np.diag([0.,1/(x*x*eos.state(x)[3])])
        self.R,self.C=symmetric_root(H)

    def matrices(self,p):
        T=self.equilibrium_T
        S=np.diag([2.,T])
        _,J,_,_=integrate(p,T,0.,0.)
        return self.R@S@J@S@self.R

    def step_flux(self,p,u,dt):
        _,_,g,T,_=self.point(u)
        S=np.diag([2.,T])
        f,_,_,_,M=integrate(p,T,2*g[0],T*g[1],return_secant=True)
        B=self.R@S@M@S@self.R
        k,V=np.linalg.eigh(B)
        if min(k)<=0. or dt<=0.:
            raise ArithmeticError('Positive finite-bath step mobility required')
        H=self.C@V
        z=H@((2./dt*np.tanh(dt*k/2.))*(H.T@g))
        s=g@z
        if not 0.<s<=f[2]*(1.+2e-13):
            raise ArithmeticError('Modified finite-bath entropy exceeded the physical rate')
        return z,s

    def state_at_force(self,g):
        l=np.log(self.equilibrium_T)
        z=brentq(lambda x:self.point([g[0],x])[2][1]-g[1],l-.015,l+.015,xtol=3e-15)
        return np.array([g[0],z])


def sample(q,eos,T,n):
    b=Filtered(T,[0.,np.log(T)],n,eos)
    z,V=np.linalg.eigh(b.matrices(q))
    base=b.R@V
    base=base/np.linalg.norm(base,axis=0)
    L=np.linalg.inv(base)
    rows=[]
    for k in (.25,4.,20.,80.):
        dt=k/z[-1]
        for mode in (0,1):
            h=2.5e-4
            results=[]
            for sign in (-1.,1.):
                x=b.state_at_force(sign*h*base[:,mode])
                u,e,s=b.step(q,x,dt)
                results.append((b.point(u)[2],e))
            v=L@(results[1][0]-results[0][0])/(2*h)
            target=np.exp(-dt*z[mode])
            if abs(v[mode]-target)>2e-6 or abs(v[1-mode])>2e-6:
                raise ArithmeticError('Closed nonlinear map did not match its exact linear modal load')
            rows.append(dict(fast_dt_times_rate=k,initial_mode=mode,measured_factors=v.tolist(),exact_factor=target,
                             maximum_discrete_entropy_defect=max(r[1] for r in results)))
    return dict(T_SM=T,eigenvalues=z.tolist(),checks=rows)


def run(path):
    eos=StandardModel(path)
    p,q=[PositivePortal(o) for o in ((13,36,26),(17,44,32))]
    rows=[]
    for T in (.3,.6882487043):
        a,b=[sample(j,eos,T,n) for j,n in ((p,112),(q,168))]
        err=max(np.max(abs(np.asarray(x['measured_factors'])-y['measured_factors']))
                for x,y in zip(a['checks'],b['checks'],strict=True))
        if err>2e-6 or np.max(abs(np.asarray(a['eigenvalues'])/b['eigenvalues']-1.))>2e-6:
            raise ArithmeticError('Independent mesh comparison failed')
        rows.append(dict(max_mesh_amplitude_difference=err,**b))
    return dict(rows=rows,table_input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                scope='Near-equilibrium numerical damping in a fixed-energy tabulated FD bath. Modified discrete entropy, not the exact finite-step physical rate.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('eos',type=Path)
    print(json.dumps(run(p.parse_args().eos),indent=2,allow_nan=False))
