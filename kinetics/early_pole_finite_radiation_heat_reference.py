"""Conservative finite visible radiation coupled by schematic thermal mixing."""
import argparse
import json

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.special import expit, roots_legendre

from early_pole_decay_mode_clock_reference import X
from early_pole_live_parent_energy_reference import state_moments
from early_pole_massive_daughter_energy_stencil_reference import MassiveFermion
from early_pole_resolved_parent_reference import ResolvedParent, S


class FiniteRadiation(MassiveFermion):
    def __init__(self,ns,nf,nparent,ratio,rate,thermostat):
        super().__init__(ns,nf,nparent,ratio,0.,True)
        self.thermal_rate=rate
        self.thermostat=thermostat

    def evaluation(self,t,z,jac=False):
        x=X*np.exp(t)
        y,dR=z[:-1],z[-1]
        E,D,P,m=state_moments(self,x,y)
        bg=self.radiation if self.thermostat else self.radiation+dR
        theta=(bg/self.radiation)**.25
        total=bg+D+P
        scale=np.sqrt(self.radiation/total)
        v=self.at(x,m,y)
        col=np.zeros(len(z))
        row=np.zeros(len(z))
        J=ResolvedParent.jac(self,t,y) if jac else None
        if self.thermal_rate:
            ep=E[self.iP:]
            p=y[self.iP:]
            N=self.pw@p
            u=ep/theta
            alpha=brentq(lambda q:self.pw@expit(q-u)-N,-70.,100.,xtol=2e-14)
            target=expit(alpha-u)
            k=self.thermal_rate
            delta=k*(target-p)
            energy=self.pw*ep
            v[self.iP:]+=delta
            heat=-energy@delta
            if jac:
                b=target*expit(u-alpha)
                norm=self.pw@b
                A=k*(np.outer(b,self.pw)/norm-np.eye(self.np))
                J[self.iP:,self.iP:]+=A
                row[self.iP:-1]=-energy@A
                if not self.thermostat:
                    col[self.iP:-1]=k*b*(ep-(energy@b)/norm)*theta/(theta**2*4*bg)
                    row[-1]=-energy@col[self.iP:-1]
        else:
            heat=0.
        Q=np.r_[v,heat]
        if not jac:
            return scale*Q
        aug=np.zeros((len(z),len(z)))
        aug[:-1,:-1]=J
        aug[:,-1]=col
        aug[-1,:]=row
        grad=np.r_[self.w*E,0. if self.thermostat else 1.]
        return scale*(aug-np.outer(Q,grad)/(2*total))

    def temperature_jacobian_test(self):
        x=.02
        f=np.r_[np.geomspace(2e5,.1,self.ns),np.linspace(.09,.01,self.nf),
                expit(-self.parent_energies(x))*.95,.012]
        t=np.log(x/X)
        J=self.evaluation(t,f,True)
        v=np.zeros_like(f)
        v[-1]=1.
        h=5e-4
        num=(self.evaluation(t,f+h*v)-self.evaluation(t,f-h*v))/(2*h)
        m=self.iP
        err=np.max(abs(num[m:]-J[m:,-1]))/max(np.max(abs(num[m:])),1.)
        assert err<1e-8
        return float(err)

    def run(self,xend):
        y=np.r_[np.zeros(self.iP),expit(-self.parent_energies(X))]
        atol=np.r_[np.full(self.ns,1e-7),np.full(self.nf,8e-12),
                   np.maximum(y[self.iP:]*5e-10,1e-19),1e-11]
        eta=np.log(xend/X)
        sol=solve_ivp(lambda t,q:self.evaluation(t,q),(0.,eta),np.r_[y,0.],
                      jac=lambda t,q:self.evaluation(t,q,True),
                      atol=atol,rtol=2e-8,max_step=.02,method='BDF',dense_output=True)
        assert sol.success,sol.message
        f,heat=sol.y[:-1,-1],sol.y[-1,-1]
        _,D,P,m=state_moments(self,xend,f)
        _,D0,P0,_=state_moments(self,X,y)
        dt=np.diff(sol.t)/2
        q,w=roots_legendre(4)
        t=(sol.t[:-1,None]+dt[:,None]*(q+1)).ravel()
        W=(dt[:,None]*w).ravel()
        trace=W@np.array([self.trace(s,sol.sol(s)[:-1]) for s in t])
        residual=((D-D0)+(P-P0)+heat-trace)/D
        assert abs(residual)<1e-6,residual
        n=self.d.w@f[:self.ns]
        pair=abs(n-(self.d.wf@f[self.ns:self.iP]))
        totalN=abs(n+self.pw@(f[self.iP:]-y[self.iP:]))
        assert max(pair,totalN)<2e-12
        theta=1. if self.thermostat else (1.+heat/self.radiation)**.25
        return dict(thermostat=self.thermostat,thermal_rate_over_H_ref=self.thermal_rate,
                    cells=[self.ns,self.nf,self.np],daughter_mass_ratio=self.ratio,
                    scalar_yield=float(n/S),added_heat_to_radiation_over_Tref4=float(heat),
                    scalar_yield_over_actual_visible_entropy=float(n/(S*theta**3)),
                    radiation_temperature_over_Tref=float(theta),
                    temperature_jacobian_error=self.temperature_jacobian_test(),
                    final_daughter_energy_over_Tref4=D,final_parent_energy_over_Tref4=P,
                    total_energy_residual_over_final_daughter_energy=float(residual),
                    absolute_final_number_invariant_over_Tref3=float(max(pair,totalN)))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--scalar',type=int,default=35)
    p.add_argument('--fermion',type=int,default=40)
    p.add_argument('--parent',type=int,default=40)
    p.add_argument('--ratio',type=float,default=.06)
    p.add_argument('--rate',type=float,default=100.)
    p.add_argument('--thermostat',action='store_true')
    p.add_argument('--xend',type=float,default=.05)
    a=p.parse_args()
    assert 0<=a.ratio<.7 and a.rate>=0 and X<a.xend<=.05
    d=FiniteRadiation(a.scalar,a.fermion,a.parent,a.ratio,a.rate,a.thermostat)
    print(json.dumps(d.run(a.xend),indent=2,allow_nan=False))
