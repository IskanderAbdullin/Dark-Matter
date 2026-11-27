"""Separate conserved visible radiation and early dark budget in H."""
import argparse
import json

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import expit, roots_legendre

from early_pole_decay_mode_clock_reference import X
from early_pole_exact_energy_parent_stencil_reference import EnergyStencil, S
from early_pole_live_parent_energy_reference import state_moments


class FriedmannStencil(EnergyStencil):
    def __init__(self,ns,nf,nparent,gstar,correct):
        super().__init__(ns,nf,nparent)
        self.gstar=gstar
        self.radiation=np.pi**2*gstar/30
        self.correct=correct

    def rhs(self,t,y):
        x=X*np.exp(t)
        _,d,p,m=state_moments(self,x,y)
        scale=np.sqrt(self.radiation/(self.radiation+d+p)) if self.correct else 1.
        return scale*self.at(x,m,y)

    def jac(self,t,y):
        j=super().jac(t,y)
        if not self.correct:
            return j
        x=X*np.exp(t)
        E,d,p,m=state_moments(self,x,y)
        total=self.radiation+d+p
        scale=np.sqrt(self.radiation/total)
        return scale*(j-np.outer(self.at(x,m,y),self.w*E)/(2*total))

    def trace(self,t,y):
        x=X*np.exp(t)
        scalar=y[:self.ns]
        m=self.d.gap(x,scalar)
        # The stationary Hartree stress removes the self-mass contribution
        # from rho-3P; the actual bare vacuum mass supplies the scalar trace.
        trS=(1e-9*x)**2*((self.d.w*scalar)@(1/np.hypot(self.d.p,m)))
        trP=x*x*np.sum(self.pqw/np.hypot(x,self.pq),axis=1)@y[self.iP:]
        return trP+trS

    def run(self,xend):
        z0=np.r_[np.zeros(self.iP),expit(-self.parent_energies(X))]
        atol=np.r_[np.full(self.ns,1e-7),np.full(self.nf,8e-12),
                   np.maximum(z0[self.iP:]*5e-10,1e-19)]
        t1=np.log(xend/X)
        sol=solve_ivp(self.rhs,(0.,t1),z0,method='BDF',jac=self.jac,
                      rtol=2e-8,atol=atol,max_step=.020,dense_output=True)
        assert sol.success,sol.message
        x0,xf=X,xend
        _,di,pi,_=state_moments(self,x0,z0)
        _,df,pf,m=state_moments(self,xf,sol.y[:,-1])
        check=[]
        for order in (3,6):
            nodes,weights=roots_legendre(order)
            widths=np.diff(sol.t)/2
            t=(sol.t[:-1,None]+widths[:,None]*(nodes+1)).ravel()
            w=(widths[:,None]*weights).ravel()
            check.append(w@np.array([self.trace(a,sol.sol(a)) for a in t]))
        defect=((df-di)+(pf-pi)-check[-1])/df
        assert abs(defect)<4e-6,defect
        assert abs(check[0]-check[1])/df<3e-7
        f=sol.y[:,-1]
        n=float(self.d.w@f[:self.ns])
        pair=abs(n-(self.d.wf@f[self.ns:self.iP]))
        parent=abs(n+self.pw@(f[self.iP:]-z0[self.iP:]))
        assert max(pair,parent)<3e-12
        return dict(visible_gstar=self.gstar,self_consistent_H=self.correct,
                    cells=[self.ns,self.nf,self.np],xend=xend,
                    dark_to_radiation_initial=float((di+pi)/self.radiation),
                    dark_to_radiation_final=float((df+pf)/self.radiation),
                    initial_H_over_H_radiation=float(np.sqrt(1+(di+pi)/self.radiation)),
                    final_H_over_H_radiation=float(np.sqrt(1+(df+pf)/self.radiation)),
                    scalar_yield_from_visible_entropy=n/S*(106.75/self.gstar),
                    final_daughter_energy_over_T4=df,initial_parent_energy_over_T4=pi,
                    final_parent_energy_over_T4=pf,final_scalar_pole_over_T=m,
                    integrated_trace_over_T4=float(check[-1]),
                    continuity_residual_over_daughter_energy=float(defect),
                    largest_absolute_endpoint_number_invariant=float(max(pair,parent)))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--scalar',type=int,default=35)
    p.add_argument('--fermion',type=int,default=40)
    p.add_argument('--parent',type=int,default=40)
    p.add_argument('--gstar',type=float,default=106.75)
    p.add_argument('--radiation-only-H',action='store_true')
    p.add_argument('--xend',type=float,default=.05)
    a=p.parse_args()
    assert a.gstar>0 and X<a.xend<=.05
    d=FriedmannStencil(a.scalar,a.fermion,a.parent,a.gstar,not a.radiation_only_H)
    print(json.dumps(d.run(a.xend),indent=2,allow_nan=False))
