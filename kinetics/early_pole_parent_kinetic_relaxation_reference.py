"""Separate parent momentum relaxation from chemical number replenishment."""

import argparse
import json

import numpy as np
from scipy.optimize import brentq
from scipy.special import expit

from early_pole_resolved_parent_reference import ResolvedParent


class KineticBath(ResolvedParent):
    """Illustrative fixed-temperature, exactly number-conserving BGK target."""

    def __init__(self, ns, nf, npmid, rate):
        super().__init__(ns,nf,npmid,gamma=0.,held=False)
        self.rate=rate
        self.cap=np.sum(self.pw)

    def thermal_target(self,x,p):
        target=self.pw@p
        if target<0 or target>self.cap:
            raise ValueError("nonrealizable total parent count")
        if target==0:
            return np.zeros(self.np),np.zeros(self.np),0.
        if target==self.cap:
            return np.ones(self.np),np.zeros(self.np),0.
        E=np.hypot(self.pq,x)
        def f(mu):
            q=expit(mu-E)
            n=np.sum(self.pqw*q,axis=1)
            df=np.sum(self.pqw*q*expit(E-mu),axis=1)
            return n,df
        mu=brentq(lambda u:f(u)[0].sum()-target,-80.,110.,xtol=3e-14)
        n,df=f(mu)
        return n/self.pw,df/self.pw,df.sum()

    def at(self,x,m,y,jac=False):
        v=super().at(x,m,y,jac)
        if not self.rate:
            return v
        f,df,χ=self.thermal_target(x,y[self.iP:])
        if jac:
            # Implicit derivative of the fixed-number chemical potential.
            M=np.outer(df,self.pw)/χ if χ else np.zeros((self.np,self.np))
            v[self.iP:,self.iP:]+=self.rate*(M-np.eye(self.np))
        else:
            v[self.iP:]+=self.rate*(f-y[self.iP:])
        return v

    def run(self,xend):
        z=super().run(xend)
        z['kinetic_relaxation_per_log_a']=self.rate
        return z

    def target_checks(self):
        x=.02
        p=self.parent_fd(x)
        rng=np.random.default_rng(55756)
        p*=rng.uniform(.25,.75,self.np)
        target,df,dn=self.thermal_target(x,p)
        assert np.max(target)<1 and np.min(target)>=0
        number=self.pw@p
        assert abs(self.pw@target-number)<number*3e-14
        v=rng.uniform(-1,1,self.np)*p
        hh=3e-5
        delta=(self.thermal_target(x,p+hh*v)[0]-self.thermal_target(x,p-hh*v)[0])/(2*hh)
        analytic=df*(self.pw@v)/dn
        err=np.max(abs(delta-analytic))/np.max(abs(analytic))
        assert err<2e-8
        return dict(target_number_defect=float((self.pw@target)/number-1),
                    directional_derivative_error=float(err))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--scalar',type=int,default=35)
    p.add_argument('--fermion',type=int,default=40)
    p.add_argument('--parent',type=int,default=40)
    p.add_argument('--rate',type=float,default=3.)
    p.add_argument('--xend',type=float,default=.05)
    p.add_argument('--static-only',action='store_true')
    a=p.parse_args()
    assert a.rate>=0.
    d=KineticBath(a.scalar,a.fermion,a.parent,a.rate)
    print(json.dumps(d.target_checks() if a.static_only else d.run(a.xend),
                     indent=2,allow_nan=False))
