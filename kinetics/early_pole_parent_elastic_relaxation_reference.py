"""Number- and projected-energy-conserving parent momentum relaxation."""

import argparse
import json

import numpy as np
from scipy.optimize import root
from scipy.special import expit

from early_pole_resolved_parent_reference import ResolvedParent


class ElasticBath(ResolvedParent):
    """Illustrative isolated parent mixing, not an elastic matrix element."""
    def __init__(self, ns,nf,npmid,rate):
        super().__init__(ns,nf,npmid,gamma=0.,held=False)
        self.rate=rate

    def elastic_target(self,x,p,jac=False):
        # This quadrature moment defines parent energy consistently for
        # inventory and target; the decays can still have centroid defects.
        e=np.sum(self.pqw*np.hypot(x,self.pq),axis=1)/self.pw
        U=np.column_stack((np.ones_like(e),e))
        B=(self.pw[:,None]*U).T
        if np.all(p==0.) or np.all(p==1.):
            return np.zeros((self.np,self.np)) if jac else p.copy()
        moments=B@p
        def value(z):
            return expit(U@z)
        def deriv(z):
            q=U@z
            return U*(expit(q)*expit(-q))[:,None]
        sol=root(lambda z:B@value(z)-moments,np.array([0.,-1.]),
                 jac=lambda z:B@deriv(z),tol=2e-11)
        f=value(sol.x)
        if np.max(abs((B@f)/moments-1))>2e-12:
            raise RuntimeError(str(sol.message))
        if jac:
            D=deriv(sol.x)
            return D@np.linalg.solve(B@D,B)
        return f

    def at(self,x,m,y,jac=False):
        z=super().at(x,m,y,jac)
        if self.rate:
            if jac:
                z[self.iP:,self.iP:]+=self.rate*(self.elastic_target(x,y[self.iP:],True)-np.eye(self.np))
            else:
                z[self.iP:]+=self.rate*(self.elastic_target(x,y[self.iP:])-y[self.iP:])
        return z

    def check(self):
        # The strict base check isolates the collision operator. The BGK
        # moment error has its own scaled check before it is multiplied by
        # the requested rate; the trajectory still checks total count.
        rate=self.rate
        try:
            self.rate=0.
            counts=super().check()
        finally:
            self.rate=rate
        self.elastic_checks()
        return counts

    def run(self,xend):
        z=super().run(xend)
        z['energy_conserving_parent_mixing_per_log_a']=self.rate
        return z

    def elastic_checks(self):
        x=.02
        p=self.parent_fd(x)*(.53+.33*np.exp(-self.pe[:-1]/2))
        v=np.sin(self.pe[1:])*p
        f=self.elastic_target(x,p)
        J=self.elastic_target(x,p,True)
        eps=1e-5
        numerical=(self.elastic_target(x,p+eps*v)-self.elastic_target(x,p-eps*v))/(2*eps)
        e=np.sum(self.pqw*np.hypot(x,self.pq),axis=1)/self.pw
        B=np.vstack((self.pw,self.pw*e))
        residual=(B@f)/(B@p)-1
        err=np.max(abs(numerical-J@v))/np.max(abs(numerical))
        assert np.max(abs(residual))<3e-13
        assert err<3e-8
        assert np.all((f>=0)&(f<=1))
        return dict(number_relative_defect=float(residual[0]),
                    energy_relative_defect=float(residual[1]),
                    directional_derivative_error=float(err))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--scalar',type=int,default=35)
    p.add_argument('--fermion',type=int,default=40)
    p.add_argument('--parent',type=int,default=40)
    p.add_argument('--xend',type=float,default=.05)
    p.add_argument('--rate',type=float,default=3.)
    p.add_argument('--static-only',action='store_true')
    a=p.parse_args()
    assert a.rate>=0.
    d=ElasticBath(a.scalar,a.fermion,a.parent,a.rate)
    print(json.dumps(d.elastic_checks() if a.static_only else d.run(a.xend),
                     indent=2,allow_nan=False))
