"""Collisionless energy and realizability bounds for unresolved source tails.

One mean lab energy per species/side/injection time. Count and injected
energy are exact, the subsequent redshifted energy is a Jensen lower bound
for arbitrary positive distributions on each tail.
"""
import numpy as np
from scipy.special import roots_legendre


class Exterior:
    def __init__(self):
        self.cohorts = []
        self.a = 1.
        self.trace = 0.
        self.injected = 0.

    def add(self, a, mass, channel, side, boundary_p, n_a3, k_a4):
        """Add positive number and kinetic-energy increments (comoving units)."""
        assert a==self.a and n_a3>=0 and k_a4>=0
        if n_a3==0:
            assert k_a4==0
            return
        K=k_a4/(a*n_a3)
        mean=mass+K
        bound=np.hypot(mass,boundary_p)
        if side==0: assert mean<=bound*(1+2e-12)
        else: assert mean>=bound*(1-2e-12)
        self.cohorts.append((a,mass,channel,side,bound,n_a3,mean))
        self.injected+=k_a4+a*mass*n_a3

    def array(self):
        return np.asarray(self.cohorts,dtype=float).reshape(-1,7).T

    def energy(self, a):
        ai,m,_,_,_,n,e=self.array()
        r=ai/a
        return float(a*np.sum(n*np.hypot(m*np.sqrt(np.maximum(0.,1-r*r)),r*e)))

    def upper_bound(self, a):
        ai,m,_,side,bound,n,e=self.array()
        r=ai/a
        f=lambda x:np.hypot(m*np.sqrt(np.maximum(0.,1-r*r)),r*x)
        # The finite low tail is bounded by a chord. The semi-infinite
        # high tail is bounded by the tangent slope limit r at infinity.
        y=f(bound)
        dy=np.divide(y-m,bound-m,out=np.zeros_like(y),where=bound>m)
        hi=np.where(side==0,m+dy*(e-m),y+r*(e-bound))
        return float(a*np.sum(n*hi))

    def advance(self,a):
        assert a>=self.a
        if a==self.a: return
        ai,m,_,_,_,n,e=self.array()
        l=np.log(a/self.a)
        g,w=roots_legendre(12)
        x=self.a*np.exp(l*(g+1)/2)
        r=ai[:,None]/x
        mass=m[:,None]
        Eval=np.hypot(mass*np.sqrt(np.maximum(0.,1-r*r)),r*e[:,None])
        de=l/2*np.sum((n[:,None]*mass*mass*x/Eval)*w)
        before=self.energy(self.a)
        assert abs(self.energy(a)-before-de)<=3e-12*max(before,self.energy(a),1e-280)
        self.trace+=de
        self.a=a

    def result(self):
        e=self.energy(self.a)
        u=self.upper_bound(self.a)
        defect=abs(e-self.injected-self.trace)/max(self.injected,1e-280)
        assert defect<3e-12
        assert u>=e-3e-13*max(e,1e-280)
        return dict(exterior_cohorts=len(self.cohorts),
             exterior_a4_energy_one_mean_per_cohort=e,
             exterior_a4_energy_distribution_independent_upper_bound=u,
             exterior_a4_pressure_trace=self.trace,
             exterior_a4_portal_injection=self.injected,
             exterior_relative_energy_ledger_defect=defect)


def test():
    t=Exterior()
    rng=np.random.default_rng(11)
    full=[]
    for a in (1.,1.02,1.05):
        t.advance(a)
        for mass in (.02,8.6e-8):
            b=.015
            for side in (0,1):
                p=b*(rng.random(100) if side==0 else 1.+10*rng.random(100))
                w=rng.random(100)*.0001
                k=p*p/(np.hypot(p,mass)+mass)
                t.add(a,mass,0,side,b,w.sum(),a*np.sum(w*k))
                full.append((a,mass,w,p))
    t.advance(1.11)
    r=t.result()
    exact=1.11*sum(w@np.hypot(m,p*a/1.11) for a,m,w,p in full)
    assert r["exterior_a4_energy_one_mean_per_cohort"]<=exact<=r["exterior_a4_energy_distribution_independent_upper_bound"]
    return dict(exact_redshifted_positive_samples=exact,**r)


if __name__=="__main__":
    import json
    print(json.dumps(test(),indent=2,allow_nan=False))
