"""Conserve photon number plus both fermion counts in a pure pair-only box."""
import argparse
import gzip
import json
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.special import expit, zeta

from early_pole_finite_shared_photon_chemistry_reference import Model


def coefficients(alpha):
    if alpha==0:
        return 2*zeta(3)/np.pi**2,np.pi**2/15
    x=np.exp(alpha)
    return 2*float(mp.polylog(3,x))/np.pi**2,6*float(mp.polylog(4,x))/np.pi**2


def photon_fit(N,E):
    A,B=coefficients(0.)
    T=(E/B)**.25
    Nmax=A*T**3
    if N>=Nmax:
        return T,0.,N-Nmax
    def find(a):
        c,d=coefficients(a)
        return c*(E/d)**.75-N
    alpha=brentq(find,-90.,0.,xtol=4e-13)
    c,d=coefficients(alpha)
    T=(E/d)**.25
    assert abs(c*T**3/N-1)<2e-10
    return T,alpha,0.


class NumberModel(Model):
    def __init__(self,s,theta,g,order):
        super().__init__(s,theta,g,order,False)
        self.N0=self.y0[:2].sum()+coefficients(0.)[0]*theta**3

    def rates(self,y):
        a,tf=self.targets(y)
        photon_e=self.E0-y[2]
        ph_n=self.N0-y[:2].sum()
        tg,ag,cond=photon_fit(ph_n,photon_e)
        p=[q.rate(tg,tf,mu,ag) for q,mu in zip(self.pairs,a,strict=True)]
        rates=self.g**4/self.H0*np.array([p[0][0],p[1][0],p[0][1]+p[1][1]])
        return rates, (tg,tf,ag,cond,a)

    def stationary(self):
        def at(a):
            C,D=coefficients(a)
            def fermion(T):
                return self.B@expit(a-self.features[:,2]/T)
            T=brentq(lambda t:fermion(t)[2]+D*t**4-self.E0,.004,1.001*(self.E0/D)**.25,xtol=1e-14)
            n=fermion(T)
            return T,C*T**3+n[:2].sum(),n
        t,n,v=at(0.)
        if n<=self.N0:
            return dict(temperature=t,common_mu_over_T=0.,condensate=self.N0-n,daughter_number=v[0],parent_number=v[1])
        q=brentq(lambda z:at(z)[1]-self.N0,-12.,0.,xtol=1e-12)
        t,n,v=at(q)
        assert abs(n/self.N0-1)<1e-10
        return dict(temperature=t,common_mu_over_T=q,condensate=0.,daughter_number=v[0],parent_number=v[1])

    def go(self,end):
        sol=solve_ivp(lambda t,y:self.rates(y)[0],(0,end),self.y0,t_eval=np.linspace(0,end,25),
                      atol=3e-11,rtol=3e-8,method='DOP853')
        assert sol.success
        rows=[]
        for t,y in zip(sol.t,sol.y.T):
            dr,(tg,tf,ag,cond,a)=self.rates(y)
            entropy=(1/tf-1/tg)*dr[-1]-(a-ag)@dr[:2]
            assert entropy>=-2e-13
            phN=self.N0-y[:2].sum()
            phE=self.E0-y[2]
            c,d=coefficients(ag)
            assert abs(d*tg**4/phE-1)<1e-11
            assert abs((c*tg**3+cond)/phN-1)<2e-10
            rows.append(dict(reference_time=t,photon_T=tg,photon_alpha=ag,photon_condensate=cond,
                             photon_number=phN,photon_energy=phE,fermion_T=tf,
                             fermion_energy=y[-1],daughter_number=y[0],parent_number=y[1],
                             fermion_mu_difference=(a[1]-a[0])*tf,entropy_production=entropy))
        return dict(gauge_coupling=self.g,initial_photon_temperature=self.theta0,
                    conserved_total_particle_number=self.N0,conserved_total_energy=self.E0,
                    analytic_stationary_state=self.stationary(),rows=rows)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('spec')
    p.add_argument('--theta',type=float,default=1.)
    p.add_argument('--end',type=float,default=.8)
    p.add_argument('--order',type=int,default=14)
    a=p.parse_args()
    f=Path(a.spec)
    raw=gzip.decompress(f.read_bytes()) if f.suffix=='.gz' else f.read_bytes()
    d=NumberModel(json.loads(raw),a.theta,.001,(a.order,10,6,6))
    print(json.dumps(d.go(a.end),indent=2,allow_nan=False))
