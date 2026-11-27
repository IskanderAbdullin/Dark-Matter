"""Resolved hard photon occupations against two thermal fermion moment baths."""
import argparse
import gzip
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import expit, roots_legendre

from early_pole_finite_shared_photon_chemistry_reference import Model
from early_pole_shared_photon_number_budget_reference import coefficients, photon_fit


class SpectralBox(Model):
    def __init__(self,s,theta,g,nr):
        super().__init__(s,theta,g,(nr,10,6,6),False)
        z,w=roots_legendre(nr)
        self.node=.5+19.75*(z+1)
        self.weight=19.75*w*self.node**2/(np.pi*np.pi)
        self.ij=[]
        for pair in self.pairs:
            i=np.searchsorted(self.node,pair.g1)
            j=np.searchsorted(self.node,pair.g2)
            assert np.array_equal(self.node[i],pair.g1) and np.array_equal(self.node[j],pair.g2)
            self.ij.append((i,j))
        self.initial=np.r_[self.y0,1/np.expm1(self.node/theta)]
        a,b=coefficients(0.)
        self.N_static=a*theta**3-self.weight@self.initial[3:]
        self.E_static=b*theta**4-(self.weight*self.node)@self.initial[3:]
        assert min(self.N_static,self.E_static)>0
        self.N_total=self.initial[:2].sum()+self.N_static+self.weight@self.initial[3:]
        assert abs((self.initial[2]+self.E_static+(self.weight*self.node)@self.initial[3:])/self.E0-1)<2e-15

    def flux(self,y):
        a,tf=self.targets(y[:3])
        f=np.maximum(y[3:],0.)
        dgamma=np.zeros_like(f); nb=[]
        rateheat=0.
        for A,(i,j),alpha in zip(self.pairs,self.ij,a,strict=True):
            ga,gb=f[i],f[j]
            da=expit(alpha-A.f1/tf);db=expit(alpha-A.f2/tf)
            net=(1-da)*(1-db)*ga*gb-da*db*(1+ga)*(1+gb)
            event=.5*A.w*net*self.g**4/self.H0
            nb.append(2*event.sum())
            rateheat+=event@A.E
            dg=-np.bincount(i,weights=event,minlength=len(f))-np.bincount(j,weights=event,minlength=len(f))
            dgamma+=dg
        assert abs(rateheat+dgamma@self.node)<3e-12
        assert abs(sum(nb)+dgamma.sum())<3e-12
        return np.r_[nb,rateheat,dgamma/self.weight]

    def go(self,end):
        t=np.linspace(0,end,25)
        y=self.initial
        # Strict modal error control also applies in the dilute high-p tail.
        atol=np.r_[np.full(3,2e-12),np.maximum(1e-25,2e-10*y[3:])]
        sol=solve_ivp(lambda t,u:self.flux(u),(0,end),y,method='DOP853',t_eval=t,
                      rtol=2e-8,atol=atol,max_step=.025)
        assert sol.success and min(sol.y.ravel())>=-1e-20
        state=[]
        for x,u in zip(sol.t,sol.y.T,strict=True):
            E=self.E_static+(self.weight*self.node)@u[3:]
            N=self.N_static+self.weight@u[3:]
            errE=E+u[2]-self.E0;errN=N+u[:2].sum()-self.N_total
            assert max(abs(errE),abs(errN))<3e-14
            t0,a,c=photon_fit(N,E)
            fit=1/np.expm1(self.node/t0-a)
            f=np.maximum(u[3:],0.)
            violation=(self.weight@abs(f-fit))/N
            mu,tf=self.targets(u[:3])
            state.append(dict(t=x,daughter_count=u[0],parent_count=u[1],fermion_energy=u[2],
                              photon_number=N,photon_energy=E,fit_T=t0,fit_alpha=a,fit_excess_count=c,
                              not_BE_absolute_shape_mismatch_over_photon_count=violation,
                              fermion_T=tf,fermion_mu_difference=(mu[1]-mu[0])*tf,
                              energy_error=errE,count_error=errN))
        teq=.9
        ytest=np.r_[self.B@expit(-self.features[:,2]/teq),1/np.expm1(self.node/teq)]
        t=self.flux(ytest)
        assert max(abs(t[:3]))<1e-9
        return dict(radial_modes=len(self.node),g=self.g,x=self.x,initial_photon_T=self.theta0,
                    photon_node_over_T=self.node.tolist(),photon_initial_occupation=self.initial[3:].tolist(),
                    photon_final_occupation=sol.y[3:,-1].tolist(),fixed_soft_and_far_tail_number=self.N_static,
                    fixed_soft_and_far_tail_energy=self.E_static,initial_rhs=self.flux(y)[:3].tolist(),rows=state)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('spec')
    p.add_argument('--theta',type=float,default=1.)
    p.add_argument('--g',type=float,default=.001)
    p.add_argument('--order',type=int,default=14)
    p.add_argument('--end',type=float,default=.8)
    a=p.parse_args()
    f=Path(a.spec)
    raw=gzip.decompress(f.read_bytes()) if f.suffix=='.gz' else f.read_bytes()
    s=json.loads(raw)
    d=SpectralBox(s,a.theta,a.g,a.order)
    print(json.dumps(d.go(a.end),indent=2,allow_nan=False))
