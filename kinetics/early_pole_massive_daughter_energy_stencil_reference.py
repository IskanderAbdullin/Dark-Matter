"""A genuinely massive fermion daughter in the conservative early stencil."""
import argparse
import json

import numpy as np
from scipy.special import expit, roots_legendre

from early_pole_decay_mode_clock_reference import KAPPA, X, kernel
from early_pole_elastic_energy_combined_reference import ConservativeMix


class MassiveFermion(ConservativeMix):
    def __init__(self,ns,nf,nparent,ratio,rate,full_H):
        super().__init__(ns,nf,nparent,rate,full_H)
        self.ratio=ratio
        q,w=roots_legendre(16)
        edge=self.d.fed
        change=np.diff(edge)
        self.daughter_q=edge[:-1,None]+change[:,None]*(q+1)/2
        self.daughter_w=4*change[:,None]*w*self.daughter_q**2/(4*np.pi**2)

    def daughter_mass(self,x):
        return x*self.ratio

    def daughter_energies(self,x):
        if self.ratio==0:
            return super().daughter_energies(x)
        return np.sum(self.daughter_w*np.hypot(self.daughter_mass(x),self.daughter_q),axis=1)/self.d.wf

    def trace(self,t,y):
        tr=super().trace(t,y)
        x=X*np.exp(t)
        b=self.daughter_mass(x)
        if b:
            tr+=b*b*np.sum(self.daughter_w/np.hypot(b,self.daughter_q),axis=1)@y[self.ns:self.iP]
        return tr

    def stencil_checks(self,x=.008,m=8e-5):
        rate=self.rate
        try:
            self.rate=0.  # test bare collisions, not the independent BGK target
            return super().stencil_checks(x,m)
        finally:
            self.rate=rate

    def seed_control(self,x=.02,m=.0005):
        b=self.daughter_mass(x)
        E,P,i,j,k,l,a,h,W,miss=self.channels(x,m)
        parent=expit(-P)
        A,_=self.parent_products(parent[k-self.iP],parent[l-self.iP],a,h,False)
        projected=np.bincount(i,weights=W*A,minlength=self.ns)/self.d.w
        exact=[]
        for p in self.d.p:
            q=kernel(x,m,b,p)
            exact.append(KAPPA*x*x*q['seed_over_width']*q['vacuum_matched_yukawa_width_ratio'])
        exact=np.array(exact)
        active=exact>exact.max()*1e-5
        return dict(max_relative_active_seed_error=float(np.max(abs(projected[active]/exact[active]-1))),
                    signed_weighted_seed_error=float((self.d.w@projected)/(self.d.w@exact)-1),
                    omitted_thermal_seed_fraction=miss,
                    x=x,scalar_pole_over_T=m,daughter_pole_over_T=b,active_rows=int(active.sum()))

    def run(self,xend):
        checks=self.stencil_checks()
        out=super().run(xend)
        out['daughter_vacuum_mass_over_parent']=self.ratio
        out['massive_daughter_discrete_collision_checks']=checks
        out['independent_continuum_seed_check']=self.seed_control()
        return out


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--scalar',type=int,default=35)
    p.add_argument('--fermion',type=int,default=40)
    p.add_argument('--parent',type=int,default=40)
    p.add_argument('--rate',type=float,default=0.)
    p.add_argument('--ratio',type=float,default=.06)
    p.add_argument('--xend',type=float,default=.05)
    p.add_argument('--static-only',action='store_true')
    p.add_argument('--seed-only',action='store_true')
    p.add_argument('--radiation-only-H',action='store_true')
    a=p.parse_args()
    assert 0.<=a.ratio<.7 and a.rate>=0.
    d=MassiveFermion(a.scalar,a.fermion,a.parent,a.ratio,a.rate,not a.radiation_only_H)
    print(json.dumps(d.seed_control() if a.seed_only else d.stencil_checks() if a.static_only else d.run(a.xend),
                     indent=2,allow_nan=False))
