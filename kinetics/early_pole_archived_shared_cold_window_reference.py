"""Matched cold scalar support across archives with shared gauge and Hartree energy."""
import argparse
import json
from types import SimpleNamespace

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.special import expit

from early_pole_archived_extended_pair_decay_window_reference import mount
from early_pole_live_pauli_gauge_moment_decay_reference import Tool
from early_pole_decay_mode_clock_reference import G,KAPPA


class Cold:
    def __init__(self,grid,wide,level,n,aperture_refined=False,charge_floor=.45,field_floor=.45,exposure=8e-5,reference_aperture_grid=None,archive_projection="steps",aperture_subdivision=1,shape_amplitude=0.,shape_degree=3):
        self.exposure=exposure
        self.M,d,quad,remap=mount(grid,wide,level,aperture_refined,charge_floor,field_floor,reference_aperture_grid,archive_projection,aperture_subdivision,shape_amplitude,shape_degree)
        self.info=dict(grid=grid,wide=wide,quad=quad,initial_hard_moment_difference=remap,
                       aperture_refined=aperture_refined,charge_floor=charge_floor,field_floor=field_floor,hard_cells=[k.stop-k.start for k in self.M.r])
        if reference_aperture_grid is not None:
            self.info['reference_aperture_grid']=reference_aperture_grid
        if archive_projection!='steps':
            self.info['archive_projection']=archive_projection
        if aperture_subdivision>1:
            self.info['aperture_subdivision']=aperture_subdivision
        if shape_amplitude:
            self.info.update(shape_amplitude=shape_amplitude,shape_degree=shape_degree)
        old=np.asarray(d['scalar_momentum_edges'])
        counts=np.array(d['endpoint_scalar_signed_shell_number_over_T3'])
        assert -np.minimum(counts,0.).sum()<1e-10*counts.sum()
        counts=np.maximum(counts,0.)
        edge=np.geomspace(1e-6,1e-5,n+1)
        cube=lambda e:e**3
        ov=np.maximum(0.,np.minimum(cube(old[1:])[:,None],cube(edge[1:]))-
                           np.maximum(cube(old[:-1])[:,None],cube(edge[:-1])))
        weights=ov/np.diff(cube(old))[:,None]
        self.ns=counts@weights
        self.spect=counts*(1-weights.sum(axis=1))
        self.spq=np.asarray(d['scalar_momentum_nodes'])
        self.q=np.cbrt((cube(edge[1:])+cube(edge[:-1]))/2)
        self.sw=np.diff(cube(edge))/(6*np.pi**2)
        self.n=n
        self.m0=d['checkpoints'][-1]['scalar_pole_over_T']
        pole=np.hypot(self.q,self.m0)
        self.lam=2*self.m0**2/(self.ns@np.reciprocal(pole)+self.spect@np.reciprocal(np.hypot(self.spq,self.m0)))
        self.j0=self.M.init.size
        assert min(self.ns)>0. and abs(self.ns.sum()+self.spect.sum()-counts.sum())<1e-15
        self.lambdas=[SimpleNamespace(n=self.M.decay.n,cells=self.M.decay.cells,
                                      q=q,nod=self.M.z,w=self.M.w) for q in self.q]

    def scalar(self,N):
        def gap(m):
            return m*m-.5*self.lam*(N@np.reciprocal(np.hypot(self.q,m))+
                                   self.spect@np.reciprocal(np.hypot(self.spq,m)))
        m=brentq(gap,1e-8,.03,xtol=1e-17)
        E=np.hypot(self.q,m)
        rho=N@E+self.spect@np.hypot(self.spq,m)-m**4/(2*self.lam)
        f=N/self.sw
        assert min(f)>0.
        S=self.sw@(np.log1p(f)+f*np.log1p(1/f))
        return E,rho,S,f,m

    def rhs(self,y,pair):
        Y=y[:self.j0].reshape(-1,2)
        N=y[self.j0:self.j0+self.n]
        e,f,grad,charged_rho,_=self.M.weights(Y,True)
        E,rho,S,fs,m=self.scalar(N)
        self.last_total_energy=charged_rho+rho
        out=np.zeros_like(Y)
        dN=np.zeros_like(N);S=0.
        x=.05*np.exp(getattr(self.M,"eta",0.))
        coeff=self.sw*G*KAPPA*x*((1.06*x)**2-m*m)/(self.q*E)
        for j,es in enumerate(E):
            model=self.lambdas[j]
            for i,k,u,v,w in Tool.events(model,dict(coeff=e[self.M.start:],es=es)):
                i+=self.M.start;k+=self.M.start
                a=f[i,0]+f[i,1]*u
                b=f[k,0]+f[k,1]*v
                F=expit(a);B=expit(b)
                R=w*(F*(1-B)+fs[j]*(F-B))*coeff[j]
                z=R.sum()
                out[i]-=np.r_[z,R@u]
                out[k]+=np.r_[z,R@v]
                dN[j]+=z
                S+=R@(a-b+np.log1p(1/fs[j]))
        assert S>=-1e-13
        residual=np.sum(grad*out)+E@dN
        assert abs(residual)<2e-13*max(np.sum(abs(grad*out)),1.)
        if pair:
            v=self.M.pair_rhs(Y,True)
            out+=v[:-2].reshape(-1,2)
            S+=v[-2]
            residual+=v[-1]
        return np.r_[out.ravel(),dN,S,residual]

    def one(self,pair):
        z=np.r_[self.M.init.ravel(),self.ns,0.,0.]
        def totals(a):
            _,_,_,r,s=self.M.weights(a[:self.j0].reshape(-1,2),True)
            _,rs,ss,_,pole=self.scalar(a[self.j0:self.j0+self.n])
            return r+rs,s+ss,pole
        e0,s0,_=totals(z)
        Y=solve_ivp(lambda _,x:self.rhs(x,pair),(0.,self.exposure),z,method='DOP853',
                    rtol=2e-10,atol=3e-13,max_step=3e-5,
                    t_eval=np.linspace(0.,self.exposure,13))
        assert Y.success,Y.message
        a=Y.y[:,-1]
        e1,s1,pole=totals(a)
        assert abs(e1-e0)<2e-13 and abs(e1-e0-a[-1])<2e-13
        assert abs(s1-s0-a[-2])<3e-11
        d=np.sum((a[:self.j0]-z[:self.j0]).reshape(-1,2)[:,0])
        assert abs(d)<2e-13
        dn=(a-z)[self.j0:self.j0+self.n]
        tot=np.sum(Y.y[self.j0:self.j0+self.n],axis=0)
        return dict(pair=pair,absolute_scalar_gain=dn.sum(),scalar_gain_by_q=dn.tolist(),
                    net_gain_timeline=(tot-tot[0]).tolist(),
                    scalar_initial_bin_counts=self.ns.tolist(),final_scalar_pole=pole,
                    gauge_mixed_count_drift=d,energy_drift=e1-e0,entropy_gain=s1-s0,
                    entropy_ledger=a[-2])

    def run(self):
        nom,pair=self.one(False),self.one(True)
        return dict(scope='Same [1e-6,1e-5] scalar support, vacuum hard rate proxies in a fixed box; '
                          'no full rate, scale evolution, or late relic.',**self.info,
                    scalar_q=self.q.tolist(),exposure=self.exposure,rows=[nom,pair],remap_rejection=self.M.reject)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--grid',choices=['fine','third'],default='third')
    p.add_argument('--wide',type=int,choices=[0,1,2],default=0)
    p.add_argument('--level',type=int,default=0)
    p.add_argument('--n',type=int,default=4)
    p.add_argument("--refine-aperture",action="store_true")
    p.add_argument("--floor",type=float,default=.45)
    p.add_argument("--field-floor",type=float,default=.45)
    p.add_argument("--exposure",type=float,default=8e-5)
    o=p.parse_args()
    print(json.dumps(Cold(o.grid,o.wide,o.level,o.n,o.refine_aperture,o.floor,o.field_floor,o.exposure).run(),indent=2,allow_nan=False))
