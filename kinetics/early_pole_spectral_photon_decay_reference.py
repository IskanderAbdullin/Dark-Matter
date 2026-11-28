"""Conservative short expanding decay with resolved hard photon populations."""
import argparse
import gzip
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, root
from scipy.special import expit, roots_legendre

from early_pole_joint_photon_decay_reference import Coupled
from early_pole_shared_photon_number_budget_reference import coefficients, photon_fit
from early_pole_spectral_photon_pair_box_reference import SpectralBox


class Live:
    def __init__(self,s,theta,g,order,mix=0.):
        self.base=Coupled(s,0.,theta,order)
        self.box=SpectralBox(s,theta,g,order)
        self.g=g
        self.mix=mix
        self.x0=s['x_end']
        self.ns=self.base.ns
        self.weights=self.base.weights
        self.gamma_weight=self.box.weight
        self.gamma_energy=self.box.weight*self.box.node
        self.y0=np.r_[self.base.y0[:-1],self.box.initial[3:]]
        self.n0_gauge=self.box.N_total

    def original(self,z):
        return np.r_[z[:self.ns+3],self.box.E_static+self.gamma_energy@z[self.ns+3:]]

    def relax(self,f):
        p=self.box.node
        B=np.array([self.gamma_weight,self.gamma_energy])
        moments=B@f
        assert min(moments)>0
        z=p-p.min()
        # Positive inverse temperature and positive distance to the pole
        # of the lowest retained Bose node.
        def dist(u):
            q,e=np.exp(u)
            return q+e*z,np.column_stack((q*np.ones_like(z),e*z))
        def f0(u):
            a,_=dist(u)
            return 1/np.expm1(np.minimum(a,700.))
        def J(u):
            f=f0(u)
            return (B@((-f*(1+f))[:,None]*dist(u)[1]))/moments[:,None]
        guess=np.log([max(p.min(),.05),1.])
        sol=root(lambda q:(B@f0(q)-moments)/moments,guess,jac=J,tol=1e-10)
        a=f0(sol.x)
        assert max(abs(B@a/moments-1))<2e-11,(sol.message,moments,sol.x)
        d=self.mix*(a-f)
        assert max(abs(B@d))<3e-10*max(self.mix,1.)
        return d

    def state(self,t,y):
        x=self.x0*np.exp(t)
        m=self.base.d.d.gap(x,y[:self.ns])
        u=self.original(y)
        rates,H,d,en,base,info=self.base.physics(x,m,u)
        tf,Tg,mu,ag,ct=info
        fj=np.maximum(y[self.ns+3:],0.)
        photons=np.zeros_like(fj)
        pN=[]; heat=0.
        coef=self.g**4/(np.sqrt(np.pi**2*106.75/90)*(86./x)/2.435e18)
        for p,(i,j),a in zip(self.box.pairs,self.box.ij,mu,strict=True):
            A=expit(a-p.f1/tf); B=expit(a-p.f2/tf)
            net=(1-A)*(1-B)*fj[i]*fj[j]-A*B*(1+fj[i])*(1+fj[j])
            event=.5*p.w*net*coef
            photons-=np.bincount(i,weights=event,minlength=len(fj))+np.bincount(j,weights=event,minlength=len(fj))
            pN.append(2*event.sum())
            heat+=event@p.E
        new=np.r_[np.zeros(self.ns),pN,heat,photons/self.gamma_weight]
        if self.mix:
            new[self.ns+3:]+=self.relax(fj)
        v=np.r_[rates[:-1],np.zeros_like(photons)]+new
        assert abs(photons.sum()+sum(pN))<3e-12*max(1.,abs(sum(pN)))
        assert abs(photons@self.box.node+heat)<3e-12*max(1.,abs(heat))
        led=np.array([base[0],*pN,base[-1],heat])
        return v*H,led*H,new*H,en,float((mu[1]-mu[0])*tf-np.hypot(4.1e-6,m))

    def rhs(self,t,y):return self.state(t,y)[0]

    def jac(self,t,y):
        ns=self.ns
        x=self.x0*np.exp(t)
        u=self.original(y)
        # Coupled.jac uses a local time coordinate around its input snapshot.
        jb=self.base.jac(t,u)
        v,ledger,new,en,gate=self.state(t,y)
        J=np.zeros((len(y),len(y)))
        J[:ns+3,:ns]=jb[:-1,:ns]
        J[:,:ns]-=np.outer(new,self.weights*en[:-1])/(2*en[-1])
        for i in range(ns,len(y)):
            h=max(abs(y[i])*6e-6,1e-22)
            a=y.copy(); b=y.copy();a[i]+=h;b[i]-=h
            J[:,i]=(self.rhs(t,a)-self.rhs(t,b))/(2*h)
        return J

    def run(self,end):
        s=np.log(end/self.x0)
        x=self.y0
        cut=self.ns+3
        eps=np.r_[np.full(self.ns,3e-8),np.full(3,5e-13),np.maximum(x[cut:]*8e-10,1e-25)]
        sol=solve_ivp(self.rhs,(0.,s),x,method='BDF',jac=self.jac,rtol=2e-9,atol=eps,
                      max_step=.003,dense_output=True)
        assert sol.success,sol.message
        y=sol.y[:,-1]
        assert np.min(y[cut:])>=0.
        energy,part,m=self.base.energy(end,self.original(y))
        E0,_,_=self.base.energy(self.x0,self.original(x))
        t,w=roots_legendre(6)
        h=np.diff(sol.t)*.5
        xx=(sol.t[:-1,None]+h[:,None]*(t+1)).ravel()
        ww=(h[:,None]*w).ravel()
        led=ww@np.array([self.state(z,sol.sol(z))[1] for z in xx])
        count0=self.weights@x[:self.ns]
        count=self.weights@y[:self.ns]
        ng=self.box.N_static+self.gamma_weight@y[cut:]
        nf=y[self.ns:self.ns+2]
        diff=np.array([(energy-E0-led[3])/part,
               ((nf[0]-x[self.ns])-led[0]-led[1])/count0,
               ((nf[1]-x[self.ns+1])+led[0]-led[2])/count0])
        assert max(abs(diff))<7e-6
        assert abs(ng+nf.sum()-self.n0_gauge)<5e-11
        ph=self.original(y)[-1]
        te,ag,cond=photon_fit(ng,ph)
        target=1/np.expm1(self.box.node/te-ag)
        mismatch=self.gamma_weight@abs(target-y[cut:])/ng
        def gate(t):
            u=sol.sol(t);x=self.x0*np.exp(t);m=self.base.d.d.gap(x,u[:self.ns])
            _,mu,tf=self.base.target(x,u[self.ns:self.ns+3])
            return (mu[1]-mu[0])*tf-np.hypot(4.1e-6,m)
        q=np.array([gate(z) for z in sol.t])
        ids=np.flatnonzero(q<0.)
        crossing=None
        if len(ids):
            i=ids[0]
            crossing=self.x0*np.exp(brentq(gate,sol.t[i-1],sol.t[i])) if i else self.x0
        return dict(g=self.g,xend=end,radial_photon_nodes=len(self.gamma_weight),initial_scalar_count=count0,
                    final_scalar_count=count,final_daughter_number=nf[0],final_parent_number=nf[1],
                    final_photon_number=ng,final_photon_energy=ph,photon_BE_shape_absolute_mismatch=mismatch,
                    assuming_hard_photon_momentum_mixing_rate=self.mix,
                    final_reference_stimulus_gate=q[-1],first_resolved_crossing_epoch=crossing,
                    relative_collision_and_energy_ledgers=diff.tolist(),
                    gauge_particle_count_error=ng+nf.sum()-self.n0_gauge,solver_calls=sol.nfev)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('spec')
    p.add_argument('--g',type=float,default=.005)
    p.add_argument('--order',type=int,default=12)
    p.add_argument('--theta',type=float,default=1.)
    p.add_argument('--xend',type=float,default=.052)
    p.add_argument('--photon-mix',type=float,default=0.)
    a=p.parse_args()
    f=Path(a.spec)
    raw=gzip.decompress(f.read_bytes()) if f.suffix=='.gz' else f.read_bytes()
    d=json.loads(raw)
    print(json.dumps(Live(d,a.theta,a.g,a.order,a.photon_mix).run(a.xend),indent=2,allow_nan=False))
