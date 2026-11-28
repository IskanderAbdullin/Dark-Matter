"""Coupled scalar production and shared photon chemistry in a fast-FD comparator."""
import argparse
import gzip
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, root
from scipy.special import expit, roots_legendre

from early_pole_finite_shared_photon_chemistry_reference import Pair
from early_pole_joint_fermion_thermal_reference import Joint
from early_pole_shared_photon_number_budget_reference import coefficients, photon_fit


class Coupled:
    def __init__(self,s,g,theta,order,number_budget=False):
        self.spec=s
        self.x0=s['x_end']; self.g=g
        self.number_budget=number_budget
        self.gamma_const=np.pi**2/15
        self.visible=np.pi**2*106.75/30
        self.ratio=s['mass_ratio']
        self.shapes={(67,112,122):(65,90,72),(97,156,176):(95,125,108)}[tuple(s['cells'])]
        self.d=Joint(self.shapes,self.ratio,0.,'none')
        self.bath=[Pair(m,order,8,4,4) for m in (self.ratio*self.x0,self.x0)] if g else []
        self.ns=self.d.ns
        self.weights=self.d.w[:self.ns]
        self.w=np.r_[self.d.d.wf,self.d.pw]
        f=np.array(s['endpoint_scalar_occupations'])
        b=np.array(s['endpoint_fermion_occupations']);p=np.array(s['endpoint_parent_occupations'])
        nB=self.d.d.wf@b;nP=self.d.pw@p
        heat=(self.d.d.wf*b)@self.d.daughter_energies(self.x0)+(self.d.pw*p)@self.d.parent_energies(self.x0)
        self.y0=np.r_[f,nB,nP,heat,self.gamma_const*theta**4]
        self.N_gauge=nB+nP+coefficients(0.)[0]*theta**3

    def target(self,x,mom):
        d=self.d
        Ef=np.r_[d.daughter_energies(x),d.parent_energies(x)]
        b=np.arange(len(Ef))<d.nf
        U=np.column_stack((b,~b,Ef))
        B=(self.w[:,None]*U).T
        sc=np.maximum(abs(mom),1e-20)
        init=np.array([np.log(max(mom[0],1e-15)/.365),np.log(max(mom[1],1e-15)/.365),-1.])
        def f(z):return expit(U@z)
        def df(z):return U*(f(z)*expit(-U@z))[:,None]
        ans=root(lambda v:(B@f(v)-mom)/sc,init,jac=lambda v:(B@df(v))/sc[:,None],tol=2e-11)
        z=ans.x
        new=f(z)
        assert max(abs((B@new-mom)/sc))<1e-10,(x,mom,z)
        assert z[2]<0
        return new,z[:2],-1/z[2]

    def physics(self,x,m,u):
        d=self.d
        f=u[:self.ns]; NB,NP,RF,RG=u[self.ns:]
        n=self.weights*f
        E=np.hypot(m,d.d.p)
        tr_s=(1e-9*x)**2*(n@(1/E))
        R=n@E-(m*m-(1e-9*x)**2)**2/(2*d.d.lam)
        fb,mu,tf=self.target(x,np.array([NB,NP,RF]))
        b=fb[:d.nf];p=fb[d.nf:]
        Tg=(max(RG,1e-30)/self.gamma_const)**.25
        ag=condensate=0.
        if self.number_budget:
            Tg,ag,condensate=photon_fit(self.N_gauge-NB-NP,RG)
        parent_trace=x*x*np.sum(d.pqw/np.hypot(x,d.pq),axis=1)@p
        mm=d.daughter_mass(x)
        daughter_trace=mm*mm*np.sum(d.daughter_w/np.hypot(mm,d.daughter_q),axis=1)@b
        Es,Ep,i,j,k,l,a,c,W,_=d.channels(x,m)
        yy=np.r_[f,fb]
        filled,empty=d.parent_products(yy[k],yy[l],a,c,False)
        B=yy[j]
        loss=W*(filled*(1-B)-empty*B)
        V=np.bincount(i,weights=W*filled*(1-B)+loss*f[i],minlength=self.ns)
        deriv=np.bincount(i,weights=loss,minlength=self.ns)/self.weights
        D=V.sum()
        Q=V@E
        Href=np.sqrt(np.pi**2*106.75/90)*(86./x)/2.435e18
        if self.bath:
            gamma=[self.g**4/Href*pair.rate(Tg,tf,mi,ag) for pair,mi in zip(self.bath,mu,strict=True)]
            srcB,srcP=gamma[0][0],gamma[1][0]
            heat=gamma[0][1]+gamma[1][1]
        else:
            srcB=srcP=heat=0.
        solar=self.visible+R+RF+RG
        Hinv=np.sqrt(self.visible/solar)
        value=np.r_[V/self.weights,D+srcB,-D+srcP,-Q+heat+parent_trace+daughter_trace,-heat]
        return value,Hinv,deriv,np.r_[E,solar],np.array([D,srcB,srcP,heat,tr_s+parent_trace+daughter_trace]),(tf,Tg,mu,ag,condensate)

    def local(self,t,u):
        x=self.x0*np.exp(t)
        m=self.d.d.gap(x,u[:self.ns])
        base,h,*_=self.physics(x,m,u)
        return base*h

    def jac(self,t,u):
        x=self.x0*np.exp(t)
        y=u[:self.ns]
        m=self.d.d.gap(x,y)
        s,c,der,pole,_rv,_=self.physics(x,m,u)
        E=pole[:-1]
        J=np.zeros((len(u),len(u)))
        J[:self.ns,:self.ns]=c*np.diag(der)
        V=self.weights*der
        J[self.ns,:self.ns]=c*V
        J[self.ns+1,:self.ns]=-c*V
        J[self.ns+2,:self.ns]=-c*V*E
        J[:,:self.ns]-=np.outer(c*s,self.weights*E)/(2*pole[-1])
        h=max(2e-12,m*2e-5)
        p1,c1,*_=self.physics(x,m+h,u)
        p0,c0,*_=self.physics(x,m-h,u)
        A=.5*self.d.d.lam
        dN=A*self.weights/E/(2*m+A*((self.weights*np.maximum(y,0.))@(m/E**3)))
        J[:,:self.ns]+=np.outer((p1*c1-p0*c0)/(2*h),dN)
        for i in range(self.ns,len(u)):
            h=max(abs(u[i])*5e-6,2e-11)
            vp=u.copy(); vm=u.copy();vp[i]+=h;vm[i]-=h
            p1,c1,*_=self.physics(x,m,vp)
            p0,c0,*_=self.physics(x,m,vm)
            J[:,i]=(p1*c1-p0*c0)/(2*h)
        return J

    def energy(self,x,u):
        f=u[:self.ns]
        m=self.d.d.gap(x,f)
        bulk=(self.weights*f)@np.hypot(m,self.d.d.p)-(m*m-(1e-9*x)**2)**2/(2*self.d.d.lam)
        return bulk+u[-2]+u[-1],bulk,m

    def go(self,end,rtol=2e-9):
        eta=np.log(end/self.x0)
        y0=self.y0
        eps=np.r_[np.full(self.ns,3e-8),np.full(4,8e-13)]
        sol=solve_ivp(self.local,(0.,eta),y0,method='BDF',jac=self.jac,atol=eps,rtol=rtol,
                      max_step=.004,dense_output=True)
        assert sol.success,sol.message
        n0=self.weights@y0[:self.ns]
        q=sol.y[:,-1]
        n1=self.weights@q[:self.ns]
        whole,part,m=self.energy(end,q)
        whole0,_,m0=self.energy(self.x0,y0)
        def witness(t):
            x=self.x0*np.exp(t)
            u=sol.sol(t)
            m=self.d.d.gap(x,u[:self.ns])
            v,H,d,P,F,_=self.physics(x,m,u)
            return F*H
        integ=[]
        for n in (4,8):
            w,ww=roots_legendre(n)
            h=np.diff(sol.t)/2
            ti=(sol.t[:-1,None]+h[:,None]*(w+1)).ravel()
            wt=(h[:,None]*ww).ravel()
            integ.append(wt@np.array([witness(t) for t in ti]))
        quadtest=np.max(abs(integ[1]-integ[0]))/n0
        assert quadtest<2e-6,quadtest
        integ=integ[-1]
        dS,dB,dP,heat,trace=integ
        errors=np.array([(n1-n0-dS)/n0,(q[-4]-y0[-4]-dS-dB)/n0,
                         (q[-3]-y0[-3]+dS-dP)/n0,(whole-whole0-trace)/part])
        assert np.max(abs(errors))<3e-6,errors
        _,_,_,_,_,(tf,Tg,mu,ag,cond)=self.physics(end,m,q)
        def gate(t):
            x=self.x0*np.exp(t)
            u=sol.sol(t)
            h=self.d.d.gap(x,u[:self.ns])
            _,mu,T=self.target(x,u[self.ns:self.ns+3])
            return (mu[1]-mu[0])*T-np.hypot(4.1e-6,h)
        sampled=np.array([gate(t) for t in sol.t])
        ids=np.flatnonzero(sampled<=0.)
        crossing=None
        if len(ids):
            i=int(ids[0])
            crossing=self.x0*np.exp(brentq(gate,sol.t[i-1],sol.t[i],xtol=1e-13)) if i else self.x0
        return dict(x0=self.x0,xend=end,g=self.g,
                    first_resolved_negative_stimulus_epoch_for_mode_q4p1e6=crossing,
                    final_frozen_stimulus_gate_gap_over_Tref=float(sampled[-1]),
                    scalar_count_initial_over_T3=n0,scalar_count_final_over_T3=n1,
                    final_parent_number_over_T3=q[-3],final_daughter_number_over_T3=q[-4],
                    final_fermion_T_over_Tref=tf,final_photon_T_over_Tref=Tg,
                    photon_number_budget=self.number_budget,final_photon_mu_over_Tgamma=ag,
                    final_photon_condensate_count_over_T3=cond,
                    final_mu_difference_over_Tref=(mu[1]-mu[0])*tf,
                    final_scalar_pole_over_T=m,integrated_scalar_decay_count_over_T3=dS,
                    integrated_photon_pair_counts=[dB,dP],integrated_pair_heat_over_T4=heat,
                    integrated_mass_trace_over_T4=trace,relative_ledger_errors=errors.tolist(),
                    relative_quadrature_disagreement=float(quadtest),relative_solver_tolerance=rtol,
                    solver_calls=sol.nfev)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('spec')
    p.add_argument('--g',type=float,default=.001)
    p.add_argument('--theta',type=float,default=1.)
    p.add_argument('--xend',type=float,default=.0505)
    p.add_argument('--order',type=int,default=10)
    p.add_argument('--rtol',type=float,default=2e-9)
    p.add_argument('--photon-number',action='store_true')
    a=p.parse_args()
    f=Path(a.spec)
    raw=gzip.decompress(f.read_bytes()) if f.suffix=='.gz' else f.read_bytes()
    d=json.loads(raw)
    assert d['x_end']<a.xend<=1.05*d['x_end']
    print(json.dumps(Coupled(d,a.g,a.theta,a.order,a.photon_number).go(a.xend,a.rtol),indent=2,allow_nan=False))
