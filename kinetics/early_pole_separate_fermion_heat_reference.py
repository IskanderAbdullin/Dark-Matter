"""Coupled pair/decay reference with independent fermion temperatures and finite heat exchange."""
import argparse
import gzip
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, root
from scipy.special import expit, roots_legendre

from early_pole_joint_photon_decay_reference import Coupled
from early_pole_shared_photon_number_budget_reference import photon_fit


class Separate(Coupled):
    def __init__(self,s,g,theta,order,ex,nbudget):
        super().__init__(s,g,theta,order,nbudget)
        self.ex=ex
        d=self.d
        b,p=[np.asarray(s[k]) for k in ('endpoint_fermion_occupations','endpoint_parent_occupations')]
        B=(d.d.wf*b)@d.daughter_energies(self.x0)
        P=(d.pw*p)@d.parent_energies(self.x0)
        assert abs(B+P-self.y0[-2])<2e-15
        self.y0=np.r_[self.y0[:-2],B,P,self.y0[-1]]

    @staticmethod
    def fermion(N,R,E,w):
        X=np.column_stack((np.ones_like(E),E))
        B=(w[:,None]*X).T
        target=np.array([N,R]); scale=np.maximum(abs(target),1e-20)
        def f(z):return expit(X@z)
        def df(z):return X*(f(z)*expit(-X@z))[:,None]
        init=[np.log(max(N,1e-20)/.365),-1.]
        q=root(lambda z:(B@f(z)-target)/scale,init,
               jac=lambda z:(B@df(z))/scale[:,None],tol=2e-11)
        assert q.x[1]<0 and max(abs((B@f(q.x)-target)/scale))<3e-10,(q.message,N,R,q.x)
        return f(q.x),q.x[0],-1/q.x[1]

    def state(self,x,m,y):
        d=self.d
        f=y[:self.ns];NB,NP,RB,RP,RG=y[self.ns:]
        E=np.hypot(m,d.d.p); q=self.weights*f
        s=q@E-(m*m-(1e-9*x)**2)**2/(2*d.d.lam)
        eb=d.daughter_energies(x); ep=d.parent_energies(x)
        b,ab,tb=self.fermion(NB,RB,eb,d.d.wf)
        p,ap,tp=self.fermion(NP,RP,ep,d.pw)
        Tg=(RG/self.gamma_const)**.25
        ag=cap=0.
        if self.number_budget:
            Tg,ag,cap=photon_fit(self.N_gauge-NB-NP,RG)
        Es,Ep,i,j,k,l,u,v,W,_=d.channels(x,m)
        occ=np.r_[f,b,p]
        F,B=d.parent_products(occ[k],occ[l],u,v,False)
        block=occ[j]
        slope=W*(F*(1-block)-B*block)
        R=W*F*(1-block)+slope*f[i]
        E_j=eb[j-self.ns]
        deriv=np.bincount(i,weights=slope,minlength=self.ns)
        deriv_d=np.bincount(i,weights=slope*E_j,minlength=self.ns)
        V=np.bincount(i,weights=R,minlength=self.ns)
        D=R.sum(); Q=R@E_j
        traceP=x*x*np.sum(d.pqw/np.hypot(x,d.pq),axis=1)@p
        mb=d.daughter_mass(x)
        traceB=mb*mb*np.sum(d.daughter_w/np.hypot(mb,d.daughter_q),axis=1)@b
        traceS=(1e-9*x)**2*(q@(1/E))
        Href=np.sqrt(np.pi**2*106.75/90)*(86./x)/2.435e18
        if self.bath:
            nb,hb=self.bath[0].rate(Tg,tb,ab,ag)*self.g**4/Href
            np_,hp=self.bath[1].rate(Tg,tp,ap,ag)*self.g**4/Href
        else:
            nb=hb=np_=hp=0.
        exchange=0.
        if self.ex:
            ff,_,_=self.target(x,np.array([NB,NP,RB+RP]))
            exchange=self.ex*((ff[:d.nf]*d.d.wf)@eb-RB)
            assert exchange*(tp-tb)>-1e-11
        solar=self.visible+s+RB+RP+RG
        bypass=np.sqrt(self.visible/solar)
        v=np.r_[V/self.weights,D+nb,-D+np_,Q+hb+traceB+exchange,
                -Q-V@E+hp+traceP-exchange,-hb-hp]
        entropy=(1/tb-1/Tg)*hb+(1/tp-1/Tg)*hp-(ab-ag)*nb-(ap-ag)*np_+(1/tb-1/tp)*exchange
        assert entropy>=-1e-11
        return (v,bypass,deriv,deriv_d,E,solar,
                np.array([D,nb,np_,traceS+traceP+traceB,exchange]),(Tg,ag,cap,tb,tp,ab,ap,entropy))

    def local(self,t,y):
        x=self.x0*np.exp(t)
        m=self.d.d.gap(x,y[:self.ns])
        v,b,*_=self.state(x,m,y)
        return v*b

    def jac(self,t,y):
        x=self.x0*np.exp(t)
        f=y[:self.ns]
        m=self.d.d.gap(x,f)
        v,b,d,qd,E,S,_,_=self.state(x,m,y)
        J=np.zeros((len(y),len(y)))
        J[:self.ns,:self.ns]=b*np.diag(d/self.weights)
        J[self.ns,:self.ns]=b*d
        J[self.ns+1,:self.ns]=-b*d
        J[self.ns+2,:self.ns]=b*qd
        J[self.ns+3,:self.ns]=-b*(qd+d*E)
        J[:,:self.ns]-=np.outer(v*b,self.weights*E)/(2*S)
        diff=max(2e-12,m*2e-5)
        a=.5*self.d.d.lam
        dm=a*self.weights/E/(2*m+a*((self.weights*np.maximum(f,0.))@(m/E**3)))
        def e(mm,q):
            v,h,*_=self.state(x,mm,q);return v*h
        J[:,:self.ns]+=np.outer((e(m+diff,y)-e(m-diff,y))/(2*diff),dm)
        for k in range(self.ns,len(y)):
            h=max(abs(y[k])*5e-6,2e-11)
            a=y.copy();c=y.copy();a[k]+=h;c[k]-=h
            J[:,k]=(e(m,a)-e(m,c))/(2*h)
        return J

    def energy(self,x,u):
        m=self.d.d.gap(x,u[:self.ns])
        n=self.weights*u[:self.ns]
        rho=n@np.hypot(m,self.d.d.p)-(m*m-(1e-9*x)**2)**2/(2*self.d.d.lam)
        return rho+sum(u[-3:]),rho,m

    def gate(self,x,m,conf,probe=4.1e-6):
        from early_pole_decay_mode_clock_reference import intervals
        Tg,ag,n0,tb,tp,ab,ap,entr=conf
        E,lo,hi,span,pp=intervals(x,m,self.d.daughter_mass(x),probe)
        q,w=roots_legendre(12)
        p=lo+span*(q+1)/2
        fit=expit(ap-p/tp)-expit(ab-(p-E)/tb)
        return np.asarray([.5*w@fit,min(fit),max(fit)])

    def go(self,end):
        span=np.log(end/self.x0); y0=self.y0
        tol=np.r_[np.full(self.ns,3e-8),np.full(5,6e-13)]
        sol=solve_ivp(self.local,(0,span),y0,method='BDF',jac=self.jac,atol=tol,
                      rtol=2e-9,max_step=.003,dense_output=True)
        assert sol.success,sol.message
        n0=self.weights@y0[:self.ns]; y=sol.y[:,-1]
        n1=self.weights@y[:self.ns]
        ini,e0,m0=self.energy(self.x0,y0)
        fin,e1,m1=self.energy(end,y)
        def witness(t):
            x=self.x0*np.exp(t)
            u=sol.sol(t)
            m=self.d.d.gap(x,u[:self.ns])
            _,b,*tail=self.state(x,m,u)
            return b*tail[-2]
        computed=[]
        for order in (4,8):
            p,w=roots_legendre(order)
            h=np.diff(sol.t)/2
            z=(sol.t[:-1,None]+h[:,None]*(p+1)).ravel()
            v=(h[:,None]*w).ravel()
            computed.append(v@np.array([witness(t) for t in z]))
        D,NB,NP,trace,he=computed[-1]
        checks=np.array([(n1-n0-D)/n0,(y[-5]-y0[-5]-D-NB)/n0,
                         (y[-4]-y0[-4]+D-NP)/n0,(fin-ini-trace)/e1])
        assert max(abs(checks))<4e-6,checks
        assert max(abs(computed[-1]-computed[0]))/n0<2e-6
        last=self.state(end,m1,y)[-1]
        first=self.state(self.x0,m0,y0)[-1]
        probes=(1e-6,4.1e-6,1e-5)
        gate_time={}
        def avg(t,q):
            x=self.x0*np.exp(t)
            z=sol.sol(t)
            m=self.d.d.gap(x,z[:self.ns])
            NB,NP,RB,RP,RG=z[self.ns:]
            d=self.d
            _,ab,tb=self.fermion(NB,RB,d.daughter_energies(x),d.d.wf)
            _,ap,tp=self.fermion(NP,RP,d.parent_energies(x),d.pw)
            conf=(None,None,None,tb,tp,ab,ap,None)
            return self.gate(x,m,conf,q)[0]
        times=np.linspace(0,span,100)
        for probe in probes:
            g=np.array([avg(t,probe) for t in times])
            ii=np.flatnonzero(g<=0.)
            xc=None
            if len(ii):
                i=ii[0]
                xc=self.x0*np.exp(brentq(lambda t:avg(t,probe),times[i-1],times[i],xtol=1e-13)) if i else self.x0
            gate_time[str(probe)]=xc
        return dict(g=self.g,mixing_strength=self.ex,photon_number_budget=self.number_budget,x0=self.x0,xend=end,
                    scalar_count_initial_over_T3=n0,scalar_count_final_over_T3=n1,
                    final_daughter_count_over_T3=y[-5],final_parent_count_over_T3=y[-4],
                    final_photon_T=last[0],final_photon_chemical_alpha=last[1],
                    final_daughter_T=last[3],final_parent_T=last[4],
                    initial_reference_fermi_window_difference=self.gate(self.x0,m0,first).tolist(),
                    final_reference_fermi_window_difference=self.gate(end,m1,last).tolist(),
                    first_sampled_negative_average_epoch_by_probe=gate_time,
                    reference_initial_temperatures=[first[3],first[4]],
                    integrated_pair_particle_numbers=[NB,NP],integrated_optional_exchange_heat=he,
                    relative_ledger_errors=checks.tolist())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('spec')
    p.add_argument('--g',type=float,default=.005)
    p.add_argument('--theta',type=float,default=1.)
    p.add_argument('--exchange',type=float,default=100.)
    p.add_argument('--order',type=int,default=10)
    p.add_argument('--xend',type=float,default=.052)
    p.add_argument('--photon-zero-mu',action='store_true')
    a=p.parse_args()
    f=Path(a.spec)
    raw=gzip.decompress(f.read_bytes()) if f.suffix=='.gz' else f.read_bytes()
    d=json.loads(raw)
    s=Separate(d,a.g,a.theta,a.order,a.exchange,not a.photon_zero_mu)
    print(json.dumps(s.go(a.xend),indent=2,allow_nan=False))
