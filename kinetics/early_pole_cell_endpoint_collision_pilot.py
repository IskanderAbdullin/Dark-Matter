"""Discontinuous per-cell two-node collision chart without vanishing face volumes."""
import argparse
import json
from math import fsum

import numpy as np
from scipy.special import expit,roots_legendre

from early_pole_decay_mode_clock_reference import G,KAPPA,X
from early_pole_bounded_face_extents_pilot import step
from early_pole_double_charge_stencil_pilot import event
from early_pole_endpoint_stencil_remap_pilot import (
    LABELS,SPECS,archives,load,parts,remap,MassiveFermion)
from early_pole_live_bounded_face_extents_pilot import hold_trace


def split(e,n):
    a=np.concatenate([np.linspace(l,h,n+1)[:-1]
                      for l,h in zip(e[:-1],e[1:],strict=True)]+[e[-1:]])
    if not np.array_equal(a[::n],e):
        raise ArithmeticError('Original step interfaces moved')
    return a


def geometry(e,m,deg):
    l,r=e[:-1],e[1:]
    E=np.hypot(e,m)
    z,w=roots_legendre(20)
    h=(r-l)/2
    p=(r+l)[:,None]/2+h[:,None]*z
    V=deg*(r**3-l**3)/(6*np.pi**2)
    integral=deg*h/(2*np.pi**2)*((p*p*np.hypot(p,m))@w)
    b=(integral/V-E[:-1])/np.diff(E)
    if min(b)<=0 or max(b)>=1:
        raise ArithmeticError('Exact cell energy not inside endpoint span')
    W=np.c_[V*(1-b),V*b].ravel()
    Q=np.c_[l,r].ravel()
    if np.max(abs(((W*np.hypot(Q,m)).reshape(-1,2).sum(1)-integral)/integral))>2e-14:
        raise ArithmeticError('Endpoint quadrature lost its initial energy moment')
    return Q,W


class CellChart:
    parent_products=staticmethod(MassiveFermion.parent_products)

    def __init__(self,arc,shape,n,reference_x=None):
        x=arc['x_end'] if reference_x is None else reference_x
        self.d0=MassiveFermion(*shape,arc['mass_ratio'],arc['elastic_rate'],True)
        if not all(np.array_equal(a,b) for a,b in zip(archives(arc),parts(self.d0),strict=True)):
            raise ArithmeticError('Unexpected archive grid')
        self.ds=self.d0.d
        self.ratio=self.d0.ratio
        self.radiation=self.d0.radiation
        e0,e1,e2=[np.asarray(a) for a in archives(arc)]
        self.d_edges=split(e1,n)
        self.p_edges=split(e2,n)
        self.qd,self.wd=geometry(self.d_edges,x*self.ratio,4)
        self.qp,self.wp=geometry(self.p_edges,x,G)
        self.ns=self.d0.ns
        self.nf=len(self.wd)
        self.np=len(self.wp)
        self.iP=self.ns+self.nf
        self.n=self.iP+self.np
        self.w=np.r_[self.ds.w,self.wd,self.wp]
        hist=[np.asarray(arc['endpoint_'+k+'_occupations']) for k in LABELS]
        self.initial=np.r_[hist[0],np.repeat(hist[1],2*n),np.repeat(hist[2],2*n)]
        self.seed_origin=[e0,self.d_edges,self.p_edges]
        old=np.r_[self.d0.d.w,self.d0.d.wf,self.d0.pw]*np.concatenate(hist)
        for o,p in zip(np.split(old,(self.ns,self.d0.iP)),
                     np.split(self.w*self.initial,(self.ns,self.iP)),strict=True):
            if abs(fsum(o)-fsum(p))>2e-14*max(abs(fsum(o)),1e-20):
                raise ArithmeticError('Archive count not preserved')

    def energies(self,x,m):
        return (np.hypot(self.ds.p,m),np.hypot(self.qd,x*self.ratio),np.hypot(self.qp,x))

    def state(self,x,y):
        m=self.ds.gap(x,y[:self.ns])
        E=np.concatenate(self.energies(x,m))
        u=(m*m-(1e-9*x)**2)**2/(2*self.ds.lam)
        rho=fsum(self.w*y*E)-u
        H=np.sqrt(self.radiation/(self.radiation+rho))
        return m,E,rho,H

    def trace(self,t,y):
        x=X*np.exp(t)
        m,E,_,_=self.state(x,y)
        masses=np.r_[np.full(self.ns,(1e-9*x)**2),
                     np.full(self.nf,(self.ratio*x)**2),np.full(self.np,x*x)]
        return fsum(self.w*y*masses/E)

    def energy_difference(self,x,a,z,b):
        m,E,_,_=self.state(x,a)
        u,F,_,_=self.state(z,b)
        c=(z-x)*(z+x)
        d2=np.r_[np.full(self.ns,(u-m)*(u+m)),np.full(self.nf,c*self.ratio**2),np.full(self.np,c)]
        U=((u-m)*(u+m)-1e-18*c)*(u*u+m*m-1e-18*(x*x+z*z))/(2*self.ds.lam)
        return fsum(self.w*((b-a)*E+b*d2/(E+F)))-U

    def chart(self,x,m,order):
        Es,Ed,Ep=self.energies(x,m)
        mb=self.ratio*x
        A=(x*x+m*m-mb*mb)/2
        D=np.sqrt((x*x-(m+mb)**2)*(x*x-(m-mb)**2))/2
        lower=(x*x*Es*Es+D*D)/(A*Es+self.ds.p*D)
        upper=lower+2*self.ds.p*D/(m*m)
        dE=np.hypot(self.d_edges,mb)
        pE=np.hypot(self.p_edges,x)
        a=np.maximum.reduce([lower,np.full_like(lower,pE[0]),Es+dE[0]])
        b=np.minimum.reduce([upper,np.full_like(upper,pE[-1]),Es+dE[-1]])
        pp=np.broadcast_to(pE,(self.ns,len(pE)))
        lim=np.c_[pp,Es[:,None]+dE]
        cuts=np.sort(np.c_[a,b,np.clip(lim,a[:,None],np.maximum(a,b)[:,None])],axis=1)
        left,right=cuts[:,:-1],cuts[:,1:]
        i,j=np.where((right>left)&(a<b)[:,None])
        w=(right[i,j]-left[i,j])/2
        center=(right[i,j]+left[i,j])/2
        v,t=roots_legendre(order)
        P=(center[:,None]+w[:,None]*v).ravel()
        S=np.repeat(i,order)
        B=P-Es[S]
        J=2*np.searchsorted(dE,B,side='right')-1
        K=2*np.searchsorted(pE,P,side='right')-1
        if min(J)<1 or max(J)>=self.nf or min(K)<1 or max(K)>=self.np:
            raise ArithmeticError('An event left the endpoint chart')
        ad=(B-Ed[J-1])/(Ed[J]-Ed[J-1])
        ap=(P-Ep[K-1])/(Ep[K]-Ep[K-1])
        if min(ad.min(),ap.min())<0 or max(ad.max(),ap.max())>1:
            raise ArithmeticError('Pauli deposition weights not convex')
        root=KAPPA*G*x*((x+mb)**2-m*m)/(Es*self.ds.p)*self.ds.w
        W=(root[i,None]*w[:,None]*t).ravel()
        Q=np.abs((Es[S]+(1-ad)*Ed[J-1]+ad*Ed[J]-(1-ap)*Ep[K-1]-ap*Ep[K])/P)
        if max(Q)>2e-14:
            raise ArithmeticError('Event energy projection failed')
        return (Es,Ed,Ep,S,J,K,ad,ap,W),float(max(Q))


def initial_checks(d,x):
    m,E,_,H=d.state(x,d.initial)
    ch,proj=d.chart(x,m,8)
    y=np.r_[1/np.expm1(E[:d.ns]),expit(-E[d.ns:])]
    _,f,b=event(ch,y,d)
    bal=float(max(abs(f-b)/np.maximum(f,1e-99)))
    z,l,r=event(ch,d.initial,d)
    fs,fd,fp=np.split(d.initial,(d.ns,d.iP))
    _,_,_,i,j,k,_,_,W=ch
    direct=np.array([fsum(W*fp[k]*(1+fs[i])*(1-fd[j])),fsum(W*(1-fp[k])*fs[i]*fd[j])])
    out=np.array([fsum(l),fsum(r)])
    replay=float(max(abs(direct/out-1)))
    if max(bal,replay)>3e-13 or abs(fsum(E*z))>2e-14:
        raise ArithmeticError(('Endpoint thermal or direct step control',bal,replay))
    return dict(max_event_energy_error=proj,max_thermal_defect=bal,
                direct_step_relative_defect=replay,forward_inverse=out.tolist(),
                initial_log_clock=H*fsum(z[:d.ns])/fsum(d.w[:d.ns]*fs))


def advance(d,x,y,span,steps):
    u=y.copy()
    t=np.linspace(0,span,steps+1)
    tr=[]
    for a,b in zip(t[:-1],t[1:],strict=True):
        c=(a+b)/2
        v=x*np.exp(c)
        m,_,_,H=d.state(v,u)
        A,_=d.chart(v,m,3)
        p,l=step(A,d,u,H*(b-a)/2)
        m,_,_,H=d.state(v,p)
        A,_=d.chart(v,m,3)
        z,l=step(A,d,u,H*(b-a),rate_state=p)
        tr.append(hold_trace(d,x,a,c,u)+hold_trace(d,x,c,b,z))
        u=z
    W=d.w*(u-y)
    N=[fsum(z) for z in np.split(W,(d.ns,d.iP))]
    diff=d.energy_difference(x,y,x*np.exp(span),u)-fsum(tr)
    if max(u[d.ns:])>1 or min(u)<0 or max(abs(N[0]-N[1]),abs(N[0]+N[2]))>2e-13 or abs(diff)>3e-13:
        raise ArithmeticError('Live endpoint number, energy or Pauli control failed')
    return dict(steps=steps,delta_species=N,energy_residual=diff,max_fermion=float(max(u[d.ns:])))


def run(g,nlevels):
    a=load(g)[0]
    x=a['x_end']
    results=[]
    for n in nlevels:
        d=CellChart(a,SPECS[g],n)
        tests=initial_checks(d,x)
        attempts=[advance(d,x,d.initial,5e-5,s) for s in (32,64)]
        z,q=[f['delta_species'][0] for f in attempts]
        if abs(q/z-1)>3e-4:
            raise ArithmeticError(('Live endpoint time differences not resolved',g,n,z,q))
        results.append(dict(subintervals_per_native_bin=n,counts=[d.ns,d.nf,d.np],initial=tests,live=attempts))
    return dict(grid=g,scope='Local opt-in phase-space control; no full evolving archive',results=results)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('grid',choices=SPECS)
    p.add_argument('--levels',nargs='+',type=int,default=[1,2,4])
    a=p.parse_args()
    print(json.dumps(run(a.grid,a.levels),indent=2,allow_nan=False))
