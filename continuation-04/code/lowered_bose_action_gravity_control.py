"""Canonical fixed-L residual Hamiltonian transport. Harmonic analytic control.
Coordinates e=reference energy, normalized phase theta. J=pi(e-sqrt(L)).
The discrete streamfunction face flux is exactly divergence-free.
"""
import argparse,json
import numpy as np
from numpy.polynomial.legendre import leggauss

def mm(a,b):return np.where(a*b>0,np.sign(a)*np.minimum(abs(a),abs(b)),0.)

class Grid:
 def __init__(self,NE,NT,L=.3,Emax=8.,alpha=.2,full=True):
    self.NE=NE;self.NT=NT;self.L=L;self.alpha=alpha;self.full=full;self.e=np.linspace(np.sqrt(L),Emax,NE+1);self.th=np.linspace(0,1,NT+1);self.de=np.diff(self.e);self.dt=1/NT;self.ec=(self.e[1:]+self.e[:-1])/2;self.tc=(self.th[1:]+self.th[:-1])/2
    self.mass=np.pi*self.de[:,None]*self.dt;E=self.e[:,None];T=self.th[None,:];r2=E+np.sqrt(np.maximum(E*E-L,0))*np.cos(2*np.pi*T);H=alpha*r2/2+(E if full else 0)
    self.Fe=-np.diff(H,axis=1);self.Ft=np.diff(H,axis=0);self.Fe[0]=0;self.Ft[:,-1]=self.Ft[:,0];self.cfl=.35/np.max((np.maximum(self.Fe[1:],0)-np.minimum(self.Fe[:-1],0)+np.maximum(self.Ft[:,1:],0)-np.minimum(self.Ft[:,:-1],0))/self.mass)
 def coords(self,e,th):
    A=np.sqrt(np.maximum(e*e-self.L,0));r=np.sqrt(e+A*np.cos(2*np.pi*th));vr=-A*np.sin(2*np.pi*th)/r;vt=np.sqrt(self.L)/r;return r,vr,vt
 def avg(self,fun,n=6):
    x,w=leggauss(n);e=self.ec[:,None,None,None]+self.de[:,None,None,None]*x[None,None,:,None]/2;th=self.tc[None,:,None,None]+self.dt*x[None,None,None,:]/2;r,vr,vt=self.coords(e,th);return np.sum(fun(r,vr,vt,e)*w[None,None,:,None]*w[None,None,None,:]/4,axis=(2,3))
 def initfun(self,r,vr,vt,e):
    return np.exp(-e)*np.maximum(0,1-(e/4)**4)**2*(1+.15*(r*r-vr*vr-vt*vt)/(1+e))
 def exact(self,t):
    a=self.alpha
    def fn(r,vr,vt,e):
        if self.full:
            om=np.sqrt(1+a);c=np.cos(om*t);s=np.sin(om*t);x=r*c-vr*s/om;y=-vt*s/om;vx=vr*c+om*r*s;vy=vt*c
        else:x=r;y=np.zeros_like(r);vx=vr+a*r*t;vy=vt
        E=.5*(x*x+y*y+vx*vx+vy*vy);return self.initfun(np.sqrt(x*x+y*y),vx,vy,E)
    return self.avg(fn,n=8)
 def rhs(self,f,order=2,closed=False):
    de=self.de;dt=self.dt;es=np.zeros_like(f);ts=np.zeros_like(f)
    if order==2:
        es[1:-1]=mm((f[1:-1]-f[:-2])/np.diff(self.ec)[:-1,None],(f[2:]-f[1:-1])/np.diff(self.ec)[1:,None]);ts=mm(f-np.roll(f,1,axis=1),np.roll(f,-1,axis=1)-f)/dt
        es=np.clip(es,-2*f/de[:,None],2*f/de[:,None]);ts=np.clip(ts,-2*f/dt,2*f/dt)
    el=np.vstack((f[0:1],f+es*de[:,None]/2));er=np.vstack((f-es*de[:,None]/2,np.zeros_like(f[-1:])))
    tl=np.hstack((f[:,-1:]+ts[:,-1:]*dt/2,f+ts*dt/2));tr=np.hstack((f-ts*dt/2,f[:,0:1]-ts[:,0:1]*dt/2))
    Fe=self.Fe*np.where(self.Fe>=0,el,er);Ft=self.Ft*np.where(self.Ft>=0,tl,tr)
    if closed:Fe[-1]=0
    return -(np.diff(Fe,axis=0)+np.diff(Ft,axis=1))/self.mass,float(Fe[-1].sum())
 def evolve(self,f,t,order=2):
    out=0.;T=0.;steps=0;reject=0;dmax=self.cfl
    while T<t-1e-13:
        d=min(dmax,t-T);g,o=self.rhs(f,order)
        while True:
            p=f+d*g;h,z=self.rhs(p,order);q=.5*(f+p+d*h)
            if p.min()>=-1e-14 and q.min()>=-1e-14:break
            d*=.5;reject+=1;dmax=min(dmax,d)
        f=q;out+=d*(o+z)/2;T+=d;steps+=1
    return f,out,steps,reject
 def total(self,f):return float(np.sum(f*self.mass))

def run(NE,NT,L,alpha,t,order,full):
    G=Grid(NE,NT,L,8,alpha,full);f0=G.avg(G.initfun,8);f,out,steps,reject=G.evolve(f0,t,order);ref=G.exact(t);H=G.avg(lambda r,vr,vt,e: e+alpha*r*r/2 if full else alpha*r*r/2,8);N=G.total(f0);h0=G.total(H*f0);h1=G.total(H*f);htrue=G.total(H*ref)
    div=np.diff(G.Fe,axis=0)+np.diff(G.Ft,axis=1)
    return dict(NE=NE,NT=NT,L=L,alpha=alpha,t=t,order=order,full=full,steps=steps,reject=reject,max_flux_div=float(abs(div).max()),number_ledger=float((G.total(f)+out)/N-1),out=float(out/N),min=float(f.min()),L1=float(G.total(abs(f-ref))/N),Linf=float(abs(f-ref).max()),energy_relative=float(h1/h0-1),energy_quadrature_exact_error=float(htrue/h0-1),energy_ref_difference=float((h1-htrue)/h0),initial_N=N)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--NE',type=int,default=64);p.add_argument('--NT',type=int,default=64);p.add_argument('--L',type=float,default=.3);p.add_argument('--alpha',type=float,default=.2);p.add_argument('--t',type=float,default=1);p.add_argument('--order',type=int,default=2);p.add_argument('--kick',action='store_true');a=p.parse_args();print(json.dumps(run(a.NE,a.NT,a.L,a.alpha,a.t,a.order,not a.kick)))
