"""Frozen bound orbit conservative nonuniform phase transport, paper 28."""
import argparse,json
import numpy as np
from scipy.optimize import brentq,minimize_scalar
from scipy.interpolate import PchipInterpolator
from numpy.polynomial.legendre import leggauss

def minmod(a,b):return np.where(a*b>0,np.sign(a)*np.minimum(abs(a),abs(b)),0.)
class Orbit:
    def __init__(self,psi,s,L,rmax,rc,nz=256):
        self.L=L;self.s=s;self.psi=psi
        g=lambda r:2*r*r*(psi(r)-s)-L
        lo=brentq(g,0,rc,xtol=1e-12);hi=brentq(g,rc,rmax,xtol=1e-12);self.lo=lo;self.hi=hi
        z=np.linspace(0,np.pi/2,nz+1);x,w=leggauss(6);a=(z[:-1,None]+z[1:,None])/2+np.diff(z)[:,None]*x/2;rr=lo+(hi-lo)*np.sin(a)**2
        vr=np.sqrt(np.maximum(2*(psi(rr)-s)-L/rr**2,1e-300));d=(hi-lo)*np.sin(2*a)/vr
        time=np.r_[0,np.cumsum(np.diff(z)/2*(d@w))];self.T=2*time[-1]
        self.inv=PchipInterpolator(time,lo+(hi-lo)*np.sin(z)**2);self.forward=PchipInterpolator(lo+(hi-lo)*np.sin(z)**2,time)
    def coords(self,phase):
        phase=np.mod(phase,1);time=self.T*np.minimum(phase,1-phase);r=self.inv(time);p=np.sign(.5-phase)*np.sqrt(np.maximum(self.psi(r)-self.s-self.L/(2*r*r),0));return r,p
    def exterior(self,R):
        if R<=self.lo:return 0.,1.
        if R>=self.hi:return .5,.5
        a=float(self.forward(R)/self.T);return a,1-a
def mesh(n,k):
    assert n%2==0
    u=np.linspace(0,1,n//2+1);a=.5*np.expm1(k*u)/np.expm1(k) if k else .5*u
    return np.r_[a,1-a[-2::-1]]
def slopes(edges,f,order):
    if order==1:return np.zeros_like(f)
    d=np.diff(edges);x=(edges[1:]+edges[:-1])/2;left=np.r_[x[-1]-1,x[:-1]];right=np.r_[x[1:],x[0]+1]
    s=minmod((f-np.roll(f,1))/(x-left),(np.roll(f,-1)-f)/(right-x));return np.clip(s,-2*f/d,2*f/d)
def integrator(edges,f,order):
    d=np.diff(edges);s=slopes(edges,f,order);prefix=np.r_[0,np.cumsum(f*d)];n=len(f)
    def primitive(x):
        x=np.asarray(x);turn=np.floor(x);u=x-turn;j=np.clip(np.searchsorted(edges,u,side='right')-1,0,n-1);l=u-edges[j]
        return turn*prefix[-1]+prefix[j]+l*f[j]+s[j]*l*(l-d[j])/2
    return primitive
def shift(edges,f,dt,order=2):
    d=np.diff(edges);s=slopes(edges,f,order);n=len(f);delta=dt%1;lo=edges[:-1]-delta;hi=edges[1:]-delta
    turn0=np.floor(lo);turn1=np.floor(hi);a=lo-turn0;b=hi-turn1;i=np.clip(np.searchsorted(edges,a,side='right')-1,0,n-1);j=np.clip(np.searchsorted(edges,b,side='right')-1,0,n-1);same=(i==j)&(turn0==turn1)
    A=a-edges[i];B=b-edges[j];prefix=np.r_[0,np.cumsum(f*d)];v=(d[i]-A)*(f[i]+s[i]*A/2)+B*(f[j]+s[j]*(B-d[j])/2)+prefix[j]-prefix[i+1]+(turn1-turn0)*prefix[-1]
    v[same]=(B[same]-A[same])*(f[i[same]]+s[i[same]]*(A[same]+B[same]-d[i[same]])/2)
    return v/d
def halo(nL,nphase,t,nz,steps,k,order):
    from lowered_antonov_hessian import Halo
    h=Halo(model='linear');e0=h.Qt/h.rt;s=-.9*e0;rmax=h.Qt/(e0+s)
    def psi(r):
        r=np.asarray(r);inner=h.psi(np.minimum(r,h.rt).ravel()).reshape(r.shape);return np.where(r<=h.rt,inner,h.Qt/np.maximum(r,h.rt)-e0)
    sol=minimize_scalar(lambda r:-2*r*r*(psi(r)-s),bounds=(0,rmax),method='bounded');rc=sol.x;Lmax=-sol.fun
    x,w=leggauss(nL);Ls=(x+1)*Lmax/2;Ws=w*Lmax/2;rb=.8*h.rt;sr=.03*h.rt;mb=.8;sm=.05
    x,w=leggauss(6);edges=mesh(nphase,k);d=np.diff(edges);theta=(edges[:-1,None]+edges[1:,None])/2+d[:,None]*x[None,:]/2
    initial=final=ext=initialL=finalL=0.;minf=1e100;max_mass_err=0.
    for L,W in zip(Ls,Ws):
        o=Orbit(psi,s,L,rmax,rc,nz);r,p=o.coords(theta);mu=p/np.sqrt(psi(r)-s);f=np.exp(-.5*((r-rb)/sr)**2-.5*((mu-mb)/sm)**2)@w/2;g=f.copy();mass=W*o.T/(2*np.sqrt(2))
        for _ in range(steps):
            g=shift(edges,g,t/(o.T*steps),order);minf=min(minf,float(g.min()));g=np.maximum(g,0)
        a,b=o.exterior(h.rt);prim=integrator(edges,g,order);ini=mass*(d@f);fin=mass*(d@g);ext+=mass*(prim(b)-prim(a));initial+=ini;final+=fin;initialL+=L*ini;finalL+=L*fin
        if ini>1e-20:max_mass_err=max(max_mass_err,abs(fin/ini-1))
    return dict(nL=nL,nphase=nphase,nz=nz,t=t,steps=steps,k=k,order=order,mass=initial,mass_rel=final/initial-1,L2_rel=finalL/initialL-1,outside=ext/initial,min_before_clip=minf,max_orbit_mass_rel=max_mass_err)
def oscillator(n,t,steps,k,order):
    edges=mesh(n,k);d=np.diff(edges);x,w=leggauss(8);theta=(edges[:-1,None]+edges[1:,None])/2+d[:,None]*x[None,:]/2;fun=lambda theta:1-.4*np.sqrt(.7)*np.cos(2*np.pi*theta);f=fun(theta)@w/2;g=f.copy()
    for _ in range(steps):g=shift(edges,g,t/(np.pi*steps),order)
    ref=fun(theta-t/np.pi)@w/2;return dict(n=n,k=k,t=t,steps=steps,order=order,l1=float(d@abs(g-ref)),mass_rel=float((d@g)/(d@f)-1),min=float(g.min()))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--nL',type=int,default=128);p.add_argument('--nphase',type=int,default=512);p.add_argument('--nz',type=int,default=512);p.add_argument('--t',type=float,default=10);p.add_argument('--steps',type=int,default=1);p.add_argument('--k',type=float,default=8);p.add_argument('--order',type=int,default=2);p.add_argument('--oscillator',action='store_true');a=p.parse_args();print(json.dumps(oscillator(a.nphase,a.t,a.steps,a.k,a.order) if a.oscillator else halo(a.nL,a.nphase,a.t,a.nz,a.steps,a.k,a.order)))
