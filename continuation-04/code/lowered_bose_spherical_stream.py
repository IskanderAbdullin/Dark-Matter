"""Conservative fixed-binding spherical Liouville streaming.
Closed turning surface or finite open radial outer boundary, no collisions.
Exact discrete phase-space divergence cancellation, physical angular momentum
is monitored but is not exactly preserved by the finite-volume scheme.
"""
import argparse,json
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad,solve_ivp
from scipy.optimize import brentq

def minmod(a,b):return np.where(a*b>0,np.sign(a)*np.minimum(abs(a),abs(b)),0.)

class Grid:
    def __init__(self,psi,s,edges,nmu,open_outer=False):
        self.psi=psi;self.s=s;self.r=np.asarray(edges);self.nr=len(edges)-1;self.nmu=nmu;self.open=open_outer
        self.mu=np.linspace(-1,1,nmu+1);self.mc=(self.mu[:-1]+self.mu[1:])/2;self.dm=2/nmu
        self.rc=(self.r[1:]+self.r[:-1])/2;self.dr=np.diff(self.r)
        x,w=leggauss(24);t=(x+1)/2;w=w/2
        rr=self.r[1:,None]-self.dr[:,None]*t[None,:]**2;ww=2*self.dr[:,None]*t*w
        E=np.maximum(psi(rr)-s,0);self.rquad=rr;self.Equad=E;self.weight=ww*rr*rr*np.sqrt(E)
        self.V=self.weight.sum(1);self.mass=self.V[:,None]*self.dm
        self.A=np.sqrt(2)*self.r*self.r*np.maximum(psi(self.r)-s,0);self.A[0]=0
        if not self.open:self.A[-1]=0
        self.B=np.diff(self.A)/2;self.Rflux=self.A[:,None]*self.mc[None,:]*self.dm;self.Uflux=self.B[:,None]*(1-self.mu[None,:]**2)
        outgoing=np.maximum(self.Rflux[1:],0)-np.minimum(self.Rflux[:-1],0)+np.maximum(self.Uflux[:,1:],0)-np.minimum(self.Uflux[:,:-1],0)
        self.dt=.4/np.max(outgoing/self.mass)
        self.L2=self.ravg(2*rr*rr*E)[:,None]*(1-(self.mu[1:]**3-self.mu[:-1]**3)/(3*self.dm))[None,:]
    def ravg(self,v):return (self.weight*v).sum(1)/self.V
    def total(self,f):return float(np.sum(self.mass*f))
    def rhs(self,f,order=2):
        rgrad=np.zeros_like(f);ugrad=np.zeros_like(f)
        if order==2:
            rgrad[1:-1]=minmod((f[1:-1]-f[:-2])/np.diff(self.rc)[:-1,None],(f[2:]-f[1:-1])/np.diff(self.rc)[1:,None])
            ugrad[:,1:-1]=minmod(f[:,1:-1]-f[:,:-2],f[:,2:]-f[:,1:-1])/self.dm
        rL=np.vstack((f[0:1],f+rgrad*self.dr[:,None]/2));rR=np.vstack((f-rgrad*self.dr[:,None]/2,np.zeros_like(f[-1:])))
        uL=np.hstack((f[:,0:1],f+ugrad*self.dm/2));uR=np.hstack((f-ugrad*self.dm/2,f[:,-1:]))
        Fr=self.Rflux*np.where(self.Rflux>=0,rL,rR);Fu=self.Uflux*np.where(self.Uflux>=0,uL,uR)
        return -(np.diff(Fr,axis=0)+np.diff(Fu,axis=1))/self.mass,float(Fr[-1].sum())
    def evolve(self,f,target,order=2):
        t=0.;out=0.;steps=0;reject=0;max_dt=self.dt
        while t<target-1e-13:
            d=min(max_dt,target-t);g,o=self.rhs(f,order)
            while True:
                p=f+d*g;h,z=self.rhs(p,order);new=.5*(f+p+d*h)
                if p.min()>=-1e-14 and new.min()>=-1e-14:break
                d*=.5;reject+=1;max_dt=min(max_dt,d)
            f=new;out+=d*(o+z)/2;t+=d;steps+=1
        return f,dict(steps=steps,rejects=reject,outflow=out,dt_cfl=self.dt)

def oscillator(nr,nmu,t,order):
    s=1.;psi=lambda r:2-r*r/2;G=Grid(psi,s,np.linspace(0,np.sqrt(2),nr+1),nmu)
    rr=G.rquad;E=G.Equad;A=G.ravg(rr*rr-2*E)[:,None];B=G.ravg(2*rr*np.sqrt(2*E))[:,None]*G.mc[None,:]
    exact=lambda t:1+.2*(A*np.cos(2*t)-B*np.sin(2*t))+np.zeros((nr,nmu))
    f0=exact(0);f,run=G.evolve(f0.copy(),t,order);true=exact(t)
    c=np.ones_like(f);g,o=G.rhs(c,order)
    stationary=1+.2*G.L2;fs,st=G.evolve(stationary.copy(),t,order)
    norm=G.total(np.ones_like(f));L=lambda f:G.total(f*G.L2)
    return dict(nr=nr,nmu=nmu,t=t,order=order,well_balanced=float(abs(g).max()),l1=float(G.total(abs(f-true))/norm),l2=float(np.sqrt(G.total((f-true)**2)/norm)),linf=float(abs(f-true).max()),mass_rel=(G.total(f)-G.total(f0))/G.total(f0),L2_rel=(L(f)-L(f0))/L(f0),stationary_L2_l1=float(G.total(abs(fs-stationary))/norm),stationary_L2_moment_rel=(L(fs)-L(stationary))/L(stationary),min=float(f.min()),**run)

def halo_test(model,nr,nmu,age,order,unbound=False):
    from lowered_antonov_hessian import Halo
    H=Halo(model=model);e0=H.Qt/H.rt
    def psi(r):
        r=np.asarray(r);inside=H.psi(np.minimum(r,H.rt).ravel()).reshape(r.shape);outside=H.Qt/np.maximum(r,H.rt)-e0
        return np.where(r<=H.rt,inside,outside)
    eps=(-.2 if unbound else .1)*e0;s=eps-e0;Rmax=2*H.rt
    edges=np.r_[np.linspace(0,H.rt,nr//2+1),np.linspace(H.rt,Rmax,nr//2+1)[1:]]
    G=Grid(psi,s,edges,nmu,True)
    rb=.8*H.rt;mb=.8;sigmar=.03*H.rt;sigmamu=.05
    radial=G.ravg(np.exp(-.5*((G.rquad-rb)/sigmar)**2));x,w=leggauss(12);mu=G.mc[:,None]+G.dm*x[None,:]/2;angular=np.sum(w[None,:]/2*np.exp(-.5*((mu-mb)/sigmamu)**2),axis=1)
    f0=radial[:,None]*angular[None,:];f0/=G.total(f0)
    f,run=G.evolve(f0.copy(),age,order);outside=G.total(f*(G.rc[:,None]>=H.rt));inside=G.total(f*(G.rc[:,None]<H.rt))
    x,w=leggauss(32);r0=rb+sigmar*5*x;mu0=(mb-5*sigmamu+1)/2+(1-mb+5*sigmamu)*x/2
    ri,mi=np.meshgrid(r0,mu0,indexing='ij');wi,wj=np.meshgrid(w,w,indexing='ij');ri=ri.ravel();mi=mi.ravel();v=np.sqrt(2*(psi(ri)-s));L2=(ri*v)**2*(1-mi*mi);weights=(wi*wj).ravel()*ri**2*v*np.exp(-.5*((ri-rb)/sigmar)**2-.5*((mi-mb)/sigmamu)**2);weights/=weights.sum()
    def acc(r):
        rr=np.maximum(r,1e-10);interior=-H.sol.sol(np.clip(rr,H.eps,H.rt))[1]/rr**2;interior=np.where(rr<H.eps,-rr/3,interior)
        return np.where(rr<=H.rt,interior,-H.Qt/rr**2)
    n=len(ri)
    def rhs(t,y):
        r=y[:n];vr=y[n:];return np.r_[vr,acc(r)+L2/r**3]
    sol=solve_ivp(rhs,[0,age],np.r_[ri,v*mi],rtol=1e-9,atol=1e-10,method='DOP853',max_step=2.)
    rfinal=sol.y[:n,-1];vfinal=sol.y[n:,-1];assert sol.success
    char_out=float(weights@(rfinal>=Rmax));char_ext=float(weights@((rfinal>=H.rt)&(rfinal<Rmax)));char_in=float(weights@(rfinal<H.rt))
    E0=.5*(v*mi)**2+.5*L2/ri**2-psi(ri);E1=.5*vfinal**2+.5*L2/rfinal**2-psi(rfinal)
    return dict(model=model,nr=nr,nmu=nmu,age=age,order=order,unbound=unbound,rt=H.rt,e0=e0,s=s,inside=inside,outside=outside,char_inside=char_in,char_outside=char_ext,char_outflow=char_out,ledger=inside+outside+run['outflow']-1,initial_L2=G.total(f0*G.L2),remaining_L2=G.total(f*G.L2),characteristic_energy_error=float(abs(E1-E0).max()),**run)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--nr',type=int,default=64);p.add_argument('--nmu',type=int,default=32);p.add_argument('--t',type=float,default=1.);p.add_argument('--order',type=int,default=2);p.add_argument('--model',default='oscillator');p.add_argument('--unbound',action='store_true');a=p.parse_args();print(json.dumps(oscillator(a.nr,a.nmu,a.t,a.order) if a.model=='oscillator' else halo_test(a.model,a.nr,a.nmu,a.t,a.order,a.unbound)))
