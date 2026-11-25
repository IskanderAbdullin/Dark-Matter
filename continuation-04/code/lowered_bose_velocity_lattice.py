"""Angularly resolved finite velocity lattice, exact pair momentum/energy.
Uniform lattice-shell outgoing measure is a numerical approximation,
not a proved continuum quadrature. Velocity boundary restricts outgoing
pairs symmetrically, not spatial escape.
"""
import argparse,json,time
import numpy as np
from lowered_bose_quadrupole import occ,response
from scipy.integrate import quad

class Lattice:
    def __init__(self,K,qmax):
        self.K=K;self.h=qmax/K;self.w=self.h**3/(4*np.pi)
        v=np.stack(np.meshgrid(*[np.arange(-K,K+1)]*3,indexing='ij'),axis=-1).reshape(-1,3)
        e=np.sum(v*v,axis=1);sel=e<=K*K;self.v=v[sel];self.e=e[sel];self.E=self.e*self.h**2;self.N=len(self.v)
        i,j=np.triu_indices(self.N,1);p=self.v[i]+self.v[j];s=self.e[i]+self.e[j];B=4*K+1
        key=(((p[:,0]+2*K)*B+p[:,1]+2*K)*B+p[:,2]+2*K)*(2*K*K+1)+s
        _,idx=np.unique(key,return_inverse=True);self.group=idx.astype(np.int32);self.i=i.astype(np.int32);self.j=j.astype(np.int32)
        self.rel=self.h*np.sqrt(2*s-np.sum(p*p,axis=1));self.count=np.bincount(idx);self.G=len(self.count)
        self.mult=self.count[idx];self.H=self.h*self.h*(1.5*self.v[:,2]**2-.5*self.e)
    def deriv(self,f,bose=True):
        i,j,g=self.i,self.j,self.group;a=f[i];b=f[j];prod=a*b;add=a+b if bose else np.zeros_like(a)
        P=np.bincount(g,weights=prod,minlength=self.G)/self.count;S=np.bincount(g,weights=add,minlength=self.G)/self.count
        net=self.w*self.rel*((1+add)*P[g]-prod*(1+S[g]))
        return np.bincount(i,weights=net,minlength=self.N)+np.bincount(j,weights=net,minlength=self.N)
    def report(self,f):
        return dict(N=float(self.w*f.sum()),E=float(self.w*(self.E@f)),mom=(self.w*self.h*(self.v.T@f)).tolist(),Q=float(self.w*self.H@f),S=float(self.w*np.sum((1+f)*np.log1p(f)-f*np.log(np.maximum(f,1e-300)))),min=float(f.min()))

def run(K,qmax,psi,model,eps=.01,tmax=.5,dt=.02):
    L=Lattice(K,qmax);e=L.E;z=L.v[:,2]**2*L.h**2
    def f0(e):return np.exp(-e) if model=='MB' else occ(e,psi,model)[0]
    f=f0(e);fp=f0(np.exp(-eps)*(e-z)+np.exp(2*eps)*z);fm=f0(np.exp(eps)*(e-z)+np.exp(-2*eps)*z)
    bose=model!='MB';d0=L.deriv(f,bose);dp=L.deriv(fp,bose);dm=L.deriv(fm,bose);r=L.report(f)
    rate=-float(L.H@(dp-dm)/(L.H@(fp-fm)))
    out=dict(K=K,qmax=qmax,psi=psi,model=model,nvel=L.N,npairs=len(L.i),groups=L.G,mean_shell=float(len(L.i)/L.G),singleton_fraction=float(np.mean(L.mult==1)),initial=r,rate=rate,initial_dN=float(L.w*d0.sum()),initial_dE=float(L.w*(e@d0)),initial_dP=(L.w*L.h*(L.v.T@d0)).tolist(),isotropic_Qdot=float(L.w*(L.H@d0)),eq_res=float(np.max(abs(L.deriv(1/np.expm1(.3+e))))))
    out['orientation_rates']={}
    for name,axis in [('face',[1,1,0]),('body',[1,1,1])]:
        axis=np.array(axis,dtype=float);axis/=np.linalg.norm(axis);z=(L.h*(L.v@axis))**2;H=1.5*z-.5*e
        a=f0(np.exp(-eps)*(e-z)+np.exp(2*eps)*z);b=f0(np.exp(eps)*(e-z)+np.exp(-2*eps)*z)
        out['orientation_rates'][name]=-float(H@(L.deriv(a,bose)-L.deriv(b,bose))/(H@(a-b)))
    if not bose:return out
    out['initial_outside_derivative']=float(L.w*d0[e>=psi].sum())
    f=fp;out['evolving_initial']=L.report(f);t=0;steps=0;rejects=0
    while t<tmax-1e-12:
        g=L.deriv(f);mask=g<0;lim=np.min(-f[mask]/g[mask]) if mask.any() else np.inf;d=min(dt,tmax-t,.4*lim)
        while True:
            p=f+d*g
            new=.5*(f+p+d*L.deriv(p))
            if new.min()>=0 and L.report(new)['S']>=L.report(f)['S']-1e-12:break
            d*=.5;rejects+=1
        f=new;t+=d;steps+=1
    out.update(t=t,steps=steps,rejects=rejects,final=L.report(f),outside_final=float(L.w*f[e>=psi].sum()))
    return out
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--K',type=int,default=5);p.add_argument('--qmax',type=float,default=3);p.add_argument('--psi',type=float,default=3.9);p.add_argument('--model',default='linear');p.add_argument('--tmax',type=float,default=.5);p.add_argument('--dt',type=float,default=.02);a=p.parse_args();print(json.dumps(run(a.K,a.qmax,a.psi,a.model,tmax=a.tmax,dt=a.dt)))
