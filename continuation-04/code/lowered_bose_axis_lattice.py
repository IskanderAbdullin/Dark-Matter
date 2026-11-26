"""Conservative fully nonlinear axisymmetric UU energy/axial momentum lattice.
The conditional cell-smoothed elliptic kernel is an approximation."""
import argparse,json
import numpy as np
from scipy.stats import qmc
from numpy.polynomial.legendre import leggauss
from lowered_bose_axis_kernel import kernel
from lowered_bose_quadrupole import occ

def offsets(m,seed):
    z=qmc.Sobol(6,scramble=True,seed=seed).random_base2(m)-.5
    a=z[:,:3];b=z[:,3:];a=np.c_[a,a[:,0]+a[:,1]-a[:,2]];b=np.c_[b,b[:,0]+b[:,1]-b[:,2]]
    keep=(abs(a[:,3])<.5)&(abs(b[:,3])<.5);a=a[keep];b=b[keep]
    return np.concatenate([a,a]),np.concatenate([b,-b])

def cell_moments(elo,ehi,plo,phi):
    a=np.sqrt(elo);b=np.sqrt(ehi);cuts=np.unique(np.clip([plo,phi,-b,-a,a,b],plo,phi));x,w=leggauss(4);ans=np.zeros(4)
    for l,r in zip(cuts[:-1],cuts[1:]):
        p=(l+r)/2+(r-l)*x/2;lower=np.maximum(elo,p*p);width=np.maximum(ehi-lower,0);wt=w*(r-l)/8
        ans+=np.array([wt@width,wt@(.5*(ehi+lower)*width),wt@(p*width),wt@(p*p*width)])
    return ans

class Network:
    def __init__(self,NE,K,Emax,m=4,seed=0,reg=1.,mode='all',geometry='cutcell'):
        self.NE=NE;self.K=K;self.de=de=Emax/NE;self.dp=dp=np.sqrt(Emax)/K;self.Emax=Emax
        ie,ip=np.meshgrid(np.arange(NE),np.arange(-K,K+1),indexing='ij');ie=ie.ravel();ip=ip.ravel();E=(ie+.5)*de;p=ip*dp;mom=np.array([cell_moments(e*de,(e+1)*de,(z-.5)*dp,(z+.5)*dp) for e,z in zip(ie,ip)])
        keep=mom[:,0]>0 if geometry=='cutcell' else E>p*p
        self.ie=ie=ie[keep];self.ip=ip=ip[keep];self.E=E[keep];self.p=p[keep];self.N=N=len(self.E);mom=mom[keep]
        if geometry=='cutcell':self.D=mom[:,0];self.Ebar=mom[:,1]/self.D;self.p2bar=mom[:,3]/self.D
        else:self.D=np.full(N,de*dp/4);self.Ebar=self.E;self.p2bar=self.p**2
        self.H=1.5*self.p2bar-.5*self.Ebar
        i,j=np.triu_indices(N);key=(ie[i]+ie[j])*(4*K+1)+ip[i]+ip[j]+2*K;order=np.argsort(key);i=i[order];j=j[order];key=key[order];cuts=np.r_[0,np.flatnonzero(np.diff(key))+1,len(i)]
        ids=[]
        for lo,hi in zip(cuts[:-1],cuts[1:]):
            a,b=np.triu_indices(hi-lo,1);ids.append(np.array([i[lo+a],j[lo+a],i[lo+b],j[lo+b]]))
        ids=np.concatenate(ids,axis=1).astype(np.int32);self.raw_events=ids.shape[1];allw=[];sel=[];singular=0;count=0
        a,b=offsets(m,seed);self.samples=len(a)
        for lo in range(0,ids.shape[1],100000):
            z=ids[:,lo:lo+100000];E=self.E[z];p=self.p[z];v=kernel(E,p);bad=~np.isfinite(v);singular+=int(bad.sum());count+=len(v)
            if mode=='all':mask=np.ones_like(bad)
            elif mode=='singular':mask=bad
            elif mode=='near':mask=bad|(v>1.)
            else:raise ValueError(mode)
            if np.any(mask):
                EE=E[:,mask];PP=p[:,mask];val=np.zeros(EE.shape[1])
                for aa,bb in zip(a,b):
                    w=kernel(EE+reg*de*aa[:,None],PP+reg*dp*bb[:,None]);val+=np.where(np.isfinite(w),w,0.)
                v[mask]=val/len(a)
            w=(np.pi/16)*de**3*dp**3*v*np.where(z[0]==z[1],1,2)*np.where(z[2]==z[3],1,2)
            keep=w>0;sel.append(z[:,keep]);allw.append(w[keep])
        self.ids=np.concatenate(sel,axis=1);self.w=np.concatenate(allw);self.singular_fraction=singular/count
    def deriv(self,f,bose=True):
        i,j,k,l=self.ids;a,b,c,d=f[self.ids];net=self.w*(a*b*(1+c+d)-c*d*(1+a+b)) if bose else self.w*(a*b-c*d)
        g=np.bincount(k,weights=net,minlength=self.N)+np.bincount(l,weights=net,minlength=self.N)-np.bincount(i,weights=net,minlength=self.N)-np.bincount(j,weights=net,minlength=self.N)
        return g/self.D
    def report(self,f):
        return dict(N=float(self.D@f),E=float((self.D*self.E)@f),P=float((self.D*self.p)@f),Q=float((self.D*self.H)@f),S=float(self.D@((1+f)*np.log1p(f)-f*np.log(np.maximum(f,1e-300)))),min=float(f.min()))

def run(NE,K,Emax,psi,model,m=4,seed=0,reg=1.,mode='all',tmax=.2,dt=.02,geometry='cutcell'):
    L=Network(NE,K,Emax,m,seed,reg,mode,geometry);E=L.E;p=L.p;eps=.01
    f0=lambda e:np.exp(-e) if model=='MB' else occ(e,psi,model)[0]
    f=f0(E);a=f0(E+np.expm1(-eps)*L.Ebar+(np.exp(2*eps)-np.exp(-eps))*L.p2bar);b=f0(E+np.expm1(eps)*L.Ebar+(np.exp(-2*eps)-np.exp(eps))*L.p2bar)
    bose=model!='MB';g=L.deriv(f,bose);ga=L.deriv(a,bose);gb=L.deriv(b,bose);rate=-float((L.D*L.H)@(ga-gb)/((L.D*L.H)@(a-b)))
    be=1/np.expm1(.5+E-.1*p);eq=float(np.max(abs(L.deriv(be))))
    out=dict(NE=NE,K=K,Emax=Emax,psi=psi,model=model,N=L.N,events=len(L.w),raw_events=L.raw_events,samples=L.samples,mode=mode,reg=reg,seed=seed,singular_fraction=L.singular_fraction,rate=rate,initial=L.report(f),geometry=geometry,phase_volume=float(L.D.sum()),initial_dN=float(L.D@g),initial_dE=float((L.D*E)@g),initial_dP=float((L.D*p)@g),isotropic_Qdot=float((L.D*L.H)@g),eq_res=eq,initial_outside_derivative=float(L.D[E>=psi]@g[E>=psi]))
    if not bose:return out
    f=a;out['evolving_initial']=L.report(f);t=0;steps=0;rejects=0
    while t<tmax-1e-12:
        g=L.deriv(f);neg=g<0;lim=np.min(-f[neg]/g[neg]) if neg.any() else np.inf;d=min(dt,tmax-t,.4*lim)
        while True:
            x=f+d*g;y=.5*(f+x+d*L.deriv(x))
            if x.min()>=0 and y.min()>=0 and L.report(y)['S']>=L.report(f)['S']-1e-12:break
            d*=.5;rejects+=1
        f=y;t+=d;steps+=1
    out.update(final=L.report(f),steps=steps,rejects=rejects,t=t,outside_final=float(L.D[E>=psi]@f[E>=psi]));return out
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--NE',type=int,default=16);p.add_argument('--K',type=int,default=8);p.add_argument('--Emax',type=float,default=12);p.add_argument('--psi',type=float,default=3.9);p.add_argument('--model',default='linear');p.add_argument('--m',type=int,default=4);p.add_argument('--seed',type=int,default=0);p.add_argument('--reg',type=float,default=1.);p.add_argument('--mode',default='all');p.add_argument('--tmax',type=float,default=.2);p.add_argument('--geometry',default='cutcell');a=p.parse_args();print(json.dumps(run(a.NE,a.K,a.Emax,a.psi,a.model,a.m,a.seed,a.reg,a.mode,a.tmax,geometry=a.geometry)))
