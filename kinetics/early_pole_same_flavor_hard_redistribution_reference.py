"""Reversible hard redistribution within each charge-symmetric Dirac flavor."""
import argparse
import json

import numpy as np
from scipy.special import expit,roots_legendre

from early_pole_archived_shared_cold_window_reference import Cold
from early_pole_live_joint_hard_pair_packets_reference import gauss_each


def spin(t):
    u=1-t
    pole=(1+u*u)/t**2
    moller=pole+(1+t*t)/u**2+2/(t*u)
    bhabha=pole+t*t+u*u-2*u*u/t
    return moller,bhabha


class SameFlavorKernel:
    def __init__(self,M,order,mode,floor=.15,strength=1.):
        if not 0.<floor<.25 or strength<0.:
            raise ValueError("Bad hard angular floor or strength")
        self.M=M
        self.floor=floor
        self.mode=mode
        self.strength=strength
        self.fac=.005**4/(np.sqrt(np.pi**2*106.75/90)*(86/.05)/2.435e18)
        self.graph=[self.packets(k,order) for k in (1,2)]
        self.reject=[0.,0.]

    def packets(self,k,order):
        nr,ni,nc,na=order
        M=self.M
        grids=[np.asarray(M.edge_sets[i]) for i in (k,k)]
        (p,wp),(q,wq)=[gauss_each(e,nr) for e in grids]
        z,w=roots_legendre(ni)
        a,b=np.meshgrid(p,q,indexing='ij')
        lim=4*a*b
        l=np.minimum(1.,1/lim)
        u=l[...,None]+(1-l[...,None])*(1+z)/2
        I=np.outer(2*wp*p*p/np.pi**2,2*wq*q*q/np.pi**2)[...,None]*(1-l[...,None])*u*w/2
        ok=I>0
        W=I[ok]
        s=(lim[...,None]*u)[ok]
        a,b=(np.broadcast_to(p[...,None],I.shape)[ok] for p in (a,b))
        E=a+b
        K=np.sqrt(np.maximum(E*E-s,0.))
        A=np.clip(np.divide(a-b,K,out=np.zeros_like(K),where=K>0.),-1.,1.)
        z,w=roots_legendre(na)
        X,U=gauss_each(np.array([self.floor,.25,.5,.75,1-self.floor]),nc)
        packs=[]
        for x,v in zip(X,U,strict=True):
            ct=1-2*x
            V=(E+K*A*ct)/2
            R=np.maximum(K*np.sqrt((1-ct*ct)*(1-A*A))/2,1e-30)
            cuts=np.c_[(grids[0]-V[:,None])/R[:,None],
                       (E[:,None]-grids[1]-V[:,None])/R[:,None],
                       np.broadcast_to([-1.,1.],(len(E),2))]
            cuts=np.sort(np.arccos(np.clip(cuts,-1.,1.)),axis=1)
            hh=np.diff(cuts,axis=1)/2
            i,j=np.nonzero(hh>1e-14)
            p=V[i,None]+R[i,None]*np.cos((cuts[i,j]+cuts[i,j+1])[:,None]/2+hh[i,j,None]*z)
            q=E[i,None]-p
            polar=spin(x)
            mult=(.125*polar[0] if self.mode in ("both","moller") else 0.)
            mult+=(.25*polar[1] if self.mode in ("both","bhabha") else 0.)
            f=(W*(2*v)*mult/(16*np.pi*s))[i,None]*hh[i,j,None]*w/np.pi
            ok=(p>=grids[0][0])&(p<=grids[0][-1])&(q>=grids[1][0])&(q<=grids[1][-1])
            packs.append(np.array([np.broadcast_to(y[i,None],p.shape)[ok] for y in (a,b,K)]
                                  +[p[ok],q[ok],f[ok]]))
        a,b,K,p,q,W=np.concatenate(packs,axis=1)
        ix=np.stack([np.searchsorted(g[1:],x)+s.start for x,g,s in
            zip((a,b,p,q),grids*2,[M.r[k]]*4,strict=True)])
        x=np.hypot(np.stack((a,b,p,q)),M.cells[M.r[k].start].m)
        basis=(x-M.shape[ix,0])/M.shape[ix,1]
        assert max(np.ravel(abs(basis)))<=1.+1e-12
        return ix[:3],basis[:3],p,K,W

    def collision(self,e,F,grad):
        M=self.M
        out=np.zeros_like(e)
        entropy=[]
        for k,(ix,u,p,K,W) in enumerate(self.graph,1):
            span=M.r[k]
            ev=e[ix,0]+e[ix,1]*u
            aim=ev[0]+ev[1]-ev[2]
            j=span.start+np.searchsorted(e[span,0]+e[span,1],aim)
            valid=j<span.stop
            j=np.minimum(j,span.stop-1)
            v=(aim-e[j,0])/e[j,1]
            q=M.cells[span.start].pfrombare(M.shape[j,0]+M.shape[j,1]*v)
            ok=valid&(abs(v)<=1.)&(abs(p-q)<=K)&(p+q>=K)
            v=np.where(ok,v,0.)
            self.reject[k-1]=max(self.reject[k-1],float(W[~ok].sum()/W.sum()))
            a,b,c=F[ix,0]+F[ix,1]*u
            d=F[j,0]+F[j,1]*v
            A,B,C,D=map(expit,(a,b,c,d))
            aff=a+b-c-d
            fw=A*B*(1-C)*(1-D)
            bw=C*D*(1-A)*(1-B)
            rate=np.where(abs(aff)<10.,2*np.sqrt(fw*bw)*np.sinh(np.clip(aff,-10.,10.)/2),fw-bw)
            R=self.fac*self.strength*W*ok*rate
            entropy.append(float(R@aff))
            for j,v,r in zip([*ix,j],[*u,v],[-R,-R,R,R],strict=True):
                out[:,0]+=np.bincount(j,weights=r,minlength=len(out))
                out[:,1]+=np.bincount(j,weights=r*v,minlength=len(out))
        res=[float(np.sum(grad[z]*out[z])) for z in M.r[1:]]
        scale=max(1.,float(np.sum(abs(grad*out))))
        assert max(abs(z) for z in res)<4e-12*scale and min(entropy)>=-1e-13,(res,scale)
        assert max(abs(out[r,0].sum()) for r in M.r)<4e-12*max(1.,float(np.sum(abs(out))))
        return out,entropy,res


def audit(grid,level,mode,floor):
    order=[(1,3,4,1),(2,4,5,2)][level]
    M=Cold(grid,1,0,8,True)
    K=SameFlavorKernel(M.M,order,mode,floor)
    e,F,E,*_=M.M.weights(M.M.init,True)
    D,S,W=K.collision(e,F,E)
    Y=np.r_[M.M.init.ravel(),M.ns,0.,0.]
    V=np.r_[D.ravel(),np.zeros(M.n+2)]
    def scalar(z):
        return M.rhs(z,False)[M.j0:M.j0+M.n]
    tangent=[((scalar(Y+h*V)-scalar(Y))/h).tolist() for h in (1e-10,1e-9,1e-8)]
    for i,span in enumerate(M.M.r[1:]):
        F[span]=-E[span]*(1. if i==0 else 1.8)
        F[span,0]+=[.35,-.2][i]
    e0,S0,W0=K.collision(e,F,E)
    assert np.max(abs(e0))<2e-9 and max(abs(v) for v in S0)<1e-18
    return dict(grid=grid,order=order,mode=mode,angular_floor=floor,entropy=S,
                species_energy_drift=W,separate_temperature_equilibrium_rhs_inf=float(np.max(abs(e0))),
                initial_scalar_source_tangents=tangent,max_remap_rejection=K.reject,
                packets=[h[-1].size for h in K.graph])


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('grid',choices=('fine','third'))
    p.add_argument('--level',type=int,choices=(0,1),default=0)
    p.add_argument('--mode',choices=('both','moller','bhabha'),default='both')
    p.add_argument('--floor',type=float,default=.15)
    a=p.parse_args()
    z=np.linspace(.1,.9,71)
    assert np.max(abs(spin(z)[0]-spin(1-z)[0]))<2e-12 and np.min(spin(z)[1])>0.
    print(json.dumps(audit(**vars(a)),indent=2,allow_nan=False))
