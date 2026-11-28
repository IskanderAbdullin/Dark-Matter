"""Conserving live hard Compton projection with Bose/Pauli factors."""
import argparse
import json


import numpy as np
from scipy.special import expit,roots_legendre
from early_pole_archived_shared_cold_window_reference import Cold
from early_pole_live_joint_hard_pair_packets_reference import gauss_each

class ComptonKernel:
    def __init__(self,M,order=(1,3,4,1),strength=1.):
        if strength<0.:
            raise ValueError("negative Compton rate")
        self.M=M
        self.graph=[self.packets(k,order) for k in (1,2)]
        self.strength=strength
        self.fac=.005**4/(np.sqrt(np.pi**2*106.75/90)*(86/.05)/2.435e18)
        self.reject=[0.,0.]

    def packets(self,k,order):
        nr,ni,nc,na=order
        M=self.M
        grids=[np.asarray(M.edge_sets[i]) for i in (0,k)]
        (p,wp),(q,wq)=[gauss_each(e,nr) for e in grids]
        z,w=roots_legendre(ni)
        a,b=np.meshgrid(p,q,indexing='ij')
        lim=4*a*b
        l=np.minimum(1.,1/lim)
        u=l[...,None]+(1-l[...,None])*(1+z)/2
        I=np.outer(wp*p*p/np.pi**2,2*wq*q*q/np.pi**2)[...,None]*(1-l[...,None])*u*w/2
        ok=I>0
        W=I[ok]
        s=(lim[...,None]*u)[ok]
        a,b=(np.broadcast_to(p[...,None],I.shape)[ok] for p in (a,b))
        E=a+b
        K=np.sqrt(np.maximum(E*E-s,0.))
        A=np.clip(np.divide(a-b,K,out=np.zeros_like(K),where=K>0.),-1.,1.)
        z,w=roots_legendre(na)
        X,U=gauss_each(np.array([.15,.25,.4,.7,1.]),nc)
        packs=[]
        for x,v in zip(X,U,strict=True):
            ct=2*x-1
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
            f=(W*(2*v)*(1/x+x)/(16*np.pi*s))[i,None]*hh[i,j,None]*w/np.pi
            ok=(p>=grids[0][0])&(p<=grids[0][-1])&(q>=grids[1][0])&(q<=grids[1][-1])
            packs.append(np.array([np.broadcast_to(y[i,None],p.shape)[ok] for y in (a,b,K)]
                                  +[p[ok],q[ok],f[ok]]))
        a,b,K,p,q,W=np.concatenate(packs,axis=1)
        ix=np.stack([np.searchsorted(g[1:],x)+s.start for x,g,s in
            zip((a,b,p,q),grids*2,[M.r[0],M.r[k]]*2,strict=True)])
        x=np.stack((a,np.hypot(b,M.cells[M.r[k].start].m),p,np.hypot(q,M.cells[M.r[k].start].m)))
        basis=(x-M.shape[ix,0])/M.shape[ix,1]
        assert max(np.ravel(abs(basis)))<=1.+1e-12
        return ix[:3],basis[:3],p,K,W

    def collision(self,e,F,grad):
        M=self.M
        out=np.zeros_like(e)
        ent=0.; heat=[]
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
            A,C=1/np.expm1(a),1/np.expm1(c)
            B,D=expit(b),expit(d)
            R=self.fac*self.strength*W*ok*(A*B*(1+C)*(1-D)-C*D*(1+A)*(1-B))
            ent+=R@(-a+b+c-d)
            heat.append(float(R@(aim-ev[1])))
            for j,v,r in zip([*ix,j],[*u,v],[-R,-R,R,R],strict=True):
                out[:,0]+=np.bincount(j,weights=r,minlength=len(out))
                out[:,1]+=np.bincount(j,weights=r*v,minlength=len(out))
        res=float(np.sum(grad*out))
        assert abs(res)<5e-12 and ent>=-1e-13
        assert max(abs(out[r,0].sum()) for r in M.r)<5e-12
        return out,float(ent),res,heat

    def audit(self):
        e,F,E,*_=self.M.weights(self.M.init,True)
        _,S,W,H=self.collision(e,F,E)
        for k,span in enumerate(self.M.r):
            F[span]=(1 if k==0 else -1)*E[span]
            F[span,0]+=[.1,.35,-.2][k]
        _,Q,J,K=self.collision(e,F,E)
        assert max(abs(z) for z in K)<3e-10 and abs(J)<2e-12
        return dict(heat=H,entropy=S,energy=W,common_T_heat=K,common_T_entropy=Q,
                    initial_max_rejected_weight=self.reject.copy(),
                    graph_sizes=[x[-1].size for x in self.graph])


class LiveCompton(Cold):
    def __init__(self,grid,wide=1,n=8,exposure=2e-6,order=(1,3,4,1),strength=1.):
        super().__init__(grid,wide,0,n,True,exposure=exposure)
        self.compton=ComptonKernel(self.M,order,strength)

    def rhs(self,y,pair):
        a=super().rhs(y,pair)
        if self.compton.strength:
            Y=y[:self.j0].reshape(-1,2)
            D,S,W,_=self.compton.collision(*self.M.weights(Y,True)[:3])
            a[:self.j0]+=D.ravel()
            a[-2:]+=[S,W]
        return a


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('grid',choices=('fine','third'))
    p.add_argument('--evolve',action='store_true')
    p.add_argument('--level',type=int,choices=(0,1),default=0)
    p.add_argument('--strength',type=float,default=1.)
    p.add_argument('--exposure',type=float,default=.0003)
    a=p.parse_args()
    order=[(1,3,4,1),(2,4,5,2)][a.level]
    M=LiveCompton(a.grid,exposure=a.exposure,strength=a.strength,order=order)
    raw=dict(grid=a.grid,order=order,strength=a.strength,**M.compton.audit())
    if a.evolve:
        raw['evolution']=M.run()
        raw['maximum_remap_rejection']=M.compton.reject
    print(json.dumps(raw,indent=2,allow_nan=False))
