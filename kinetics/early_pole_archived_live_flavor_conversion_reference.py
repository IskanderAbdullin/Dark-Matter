"""Conservative hard conversion of a parent charge pair into a daughter charge pair."""
import argparse
import json

import numpy as np
from scipy.special import expit, roots_legendre

from early_pole_archived_shared_cold_window_reference import Cold
from early_pole_live_joint_hard_pair_packets_reference import gauss_each


class ConversionKernel:
    def __init__(self,M,order=(1,3,4,1),strength=1.):
        if strength<0.:
            raise ValueError("negative conversion strength")
        self.M=M
        self.graph=self.packets(order)
        self.strength=strength
        self.fac=.005**4/(np.sqrt(np.pi**2*106.75/90)*(86/.05)/2.435e18)
        self.reject=0.

    def packets(self,order):
        nr,ni,nc,na=order
        M=self.M
        grids=[np.asarray(M.edge_sets[i]) for i in (1,2)]
        p,wp=gauss_each(grids[0],nr)
        q,wq=p,wp
        z,w=roots_legendre(ni)
        a,b=np.meshgrid(p,q,indexing='ij')
        lim=4*a*b
        l=np.minimum(1.,1/lim)
        u=l[...,None]+(1-l[...,None])*(1+z)/2
        I=np.outer(wp*p*p/np.pi**2,wq*q*q/np.pi**2)[...,None]*(1-l[...,None])*u*w
        ok=I>0
        W=I[ok]
        s=(lim[...,None]*u)[ok]
        a,b=(np.broadcast_to(p[...,None],I.shape)[ok] for p in (a,b))
        E=a+b
        K=np.sqrt(np.maximum(E*E-s,0.))
        A=np.clip(np.divide(a-b,K,out=np.zeros_like(K),where=K>0.),-1.,1.)
        z,w=roots_legendre(na)
        X,U=roots_legendre(nc)
        packs=[]
        for x,v in zip(X,U,strict=True):
            ct=x
            V=(E+K*A*ct)/2
            R=np.maximum(K*np.sqrt((1-ct*ct)*(1-A*A))/2,1e-30)
            cuts=np.c_[(grids[1]-V[:,None])/R[:,None],
                       (E[:,None]-grids[1]-V[:,None])/R[:,None],
                       np.broadcast_to([-1.,1.],(len(E),2))]
            cuts=np.sort(np.arccos(np.clip(cuts,-1.,1.)),axis=1)
            hh=np.diff(cuts,axis=1)/2
            i,j=np.nonzero(hh>1e-14)
            p=V[i,None]+R[i,None]*np.cos((cuts[i,j]+cuts[i,j+1])[:,None]/2+hh[i,j,None]*z)
            q=E[i,None]-p
            f=(W*v*(1+x*x)/(32*np.pi*s))[i,None]*hh[i,j,None]*w/np.pi
            ok=(p>=grids[1][0])&(p<=grids[1][-1])&(q>=grids[1][0])&(q<=grids[1][-1])
            packs.append(np.array([np.broadcast_to(y[i,None],p.shape)[ok] for y in (a,b,K)]
                                  +[p[ok],q[ok],f[ok]]))
        a,b,K,p,q,W=np.concatenate(packs,axis=1)
        tags=[M.r[1]]*2+[M.r[2]]*2
        ix=np.stack([np.searchsorted(g[1:],x)+s.start for x,g,s in
            zip((a,b,p,q),[grids[0]]*2+[grids[1]]*2,tags,strict=True)])
        ma=M.cells[M.r[1].start].m
        mb=M.cells[M.r[2].start].m
        x=np.stack((np.hypot(a,ma),np.hypot(b,ma),np.hypot(p,mb),np.hypot(q,mb)))
        basis=(x-M.shape[ix,0])/M.shape[ix,1]
        assert max(np.ravel(abs(basis)))<=1.+1e-12
        return ix[:3],basis[:3],p,K,W

    def collision(self,e,F,grad):
        M=self.M
        out=np.zeros_like(e)
        ix,u,p,K,W=self.graph
        span=M.r[2]
        ev=e[ix,0]+e[ix,1]*u
        aim=ev[0]+ev[1]-ev[2]
        j=span.start+np.searchsorted(e[span,0]+e[span,1],aim)
        valid=j<span.stop
        j=np.minimum(j,span.stop-1)
        v=(aim-e[j,0])/e[j,1]
        q=M.cells[span.start].pfrombare(M.shape[j,0]+M.shape[j,1]*v)
        ok=valid&(abs(v)<=1.)&(abs(p-q)<=K)&(p+q>=K)
        v=np.where(ok,v,0.)
        self.reject=max(self.reject,float(W[~ok].sum()/W.sum()))
        a,b,c=F[ix,0]+F[ix,1]*u
        d=F[j,0]+F[j,1]*v
        A,B,C,D=map(expit,(a,b,c,d))
        R=self.fac*self.strength*W*ok*(A*B*(1-C)*(1-D)-C*D*(1-A)*(1-B))
        ent=float(R@(a+b-c-d))
        for j,v,r in zip([*ix,j],[*u,v],[-R,-R,R,R],strict=True):
            out[:,0]+=np.bincount(j,weights=r,minlength=len(out))
            out[:,1]+=np.bincount(j,weights=r*v,minlength=len(out))
        res=float(np.sum(grad*out))
        rate=float(out[span,0].sum())
        heat=float(np.sum(grad[span]*out[span]))
        assert abs(heat-R@(ev[0]+ev[1]))<3e-11*max(abs(heat),1.)
        assert abs(rate-2*R.sum())<3e-11*max(abs(rate),1.)
        assert abs(res)<4e-12 and ent>=-1e-13
        assert abs(out[:,0].sum())<4e-12 and not np.any(out[M.r[0]])
        return out,ent,res,rate,heat

    def audit(self):
        e,F,E,*_=self.M.weights(self.M.init,True)
        _,s,w,n,h=self.collision(e,F,E)
        eq=[]
        for diff in (0.,.5):
            for j,span in enumerate(self.M.r[1:]):
                F[span]=-E[span]
                F[span,0]+=.3-j*diff
            _,S,W,N,H=self.collision(e,F,E)
            eq.append(dict(chemical_difference=diff,entropy=S,energy=W,daughter_count=N,daughter_heat=H))
            if diff==0.:
                assert max(abs(N),abs(H))<1e-10
            elif self.strength:
                assert min(N,H,S)>0.
        return dict(initial_daughter_count=n,initial_daughter_heat=h,initial_entropy=s,
                    initial_energy=w,fixture=eq,packets=self.graph[-1].size,max_reject=self.reject)


class LiveConversion(Cold):
    def __init__(self,grid,exposure,order,strength):
        super().__init__(grid,1,0,8,True,exposure=exposure)
        self.conversion=ConversionKernel(self.M,order,strength)

    def initial_number_budget(self):
        M=self.M
        Y=M.init.copy()
        y=np.r_[Y.ravel(),self.ns,0.,0.]
        conv=self.conversion.collision(*M.weights(Y,True)[:3])[0]
        pair=M.pair_rhs(Y,True)[:-2].reshape(Y.shape)
        decay=Cold.rhs(self,y,False)[:self.j0].reshape(Y.shape)
        def counts(z):
            return [float(z[r,0].sum()) for r in M.r]
        select=slice(self.j0,self.j0+self.n)
        base=Cold.rhs(self,y,False)[select]
        response={}
        for name,ray in (('photon_pair',pair),('conversion',conv)):
            v=np.r_[ray.ravel(),np.zeros(self.n+2)]
            response[name]=[((Cold.rhs(self,y+h*v,False)[select]-base)/h).tolist()
                            for h in (1e-11,1e-10,1e-9)]
        return dict(species_order=['photon','parent','daughter'],decay=counts(decay),
                    photon_pair=counts(pair),conversion=counts(conv),total=counts(decay+pair+conv),
                    initial_scalar_source=base.tolist(),initial_scalar_source_derivatives=response)

    def rhs(self,y,pair):
        a=super().rhs(y,pair)
        if self.conversion.strength:
            D,S,W,_,_=self.conversion.collision(*self.M.weights(y[:self.j0].reshape(-1,2),True)[:3])
            a[:self.j0]+=D.ravel()
            a[-2:]+=[S,W]
        return a


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('grid',choices=('fine','third'))
    p.add_argument('--level',type=int,choices=(0,1),default=0)
    p.add_argument('--strength',type=float,default=1.)
    p.add_argument('--exposure',type=float,default=.0003)
    p.add_argument('--evolve',action='store_true')
    a=p.parse_args()
    quad=[(1,3,4,1),(2,4,5,2)][a.level]
    m=LiveConversion(a.grid,a.exposure,quad,a.strength)
    o=dict(grid=a.grid,order=quad,strength=a.strength,**m.conversion.audit())
    o['initial_number_budget']=m.initial_number_budget()
    if a.evolve:
        o['evolution']=m.run()
        o['maximum_dressed_remap_rejection']=m.conversion.reject
    print(json.dumps(o,indent=2,allow_nan=False))
