"""Live joint solver with vacuum Breit-Wheeler hard packets and exact energy remap."""
import argparse
import json

import numpy as np
from scipy.special import expit,roots_legendre

from early_pole_hard_gauge_pair_source_reference import central
from early_pole_live_joint_pair_decay_moments_reference import Joint


def gauss_each(edges,n):
    z,w=roots_legendre(n)
    h=np.diff(edges)[:,None]/2
    p=(edges[:-1]+edges[1:])[:,None]/2+h*z
    return p.ravel(),(h*w).ravel()


class Packets(Joint):
    def __init__(self,n,quad,scale,masses=(.05,.003),cut=.8):
        super().__init__(.005,n,scale)
        r,inc,fo,az=quad
        self.graph=[]
        self.stats=[]
        self.reject=[0.,0.]
        self.aperture=cut
        pg,wp=gauss_each(np.r_[.5,.92,1.12,1.55,2.],r)
        ww=wp*pg*pg/(np.pi*np.pi)
        z,v=roots_legendre(inc)
        cot,cw=roots_legendre(fo)
        ph,pw=roots_legendre(abs(az))
        for idx,m in enumerate(masses,1):
            P1,P2=np.meshgrid(pg,pg,indexing='ij')
            P1=P1.ravel()[:,None]
            P2=P2.ravel()[:,None]
            outer=np.outer(ww,ww).ravel()[:,None]
            lim=4*P1*P2
            left=np.minimum(1.,1/lim)
            u=left+.5*(1-left)*(1+z)
            s=lim*u
            E=P1+P2+np.zeros_like(s)
            P=np.sqrt(np.maximum(E*E-s,0.))
            A=np.divide(P1-P2,P,out=np.zeros_like(P),where=P>0.)
            A=np.clip(A,-1.,1.)
            beta=np.sqrt(1-4*m*m/s)
            qq=np.broadcast_to(P1,s.shape)
            rr=np.broadcast_to(P2,s.shape)
            inc_weight=.5*outer*(1-left)*u*v
            vv=self.r[idx]
            high=np.array([self.cells[k].hi for k in range(vv.start,vv.stop)])
            lo=self.cells[vv.start].lo
            packets=[];ref=0.;keep=0.
            for x,w1 in zip(cut*cot,cut*cw,strict=True):
                vv1=1-(x*beta)**2
                R=4*m*m/s
                ds=beta*(4*(1+R)/vv1-2-4*R*R/(vv1*vv1))/(16*np.pi*s)
                if az>0:
                    for f,w2 in zip(np.pi*(ph+1),pw/2,strict=True):
                        angle=A*x+np.sqrt((1-x*x)*(1-A*A))*np.cos(f)
                        k1=E/2+P*beta*angle/2
                        k2=E-k1
                        p1,p2=np.sqrt(k1*k1-m*m),np.sqrt(k2*k2-m*m)
                        fw=inc_weight*ds*w1*w2
                        mask=(p1>=lo)&(p2>=lo)&(p1<=high[-1])&(p2<=high[-1])&(fw>0.)
                        ref+=sum(fw.ravel())
                        keep+=sum(fw[mask])
                        packets.append(np.array([qq[mask],rr[mask],p1[mask],P[mask],fw[mask]]))
                else:
                    # Exact azimuth splits at both outgoing cell boundaries.
                    middle=E/2+P*beta*A*x/2
                    ampli=P*beta*np.sqrt((1-x*x)*(1-A*A))/2
                    eedge=np.hypot(np.r_[lo,high],m)
                    tol=np.maximum(ampli,1e-30)
                    bc=np.concatenate([(eedge-middle[...,None])/tol[...,None],
                                       (E[...,None]-eedge-middle[...,None])/tol[...,None],
                                       np.broadcast_to(np.array([-1.,1.]),middle.shape+(2,))],axis=-1)
                    angles=np.sort(np.arccos(np.clip(bc,-1.,1.)),axis=-1)
                    lh=np.diff(angles,axis=-1)[...,None]/2
                    t=(angles[...,1:]+angles[...,:-1])[...,None]/2+lh*ph
                    ep=middle[...,None,None]+ampli[...,None,None]*np.cos(t)
                    en=E[...,None,None]-ep
                    p1,p2=np.sqrt(ep*ep-m*m),np.sqrt(en*en-m*m)
                    fw=(inc_weight*ds*w1)[...,None,None]*lh*pw/np.pi
                    mask=(p1>=lo)&(p2>=lo)&(p1<=high[-1])&(p2<=high[-1])&(fw>0.)
                    ref+=np.sum(fw)
                    keep+=np.sum(fw[mask])
                    packets.append(np.array([np.broadcast_to(qq[...,None,None],fw.shape)[mask],
                        np.broadcast_to(rr[...,None,None],fw.shape)[mask],p1[mask],
                        np.broadcast_to(P[...,None,None],fw.shape)[mask],fw[mask]]))
            exact=np.sum(inc_weight*central(s,m,1.,cut))
            assert abs(ref/exact-1)<2e-4
            a,b,p1,P,W=np.concatenate(packets,axis=1)
            ab=np.searchsorted([.92,1.12,1.55,2.],np.vstack([a,b]))
            hs=vv.start+np.searchsorted(high,p1)
            X=np.vstack([a,b,np.hypot(p1,m)])
            tags=np.vstack([ab,hs])
            e0=self.shape[tags,0]
            e1=self.shape[tags,1]
            uu=(X-e0)/e1
            assert np.max(abs(uu))<=1.+1e-12
            self.graph.append((tags,uu,p1,P,W,vv,m))
            self.stats.append(dict(species=idx,packets=W.size,unblocked_geometry_weight=keep,
                                   compared_full_particle_shape=exact,fraction_retained=keep/ref))
        # Normalize to a common arbitrary clock; not a cosmological g^4/H.
        self.norm=sum(r['unblocked_geometry_weight'] for r in self.stats)

    def pair_rhs(self,Y,closed):
        E,F,G,_,_=self.weights(Y,closed)
        out=np.zeros_like(Y);entropy=0.
        for ii,(tag,u,p1,P,W,span,m) in enumerate(self.graph):
            val=E[tag,0]+E[tag,1]*u
            Nv=val[0]+val[1]-val[2]
            hi=E[span,0]+E[span,1]
            j=span.start+np.searchsorted(hi,Nv)
            inside=j<span.stop
            j=np.clip(j,span.start,span.stop-1)
            q=(Nv-E[j,0])/E[j,1]
            hypot=self.shape[j,0]+self.shape[j,1]*q
            p2=np.sqrt(np.maximum((hypot-m)*(hypot+m),0.))
            mask=inside&(np.abs(q)<=1.)&(np.abs(p1-p2)<=P)&(P<=p1+p2)
            self.reject[ii]=max(self.reject[ii],W[~mask].sum()/W.sum())
            q=np.where(mask,q,0.)
            X=F[tag,0]+F[tag,1]*u
            Z=F[j,0]+F[j,1]*q
            g1,g2=1/np.expm1(X[:2])
            f1,f2=expit(X[2]),expit(Z)
            fw=g1*g2*(1-f1)*(1-f2)
            bw=(1+g1)*(1+g2)*f1*f2
            R=(W/self.norm)*(fw-bw)*mask
            entropy+=R@(-X.sum(axis=0)-Z)
            for ind,basis,A in [(tag[0],u[0],-R),(tag[1],u[1],-R),(tag[2],u[2],R),(j,q,R)]:
                out[:,0]+=np.bincount(ind,weights=A,minlength=len(out))
                out[:,1]+=np.bincount(ind,weights=A*basis,minlength=len(out))
        work=np.sum(G*out)
        assert entropy>=-1e-14 and abs(out[:,0].sum())<2e-13
        if closed: assert abs(work)<3e-13
        return np.r_[out.ravel(),entropy,work]


def run(level,split=False):
    quadrature=([(2,3,10,-3),(3,5,16,-4),(4,7,23,-4)] if split else
                [(2,3,10,10),(3,4,14,14),(4,6,18,24)])[level]
    rows=[]
    for fac in (0.,1.,4.):
        M=Packets(12,quadrature,fac)
        a=M.go_joint(True,1.,duration=.14)
        rows.append(dict(result=a,packet_stats=M.stats,
                         maximum_remap_rejection=M.reject))
    return dict(scope='Vacuum tree event weights on a strict hard finite subset; exact variational '
                      'energy remap and triangle test; no in-medium matrix element or physical time.',
                quad=quadrature,rows=rows)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--level',type=int,default=0)
    p.add_argument('--split-azimuth',action='store_true')
    a=p.parse_args()
    print(json.dumps(run(a.level,a.split_azimuth),indent=2,allow_nan=False))
