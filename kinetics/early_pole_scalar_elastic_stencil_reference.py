"""Isolated scalar 2->2 reversible, positive energy-projected reference.

Uses a fixed quasiparticle mass and the vacuum constant amplitude. Not a
thermal field theory computation; no claim about the extremely occupied medium.
"""
import argparse
import gzip
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import coo_matrix
from scipy.special import roots_legendre


def remap(q, m, angle_order, screen=None):
    """Sparse directions and positive weights for isotropic 2->2 events."""
    p = q
    E = np.hypot(p,m)
    K = p*p/(E+m)
    x,w=roots_legendre(angle_order)
    i,j = np.triu_indices(len(p))
    # For distinct pairs sum over both orders. The outer factor 1/4
    # corresponds to a symmetric signed weak form over ordered pairs.
    mult = np.where(i==j,.25,.5)
    po,pk = p[i],p[j]
    ea,eb = E[i],E[j]
    center = (K[i]+K[j])/2
    d0 = K[i]*eb+m*K[j]
    storage=[]
    weights=[]
    region=[0.,0.]
    for u,v in zip(x,w,strict=True):
        d=np.maximum(d0-po*pk*u,0.)
        shat=4*m*m+2*d
        lab=np.sqrt(np.maximum(po*po+pk*pk+2*po*pk*u,0.))
        vg=np.sqrt(d*(2*m*m+d))/(ea*eb)
        J=lab*np.sqrt(d/(2*shat))
        for a,b in zip(x,w,strict=True):
            lo=center-J*a
            hi=center+J*a
            A=np.searchsorted(K,lo,side='right')
            B=np.searchsorted(K,hi,side='right')
            W=mult*v*b*vg/(4*shat)
            keep=(A>0)&(B>0)&(A<len(p))&(B<len(p))&(W>0.)
            region[0]+=W.sum()
            region[1]+=W[~keep].sum()
            # Keep escaped geometry in a separate check, not in events.
            k1,k2=A[keep],B[keep]
            er1=(lo[keep]-K[k1-1])/(K[k1]-K[k1-1])
            er2=(hi[keep]-K[k2-1])/(K[k2]-K[k2-1])
            vals=np.c_[np.ones((keep.sum(),2)),1-er1,er1,1-er2,er2]
            where=np.c_[i[keep],j[keep],k1-1,k1,k2-1,k2]
            storage.append((where,vals))
            if screen is None:
                factor=1.
            else:
                factor=screen(po[keep],pk[keep],u,a,d[keep],shat[keep],
                              lab[keep],ea[keep],eb[keep],lo[keep])
                assert factor.shape==W[keep].shape
                assert np.min(factor,initial=0.)>=0. and np.all(np.isfinite(factor))
            weights.append(W[keep]*factor)
    jj=np.concatenate([a[0] for a in storage],axis=0)
    value=np.concatenate([a[1] for a in storage],axis=0)
    W=np.concatenate(weights)
    return jj,value,W,K,region[1]/region[0]


class Elastic:
    def __init__(self,p,m,edges,angular=8,lam=1.,screen=None):
        self.p=np.asarray(p)
        self.m=m
        self.j,self.v,w,self.K,self.raw_escape=remap(self.p,m,angular,screen)
        self.vol=np.diff(np.asarray(edges)**3)/(6*np.pi**2)
        self.weight=w*self.vol[self.j[:,0]]*self.vol[self.j[:,1]]*lam**2/(32*np.pi)
        nr=np.repeat(np.arange(len(w)),6)
        c=self.j.ravel()
        signed=self.v*np.array([-1.,-1.,1.,1.,1.,1.])
        self.S=coo_matrix((signed.ravel(),(nr,c)), shape=(len(w),len(p))).tocsr()
        self.A=coo_matrix((self.v.ravel(),(nr,c)), shape=(len(w),len(p))).tocsr()
        assert np.max(abs(self.S@self.K))<3e-13*max(self.K)
        assert np.max(abs(self.S@np.ones(len(p))))<1e-14

    def reaction(self,f):
        # Infinite entropy gradients at an empty supported cell are
        # handled by the one-way occupation products, not by a floor.
        assert np.min(f)>=0.
        positive=f>0
        lf=np.log(np.where(positive,f,1.))
        h=np.logaddexp(0.,-lf)
        d=self.S@h
        common=self.A@lf
        in_sum=h[self.j[:,0]]+h[self.j[:,1]]
        lo=in_sum+np.minimum(d,0.)
        close=abs(d)<.1
        R=np.zeros(len(d))
        R[close]=np.exp((common+lo)[close])*np.sign(d[close])*np.expm1(abs(d[close]))
        far=~close
        R[far]=np.exp((common+lo+abs(d))[far])*np.sign(d[far])*(-np.expm1(-abs(d[far])))
        if not np.all(positive):
            alive=positive[self.j] | (self.v==0.)
            zin=~alive[:,:2].all(axis=1)
            zout=~alive[:,2:].all(axis=1)
            forward=(~zin)&zout
            rev=zin&(~zout)
            fc,ww=f[self.j],self.v
            left=np.log1p(fc[:,:2]).sum(axis=1)
            right=(ww[:,2:]*np.log1p(fc[:,2:])).sum(axis=1)
            R[forward]=np.exp((common+right)[forward])
            R[rev]=-np.exp((common+left)[rev])
            R[zin&zout]=0.
            d[forward]=np.inf; d[rev]=-np.inf; d[zin&zout]=0.
        return self.weight*R,d

    def flux(self,f):
        R,d=self.reaction(f)
        v=self.S.T@R
        active=R!=0.
        h=float(R[active]@d[active])
        assert h>=0.
        return v,h

    def checks(self,T):
        be=1./np.expm1((self.K+.05*T)/T)
        t=self.flux(be)[0]
        y=be*(1+.8*np.exp(-((self.p-self.p[len(self.p)//3])/self.p[len(self.p)//3])**2))
        s,h=self.flux(y)
        unit=np.sum(abs(s))
        norm=np.sum(abs(t))/unit
        reln=abs(s.sum())/unit
        rele=abs(s@self.K)/sum(abs(s)*self.K)
        assert max(norm,reln,rele)<5e-11 and h>0.
        test=y.copy()
        test[len(test)//3]=0.
        empty,ent=self.flux(test)
        assert empty[len(test)//3]>0 and np.isposinf(ent)
        return dict(thermal_to_perturbed_flux=float(norm),number_rate_residual=float(reln),
                    energy_rate_residual=float(rele),entropy_production=float(h))

    def evolve(self):
        T=max(self.K)*.12
        static=self.checks(T)
        f=1./np.expm1((self.K+.05*T)/T)
        f*=1+.8*np.exp(-((self.p-self.p[len(self.p)//3])/self.p[len(self.p)//3])**2)
        n0=self.vol*f
        rate=np.max(abs(self.flux(f)[0])/n0)
        def rhs(t,y):
            return self.flux(y/self.vol)[0]/rate
        times=np.linspace(0,5,26)
        sol=solve_ivp(rhs,(0.,5.),n0,t_eval=times,method='DOP853',
                      rtol=2e-10,atol=n0*1e-11,max_step=.05)
        assert sol.success and sol.y.min()>0.
        drift_number=np.max(abs(sol.y.sum(axis=0)/(n0.sum())-1))
        drift_energy=np.max(abs((self.K@sol.y)/(self.K@n0)-1))
        z=sol.y/self.vol[:,None]
        # Stable entropy: (1+f) log(1+f)-f log(f).
        H=self.vol@(np.log1p(z)+z*np.log1p(1/z))
        assert min(np.diff(H))>=-1e-13*H[-1]
        assert max(drift_number,drift_energy)<5e-10
        return dict(dynamic_max_number_residual=float(drift_number),
                    dynamic_max_energy_residual=float(drift_energy),
                    entropy_fractional_increase=float(H[-1]/H[0]-1),
                    minimum_sampled_occupation=float(z.min()),**static)


def manufactured(n,angular):
    edges=np.linspace(0.,2.5,n+1)
    p=.75*np.diff(edges**4)/np.diff(edges**3)
    d=Elastic(p,.4,edges,angular)
    a=d.evolve()
    a.update(cells=n,angular_order=angular,projected_events=len(d.weight),
             unweighted_escape_fraction=d.raw_escape)
    return a


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cells',type=int,default=20)
    p.add_argument('--angle',type=int,default=8)
    a=p.parse_args()
    print(json.dumps(manufactured(a.cells,a.angle),indent=2,allow_nan=False))
