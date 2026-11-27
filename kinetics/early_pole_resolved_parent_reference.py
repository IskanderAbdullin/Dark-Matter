"""Early scalar/fermion decay with three separately evolved momentum inventories."""
import argparse
import json

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import expit, roots_legendre

from early_pole_decay_mode_clock_reference import G, KAPPA, X
from early_pole_paired_hartree_matrix_reference import DirectDriver, S


class ResolvedParent:
    """Massless daughter; constant fermion occupancy in each finite cell."""
    def __init__(self, ns, nf, npmid, gamma, held):
        self.d = DirectDriver(ns, nf, 7.65175e-8, 0., 1e-9)
        self.pe = np.r_[0., np.geomspace(.002,.20,13)[:-1],
                        np.linspace(.20,4.,npmid+1),
                        np.linspace(4.,28.,npmid//2+1)[1:],32.]
        Q, W = roots_legendre(12)
        widths = np.diff(self.pe)
        self.pq = self.pe[:-1,None]+widths[:,None]*(Q+1)/2
        self.pqw = G*widths[:,None]*W*self.pq**2/(4*np.pi**2)
        self.pw = G*np.diff(self.pe**3)/(6*np.pi**2)
        self.ns = self.d.n
        self.nf = len(self.d.wf)
        self.iP = self.ns + self.nf
        self.np = len(self.pw)
        self.n = self.iP+self.np
        self.w = np.r_[self.d.w,self.d.wf,self.pw]
        self.gamma = gamma
        self.held = held

    def parent_fd(self,x):
        return np.sum(self.pqw*expit(-np.hypot(self.pq,x)),axis=1)/self.pw

    def overlaps(self,x,m):
        d=self.d
        if x<=m:
            a=np.array([],dtype=int)
            return a,a,a,np.array([],dtype=float)
        E=np.hypot(d.p,m)
        delta=(x*x-m*m)/2
        lo=(x*x*E*E+delta*delta)/(((x*x+m*m)/2)*E+d.p*delta)
        hi=lo+2*d.p*delta/(m*m)
        pe=np.hypot(x,self.pe)
        edges = np.concatenate((np.broadcast_to(pe,(self.ns,len(pe))),
                                E[:,None]+d.ef),axis=1)
        cut=np.sort(edges,axis=1)
        L=np.maximum(cut[:,:-1],lo[:,None])
        H=np.minimum(cut[:,1:],hi[:,None])
        i,k=np.nonzero(H>L)
        l,h=L[i,k],H[i,k]
        c=l+(h-l)/2
        j=np.searchsorted(d.ef,c-E[i],side='right')-1
        k=np.searchsorted(pe,c,side='right')-1
        keep=(j>=0)&(j<self.nf)&(k>=0)&(k<self.np)
        i,j,k,l,h=(z[keep] for z in (i,j,k,l,h))
        factor=KAPPA*G*x*(x*x-m*m)/(E[i]*d.p[i])
        return i,j+self.ns,k+self.iP,(h-l)*factor*d.w[i]

    def at(self,x,m,y,jac=False):
        i,j,k,W=self.overlaps(x,m)
        eq=self.parent_fd(x) if self.held or self.gamma else None
        P=eq[k-self.iP] if self.held else y[k]
        u,b=y[i],y[j]
        R=W*(P*(1+u-b)-u*b)
        v=np.bincount(i,weights=R,minlength=self.n)
        v+=np.bincount(j,weights=R,minlength=self.n)
        if not self.held:
            v-=np.bincount(k,weights=R,minlength=self.n)
        v/=self.w
        if self.gamma and not self.held:
            v[self.iP:]+=self.gamma*(eq-y[self.iP:])
        if not jac:
            return v
        J=np.zeros((self.n,self.n))
        mask = slice(self.iP,None)
        incols = ((i,W*(P-b)),(j,-W*(P+u)))
        if not self.held:
            incols+=((k,W*(1+u-b)),)
        for col, val in incols:
            np.add.at(J,(i,col),val)
            np.add.at(J,(j,col),val)
            if not self.held:
                np.add.at(J,(k,col),-val)
        J/=self.w[:,None]
        if self.gamma and not self.held:
            J[mask,mask]-=self.gamma*np.eye(self.np)
        return J

    def rhs(self,t,y):
        x=X*np.exp(t)
        return self.at(x,self.d.gap(x,y[:self.ns]),y)

    def jac(self,t,y):
        x=X*np.exp(t)
        u=y[:self.ns]
        m=self.d.gap(x,u)
        J=self.at(x,m,y,True)
        E=np.hypot(self.d.p,m)
        a=self.d.lam/2
        dm=a*self.d.w/E/(2*m+a*((self.d.w*np.maximum(u,0.))@(m/E**3)))
        h=max(2e-12,m*2e-5)
        lo=max(1e-11*x,m-h)
        v=(self.at(x,m+h,y)-self.at(x,lo,y))/(m+h-lo)
        J[:,:self.ns]+=np.outer(v,dm)
        return J

    def check(self):
        x=.009
        # Arbitrary positive scalar and interior fermion distributions.
        y=np.r_[np.geomspace(1e3,.0001,self.ns),
                np.linspace(.2,.06,self.nf),self.parent_fd(x)]
        v=self.rhs(np.log(x/X),y)
        rates=[v[:self.ns]@self.w[:self.ns],
               v[self.ns:self.iP]@self.w[self.ns:self.iP],
               v[self.iP:]@self.pw]
        assert abs(rates[0]-rates[1])<1e-14*max(1.,abs(rates[0]))
        if not self.held:
            assert abs(rates[0]+rates[2])<1e-14*max(1.,abs(rates[0]))
            for span in (slice(self.ns,self.iP),slice(self.iP,None)):
                for f in (0.,1.):
                    z=y.copy()
                    z[span]=f
                    assert np.min((1-2*f)*self.rhs(np.log(x/X),z)[span])>=-1e-14
        return rates

    def run(self,xend):
        self.check()
        init=np.r_[np.zeros(self.iP),self.parent_fd(X)]
        atol=np.r_[np.full(self.ns,1e-7),np.full(self.nf,8e-12),np.full(self.np,8e-12)]
        tmax=np.log(xend/X)
        sol=solve_ivp(self.rhs,(0.,tmax),init,method='BDF',jac=self.jac,
                      atol=atol,rtol=2e-8,max_step=.020,dense_output=True)
        assert sol.success,sol.message
        v=sol.sol(np.linspace(0,tmax,51))
        assert np.min(v[:self.ns])>=-4e-7
        assert np.min(v[self.ns:])>=-5e-10 and np.max(v[self.ns:])<=1+5e-10
        pair=self.w[:self.ns]@v[:self.ns]-self.w[self.ns:self.iP]@v[self.ns:self.iP]
        parent=self.pw@(v[self.iP:]-init[self.iP:,None])
        n=self.w[:self.ns]@v[:self.ns]
        assert np.max(abs(pair))<2e-12
        if not self.gamma and not self.held:
            assert np.max(abs(n+parent))<2e-12
        f_end=v[self.iP:,-1] if not self.held else self.parent_fd(xend)
        dep=1-f_end/init[self.iP:]
        k=int(np.argmax(dep))
        return dict(xend=xend,held=self.held,gamma=self.gamma,scalar_cells=self.ns,
                    fermion_daughter_cells=self.nf,parent_cells=self.np,solver_calls=sol.nfev,
                    scalar_yield=float(n[-1]/S),final_parent_fraction=float((self.pw@f_end)/(self.pw@init[self.iP:])),
                    largest_parent_cell_depletion=float(dep[k]),largest_parent_cell_edges_over_T=self.pe[k:k+2].tolist(),
                    parent_cell_edges_over_T=self.pe.tolist(),parent_initial_occupations=init[self.iP:].tolist(),
                    parent_final_occupations=f_end.tolist(),absolute_pair_count_defect_over_T3=float(np.max(abs(pair))),
                    absolute_final_parent_plus_scalar_drift_over_T3=float(n[-1]+parent[-1]))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--scalar',type=int,default=45)
    p.add_argument('--fermion',type=int,default=60)
    p.add_argument('--parent',type=int,default=45)
    p.add_argument('--xend',type=float,default=.05)
    p.add_argument('--gamma',type=float,default=0.)
    p.add_argument('--held',action='store_true')
    p.add_argument('--static-only',action='store_true')
    a=p.parse_args()
    assert a.gamma>=0
    d=ResolvedParent(a.scalar,a.fermion,a.parent,a.gamma,a.held)
    print(json.dumps(d.check() if a.static_only else d.run(a.xend),indent=2,allow_nan=False))
