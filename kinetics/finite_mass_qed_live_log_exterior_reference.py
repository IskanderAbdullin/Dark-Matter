"""Evolve live archived moments with a frozen hard-exterior log gradient."""
import gzip
import json
from pathlib import Path

import numpy as np
from scipy.special import expit, roots_legendre

from finite_mass_qed_live_log_frw_reference import LogExpanding


def clip_cells(edges, f, bottom, top):
    low=np.maximum(edges[:-1],bottom)
    high=np.minimum(edges[1:],top)
    ok=high>low
    return zip(low[ok],high[ok],f[ok],strict=True)


class ExteriorCross:
    def __init__(self,cells,edges,occ,cut,n=10):
        self.cells=cells
        x,w=roots_legendre(n)
        t=(x+1)/2
        w=w/2
        a,b=t[:,None]**2,t[None,:]
        W=(2*t*w)[:,None]*w[None,:]
        outside=[*clip_cells(edges,occ,cut,cells[0].lo),
                 *clip_cells(edges,occ,cells[-1].hi,np.inf)]
        pieces=[]
        free_p=[]
        free_w=[]
        for l,h,f in outside:
            q=l+(h-l)*t
            free_p.extend(q)
            free_w.extend((h-l)*w*q*q*f)
        self.fixed_q=np.asarray(free_p)
        self.fixed_w=np.asarray(free_w)
        for i,C in enumerate(cells):
            L,U=C.lo,C.hi
            wi=U-L
            def add(p,r,dist,wt):
                p,r,dist,wt=np.broadcast_arrays(p,r,dist,wt)
                v=wt*p*r*(np.log(p+r)-np.log(dist))
                pieces.append(np.array([np.repeat(i,p.size),p.ravel(),r.ravel(),C.basis(p).ravel(),v.ravel()]))
            for l,h,f in outside:
                wr=h-l
                if h==L:
                    add(L+wi*a,L-wr*a*b,a*(wi+wr*b),W*a*wi*wr*f)
                    add(L+wi*a*b,L-wr*a,a*(wi*b+wr),W*a*wi*wr*f)
                elif l==U:
                    add(U-wi*a,U+wr*a*b,a*(wi+wr*b),W*a*wi*wr*f)
                    add(U-wi*a*b,U+wr*a,a*(wi*b+wr),W*a*wi*wr*f)
                else:
                    p=L+wi*a
                    r=l+wr*b
                    add(p,r,abs(p-r),W*wi*wr*f)
        i,self.p,self.q,self.bp,self.pool=np.concatenate(pieces,axis=1)
        self.i=i.astype(int)

    def moments(self,eta):
        x=self.cells[0].m*np.exp(eta)
        E=np.hypot(self.fixed_q,x)
        J=self.fixed_w@(1/E)
        Jd=-self.fixed_w@(x*x/E**3)
        rho=(self.fixed_w@E)*2/np.pi**2
        trace=(self.fixed_w@(x*x/E))*2/np.pi**2
        return J,Jd,rho,trace

    def evaluate(self,fit,eta):
        x=self.cells[0].m*np.exp(eta)
        v=expit(fit[self.i,0]+fit[self.i,1]*self.bp)
        ep,eq=np.hypot(self.p,x),np.hypot(self.q,x)
        joint=x*x*self.pool*v/(ep*eq)
        ds=joint*(1-v)
        n=len(self.cells)
        logits=np.array([np.bincount(self.i,weights=ds,minlength=n),
                         np.bincount(self.i,weights=ds*self.bp,minlength=n)]).T
        H=np.empty_like(logits)
        for i,C in enumerate(self.cells):
            f=expit(fit[i]@C.B)
            m=(C.I*(f*(1-f)))@C.B.T
            H[i]=np.linalg.solve(m.T,logits[i])
        return joint.sum(),H,joint@(2-(x/ep)**2-(x/eq)**2)


class WithExterior(LogExpanding):
    def __init__(self,data,charge_floor,log_order=14):
        super().__init__(12,(2,3,10,-3),1.,14)
        self.cross=[]
        self.cross_scale=0.
        self.cross_energy=0.
        self.cross_trace=0.
        self.cross_gradient=np.zeros_like(self.init)
        z,w=roots_legendre(14)
        for name,span in zip(('parent','fermion'),self.r[1:],strict=True):
            edges=np.asarray(data[f'endpoint_{name}_momentum_edges'])
            f=np.asarray(data[f'endpoint_{name}_occupations'])
            self.cross.append((span,ExteriorCross(self.cells[span],edges,f,charge_floor,log_order)))
            for k,C in enumerate(self.cells[span],span.start):
                moments=np.zeros(2)
                for l,h,fp in clip_cells(edges,f,C.lo,C.hi):
                    p=(l+h)/2+(h-l)*z/2
                    y=w*(h-l)/2*(2*p*p/np.pi**2)*fp
                    moments+=np.array([sum(y),y@C.basis(p)])
                self.init[k]=moments

    def spectators(self,eta):
        J=[0.];D=[0.];R=0.;T=0.
        for _,v in self.cross:
            j,d,r,t=v.moments(eta)
            J.append(j);D.append(d);R+=r;T+=t
        return np.asarray(J),np.asarray(D),R,T

    def weights(self,Y,closed):
        E,F,G,R,S=super().weights(Y,closed)
        c=-self.c*self.cross_scale/2
        R2,T2=0.,0.
        G2=np.zeros_like(G)
        for sl,obj in self.cross:
            r,g,t=obj.evaluate(F[sl],self.eta)
            R2+=r;T2+=t;G2[sl]=g
        self.cross_energy=c*R2
        self.cross_trace=c*T2
        self.cross_gradient=c*G2
        self.trace+=self.cross_trace
        H=G+self.cross_gradient
        return H if closed else E,F,H,R+self.cross_energy,S


def fixed_check(system):
    system.cross_scale=1.
    def lm(Y,eta):
        system.eta=eta
        system.weights(Y,True)
        return system.cross_energy,system.cross_trace,system.cross_gradient.copy()
    init=system.init
    q=.006
    g=lm(init,q)
    da=3e-4
    direction=.1*init*np.sin(np.arange(init.size)).reshape(init.shape)
    def ddy(f):
        a,b,c,d=[f(x) for x in (-2.,-1.,1.,2.)]
        return (a-8*b+8*c-d)/(12*da)
    e1=abs(ddy(lambda x:lm(init,q+x*da)[0])-g[1])
    e2=abs(ddy(lambda x:lm(init+x*da*direction,q)[0])-np.sum(g[2]*direction))
    if max(e1,e2)>1e-20:
        raise ValueError('Cross log first variation or scale trace failed.')
    return [float(e1),float(e2)]


def run():
    rows=[]
    for g in ('fine','third'):
        name=Path(__file__).parent/f'results/early-pole-spectrum-{g}-m06-2026-11-26.json.gz'
        with gzip.open(name) as h:d=json.load(h)
        M=WithExterior(d,.25)
        derivative_error=fixed_check(M)
        M.eta=.006
        _,F,_,_,_=M.weights(M.init,True)
        check=WithExterior(d,.25,24)
        accuracy=[]
        for (sl,C),(_,D) in zip(M.cross,check.cross,strict=True):
            a,u,b=C.evaluate(F[sl],.006)
            aa,uu,bb=D.evaluate(F[sl],.006)
            defect=max(abs(aa/a-1),abs(bb/b-1),np.max(abs(u-uu))/np.max(abs(uu)))
            if defect>1e-6:
                raise ValueError(f'Mixed quadrature did not resolve: {defect}, {aa/a-1}, {bb/b-1}, {np.max(abs(u-uu))/np.max(abs(uu))}')
            accuracy.append(float(defect))
        fixed=dict(self_log=M.log_energy,cross_log=M.cross_energy,
                   self_log_trace=M.log_trace,cross_log_trace=M.cross_trace)
        sc=[]
        for inner,outer in ((0.,0.),(1.,0.),(1.,1.)):
            M.log_weight=inner
            M.cross_scale=outer
            M.reject=[0.,0.]
            x=M.go_frw(True,1.,1.,duration=.04,etaend=.012)
            sc.append(x)
        rows.append(dict(grid=g,diagnostic=fixed,energy_trace_gradient_defects=derivative_error,
                         cross_quadrature_relative_defects=accuracy,factorized=sc[0],
                         local_log_only=sc[1],including_cross=sc[2]))
    return dict(scope='Tagged archive active moments, prescribed expansion; frozen outside hard charged occupation, excluding pure exterior log from the active ledger.',
                rows=rows)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
