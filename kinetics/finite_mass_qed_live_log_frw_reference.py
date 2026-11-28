"""Live two-moment FRW closure with one symmetric finite-mass log energy."""
import argparse
import json

import numpy as np
from scipy.special import expit, roots_legendre

from early_pole_live_joint_hard_packets_frw_reference import Expanding


class ExchangeLog:
    """Triangular quadrature of one species' log functional."""
    def __init__(self,cells,n):
        self.n=len(cells)
        x,w=roots_legendre(n)
        t=(x+1)/2
        w=w/2
        a,b=t[:,None]**2,t[None,:]
        W=(2*t*w)[:,None]*w[None,:]
        self.cells=cells
        rows=[]
        def add(i,j,p,r,gap,weight):
            p,r,gap,weight=np.broadcast_arrays(p,r,gap,weight)
            m=cells[i].m
            if m!=cells[j].m or np.any(gap<=0):
                raise ValueError('Incompatible massive log cell or singular distance.')
            klog=np.log(p+r)-np.log(gap)
            rows.append(np.array([np.full(p.size,i),np.full(p.size,j),
                                  p.ravel(),r.ravel(),cells[i].basis(p).ravel(),cells[j].basis(r).ravel(),
                                  (weight*p*r*klog).ravel()]))
        for i,C in enumerate(cells):
            L,U=C.lo,C.hi
            wi=U-L
            p=L+wi*a
            for j,B in enumerate(cells[:i+1]):
                wj=B.hi-B.lo
                if i==j:
                    gap=wi*a*b**4
                    add(i,j,p,p-gap,gap,wi*W*(wi*a)*4*b**3)
                elif i==j+1:
                    if abs(B.hi-L)>1e-14:
                        raise ValueError('The archived cell block is not contiguous.')
                    z=a*b
                    add(i,j,L+wi*a,L-wj*z,a*(wi+wj*b),W*a*wi*wj)
                    add(i,j,L+wi*z,L-wj*a,a*(wi*b+wj),W*a*wi*wj)
                else:
                    r=B.lo+wj*b
                    add(i,j,p,r,p-r,W*wi*wj)
        i,j,self.p,self.q,self.bp,self.bq,self.base=np.concatenate(rows,axis=1)
        self.i=i.astype(int)
        self.j=j.astype(int)

    def evaluate(self,fit,eta):
        i,j=self.i,self.j
        x=self.cells[0].m*np.exp(eta)
        ep,eq=np.hypot(self.p,x),np.hypot(self.q,x)
        K=x*x*self.base/(ep*eq)
        f=expit(fit[i,0]+fit[i,1]*self.bp)
        b=expit(fit[j,0]+fit[j,1]*self.bq)
        shared=K*f*b
        beta=np.zeros((self.n,2))
        for index,deriv,basis in ((i,shared*(1-f),self.bp),(j,shared*(1-b),self.bq)):
            beta[:,0]+=np.bincount(index,weights=deriv,minlength=self.n)
            beta[:,1]+=np.bincount(index,weights=deriv*basis,minlength=self.n)
        D=np.empty_like(beta)
        for ind,C in enumerate(self.cells):
            v=expit(fit[ind]@C.B)
            M=(C.I*(v*(1-v)))@C.B.T
            D[ind]=np.linalg.solve(M.T,beta[ind])
        return float(shared.sum()),D,float(shared@(2-(x/ep)**2-(x/eq)**2))


class LogExpanding(Expanding):
    def __init__(self,n,quad,scale,log_order):
        super().__init__(n,quad,scale)
        self.logs=[(sl,ExchangeLog(self.cells[sl],log_order)) for sl in self.r[1:]]
        self.log_weight=1.
        self.log_energy=0.
        self.log_trace=0.
        self.log_gradient=np.zeros_like(self.init)

    def weights(self,Y,closed):
        E,F,G,R,S=super().weights(Y,closed)
        log=red=0.
        grad=np.zeros_like(G)
        for sl,obj in self.logs:
            r,v,z=obj.evaluate(F[sl],self.eta)
            log+=r
            grad[sl]=v
            red+=z
        c=-self.c*self.log_weight/2
        self.log_energy=c*log
        self.log_trace=c*red
        self.log_gradient=c*grad
        H=G+self.log_gradient
        self.trace+=self.log_trace
        return H if closed else E,F,H,R+self.log_energy,S


def fixed_check(M):
    shape=M.init
    rows=[]
    for eta in (0.,.028):
        def energy(Y,x):
            M.eta=x
            _,_,H,R,_=M.weights(Y,True)
            return R,H,M.trace,M.log_energy,M.log_trace,M.log_gradient.copy()
        R,H,T,EL,TL,LG=energy(shape,eta)
        direction=.08*shape*np.cos(np.arange(shape.size)).reshape(shape.shape)
        h=3e-4
        def dd(v):
            return (v[0]-8*v[1]+8*v[2]-v[3])/(12*h)
        f=[energy(shape+j*h*direction,eta) for j in (-2,-1,1,2)]
        eta_vals=[energy(shape,eta+j*h) for j in (-2,-1,1,2)]
        e1=abs(dd([r[0] for r in f])-np.sum(direction*H))
        e4=abs(dd([r[3] for r in f])-np.sum(direction*LG))
        e2=abs(dd([r[0] for r in eta_vals])-T)
        e3=abs(dd([r[3] for r in eta_vals])-TL)
        if max(e1,e2)>8e-11 or max(e3,e4)>2e-21:
            raise ValueError('Log energy/trace derivatives did not reproduce.')
        rows.append(dict(log_scale=eta,log_energy=EL,log_trace=TL,
                         tested_energy_variation_error=float(e1),tested_trace_error=float(e2),
                         tested_log_trace_error=float(e3),tested_log_variation_error=float(e4)))
    return rows


def run(level):
    quad=[(2,3,10,-3),(3,5,16,-4)][level]
    M=LogExpanding(12,quad,1.,14)
    tested=fixed_check(M)
    M.eta=.028
    _,F,_,_,_=M.weights(M.init,True)
    accuracy=[]
    for sl,obj in M.logs:
        a,x,b=obj.evaluate(F[sl],M.eta)
        aa,xx,bb=ExchangeLog(M.cells[sl],24).evaluate(F[sl],M.eta)
        defect=max(abs(aa/a-1),abs(bb/b-1),np.max(abs(x-xx))/np.max(abs(xx)))
        if defect>2e-6:
            raise ValueError('Log quadrature refinement did not resolve the moment projection.')
        accuracy.append(float(defect))
    results=[]
    for eta,do_log in ((.028,0.),(.028,1.),(0.,1.)):
        M.log_weight=do_log
        M.reject=[0.,0.]
        r=M.go_frw(True,1.,1.,etaend=eta)
        results.append(dict(include_log=bool(do_log),result=r))
    return dict(scope='Live finite reciprocal log plus vacuum packets and prescribed expansion. No matched rates, Friedmann evolution or exterior exchange.',
                level=level,packet_quad=quad,log_quadrature_relative_defects=accuracy,fixed=tested,rows=results)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--level',type=int,choices=(0,1),default=0)
    p=parser.parse_args()
    print(json.dumps(run(p.level),indent=2,allow_nan=False))
