"""Finite parent count with imposed thermal shape; paired daughter event closure."""
import json
import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import roots_legendre,expit
from early_pole_paired_hartree_matrix_reference import DirectDriver
from early_pole_decay_mode_clock_reference import X
from early_pole_shell_hartree_feedback_reference import S

p,w=roots_legendre(72)
p0=(p+1)*25.
w0=25*w*p0*p0*4/(2*np.pi**2)

def neq(x):
    return w0@expit(-np.hypot(p0,x))


class Reservoir(DirectDriver):
    def __init__(self,ns,nf,lam,gamma):
        super().__init__(ns,nf,lam,0.,1e-8)
        self.gamma=gamma

    def core(self,x,f,g,m,theta,j=False):
        a,D=self.kernels(x,m)
        A=a*theta
        b=A-D*g
        R=A*(1-g)+f[:,None]*b
        phys=np.r_[R.sum(axis=1),(self.w@R)/self.wf]
        if not j:
            return phys
        t=a*(1+f[:,None]-g[None,:])
        da=np.r_[t.sum(axis=1),(self.w@t)/self.wf]
        h=A+D*f[:,None]
        J=np.block([[np.diag(b.sum(axis=1)),-h],
                    [(self.w[:,None]*b).T/self.wf[:,None],np.diag(-(self.w@h)/self.wf)]])
        n=neq(x)
        J[:,:self.n]-=np.outer(da,self.w)/n
        return phys,J,da/n

    def rhs_R(self,t,q):
        x=X*np.exp(t)
        f,g=q[:self.n],q[self.n:-1]
        N=q[-1]-self.w@f
        n=neq(x)
        return np.r_[self.core(x,f,g,self.gap(x,f),N/n),self.gamma*(n-N)]

    def jac_R(self,t,q):
        x=X*np.exp(t)
        f,g=q[:self.n],q[self.n:-1]
        N=q[-1]-self.w@f
        n=neq(x)
        m=self.gap(x,f)
        _,A,ad=self.core(x,f,g,m,N/n,True)
        if self.lam:
            E=np.hypot(self.p,m)
            la=self.lam/2
            dm=la*self.w/E/(2*m+la*((self.w*np.maximum(f,0.))@(m/E**3)))
            h=max(2e-12,m*2e-5)
            lo=max(1e-11*x,m-h)
            slope=(self.core(x,f,g,m+h,N/n)-self.core(x,f,g,lo,N/n))/(m+h-lo)
            A[:,:self.n]+=np.outer(slope,dm)
        return np.block([[A,ad[:,None]],
                         [np.r_[self.gamma*self.w,np.zeros(len(g))][None,:],np.array([[-self.gamma]])]])

    def run_R(self,xend):
        ic=np.r_[np.zeros(self.n+len(self.wf)),neq(X)]
        ep=np.r_[np.full(self.n,5e-8),np.full(len(self.wf),6e-12),1e-12]
        s=solve_ivp(self.rhs_R,(0.,np.log(xend/X)),ic,method='BDF',rtol=1e-8,
                    atol=ep,jac=self.jac_R,max_step=.024,dense_output=True)
        assert s.success
        rows=[]
        for t in np.linspace(.0005,s.t[-1],29):
            q=s.sol(t);x=X*np.exp(t)
            y=q[:self.n]@self.w
            d=q[self.n:-1]@self.wf
            h=q[-1]-y
            assert h>0 and h/neq(x)<2.
            assert abs(y-d)/max(y,1e-40)<2e-7
            if not self.gamma:
                assert abs(q[-1]/ic[-1]-1)<1e-14
            rows.append(dict(x=x,scalar_Y=y/S,physical_parent_count_over_equilibrium=h/neq(x),
                             dimensionless_external_parent_replenishment=q[-1]-ic[-1]))
        return dict(lambda_self_quartic=self.lam,parent_replenishment_rate_per_H=self.gamma,history=rows)


def run():
    res=[Reservoir(65,90,lam,g).run_R(x) for lam,x in ((7.65175e-8,.05),(0.,.025))
         for g in (0.,3.,100.)]
    return dict(scope="Fixed bath temperature and maintained parent FD momentum shape times a finite count normalization. "
                      "Relaxation replenishes count by prescription, not a derived process; no momentum-resolved parent "
                      "or complete covariant stress-energy evolution.",runs=res)

if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
