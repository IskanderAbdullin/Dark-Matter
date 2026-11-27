"""Shared scalar/fermion collision counts for early decay with positive Hartree feedback."""
import argparse
import json
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import solve_ivp
from scipy.special import expit
from early_pole_decay_mode_clock_reference import G,X,KAPPA,kernel
from early_pole_shell_hartree_feedback_reference import S,Y


class Driver:
    def __init__(self,ns,nf,lam,b,qlow):
        self.ed=np.unique(np.r_[np.geomspace(qlow,25.,ns+1),1e-5,.1])
        self.p=np.sqrt(self.ed[:-1]*self.ed[1:])
        self.w=np.diff(self.ed**3)/(6*np.pi**2)
        self.fed=np.r_[np.linspace(0.,4.,nf+1),np.linspace(4.,26.,nf//4+1)[1:]]
        self.ef=np.hypot(self.fed,b)
        self.wf=4*np.diff(self.fed**3)/(6*np.pi**2)
        self.lam=lam
        self.b=b
        self.n=len(self.w)
        self.with_inverse=True

    def gap(self,x,u):
        base=1e-9*x
        n=self.w*np.expm1(u)
        nn=self.lam/2*np.maximum(n,0.)
        if not nn.any():
            return base
        return brentq(lambda m:m*m-base*base-nn@(1/np.hypot(self.p,m)),
                      base,base+1.0001*nn.sum()**(1/3)+1e-20,
                      xtol=1e-18,rtol=1e-11)

    def kernels(self,x,m):
        p=self.p
        if x<=m+self.b:
            z=np.zeros((len(p),len(self.wf)))
            return z,z
        A=(x*x+m*m-self.b**2)/2
        D=np.sqrt((x*x-(m+self.b)**2)*(x*x-(m-self.b)**2))/2
        E=np.hypot(m,p)
        fp=A*E+p*D
        lo=(x*x*E*E+D*D)/fp
        span=2*p*D/(m*m)
        hi=lo+span
        ll=np.maximum(lo[:,None],E[:,None]+self.ef[None,:-1])
        hh=np.minimum(hi[:,None],E[:,None]+self.ef[None,1:])
        # Span can be smaller than the precision of adding endpoints.
        de=np.minimum(np.maximum(hh-ll,0.),span[:,None])
        occ=np.log1p(np.exp(-ll)*(-np.expm1(-de))/(1+np.exp(-ll-de)))
        fac=KAPPA*G*x*((x+self.b)**2-m*m)/(E*p)
        return occ*fac[:,None],de*fac[:,None]

    def rhs_at_m(self,x,u,f,m,jac=False):
        lag=np.exp(-u)
        W=self.w*np.exp(u)
        A,D=self.kernels(x,m)
        net=(A-D*f) + (D-A)*f*lag[:,None]
        f_dot=(W@net)/self.wf
        u_dot=net.sum(axis=1)
        if not jac:
            return np.r_[u_dot,f_dot]
        F=A+(D-A)*lag[:,None]
        a=np.diag(-lag*((D-A)@f))
        b=(W[:,None]*(A-D*f)).T/self.wf[:,None]
        c=np.diag(-(W@F)/self.wf)
        return np.r_[u_dot,f_dot],np.block([[a,-F],[b,c]])

    def rhs(self,t,y):
        x=X*np.exp(t)
        u,f=y[:self.n],y[self.n:]
        return self.rhs_at_m(x,u,f,self.gap(x,u))

    def jac(self,t,y):
        x=X*np.exp(t)
        u,f=y[:self.n],y[self.n:]
        m=self.gap(x,u)
        V,J=self.rhs_at_m(x,u,f,m,True)
        if self.lam:
            dn=self.w*np.exp(u)
            E=np.hypot(self.p,m)
            n=np.maximum(self.w*np.expm1(u),0.)
            a=self.lam/2
            dm=a*dn/E/(2*m+a*(n@(m/E**3)))
            h=max(2e-12,m*2e-5)
            low=max(1e-11*x,m-h)
            v=(self.rhs_at_m(x,u,f,m+h)-self.rhs_at_m(x,u,f,low))/(m+h-low)
            J[:,:self.n]+=np.outer(v,dm)
        return J

    def matrix_check(self):
        x=.006;m=5e-5
        A,D=self.kernels(x,m)
        test=np.array([kernel(x,m,self.b,p)['seed_over_width']*
            kernel(x,m,self.b,p)['vacuum_matched_yukawa_width_ratio']*KAPPA*x*x for p in self.p])
        keep=test>test.max()*1e-7
        err=float(np.max(abs(A.sum(axis=1)[keep]/test[keep]-1)))
        assert err<3e-8,(err,test.max())
        return err

    def run(self,xend):
        init=np.zeros(self.n+len(self.wf))
        eps=np.r_[np.full(self.n,5e-10),np.full(len(self.wf),4e-12)]
        sol=solve_ivp(self.rhs,(0.,np.log(xend/X)),init,method='BDF',jac=self.jac,
                      atol=eps,rtol=1.3e-8,max_step=.024,dense_output=True)
        assert sol.success,sol.message
        def row(t):
            vv=sol.sol(t)
            x=X*np.exp(t)
            k=vv[self.n:];u=vv[:self.n]
            n=self.w*np.expm1(u)
            nf=self.wf*k
            tot=n.sum()
            diff=(tot-nf.sum())/max(tot,1.e-40)
            assert np.min(u)>-2e-9 and np.min(k)>-2e-9 and np.max(k)<1+1e-9
            return dict(x=x,Y=tot/S,Y_sub_1e_5=n[self.ed[1:]<=1e-5].sum()/S,
                        scalar_pole_over_T=self.gap(x,u),max_fermion_occupation=max(k),
                        peak_scalar_occupation=np.expm1(max(u)),pair_number_relative_mismatch=diff,
                        fermion_mean_p_over_T=(nf@((self.fed[1:]+self.fed[:-1])/2))/nf.sum())
        t=np.linspace(.0005,sol.t[-1],31)
        rr=[row(z) for z in t]
        defect=max(abs(z['pair_number_relative_mismatch']) for z in rr)
        assert defect<1e-5,(self.lam,self.b,defect)
        return dict(lambda_self_quartic=self.lam,fermion_pole_over_T=self.b,
                    scalar_bins=self.n,fermion_bins=len(self.wf),solver_calls=sol.nfev,
                    empty_daughter_row_sum_control=self.matrix_check(),
                    largest_sampled_pair_number_defect=defect,checkpoints=rr,
                    scalar_p_edges=self.ed.tolist(),fermion_p_edges=self.fed.tolist(),
                    final_scalar_n_per_visible_entropy=(self.w*np.expm1(sol.y[:self.n,-1])/S).tolist(),
                    final_daughter_per_spin_and_charge_occupations=sol.y[self.n:,-1].tolist())


class DirectDriver(Driver):
    """Linear occupation variables preserve the shared number invariant under BDF."""
    def rhs_at_m(self,x,f1,f2,m,jac=False):
        A,D=self.kernels(x,m)
        base=A*(1-f2)
        if not self.with_inverse:
            D=A
        B=A-D*f2
        R=base+f1[:,None]*B
        y=np.r_[R.sum(axis=1),(self.w@R)/self.wf]
        if not jac:
            return y
        M=A+D*f1[:,None]
        return y,np.block([[np.diag(B.sum(axis=1)),-M],
                           [(self.w[:,None]*B).T/self.wf[:,None],np.diag(-(self.w@M)/self.wf)]])

    def gap(self,x,f):
        return super().gap(x,np.log1p(np.maximum(f,0.)))

    def jac(self,t,y):
        x=X*np.exp(t)
        f,g=y[:self.n],y[self.n:]
        m=self.gap(x,f)
        _,J=self.rhs_at_m(x,f,g,m,True)
        if self.lam:
            E=np.hypot(self.p,m)
            a=self.lam/2
            dm=a*self.w/E/(2*m+a*((self.w*np.maximum(f,0.))@(m/E**3)))
            h=max(2e-12,m*2e-5)
            low=max(1e-11*x,m-h)
            v=(self.rhs_at_m(x,f,g,m+h)-self.rhs_at_m(x,f,g,low))/(m+h-low)
            J[:,:self.n]+=np.outer(v,dm)
        return J

    def run(self,xend):
        ep=np.r_[np.full(self.n,5e-8),np.full(len(self.wf),6e-12)]
        sol=solve_ivp(self.rhs,(0.,np.log(xend/X)),np.zeros(len(ep)),method='BDF',jac=self.jac,
                      atol=ep,rtol=1e-8,max_step=.02,dense_output=True)
        assert sol.success,sol.message
        rows=[]
        for t in np.linspace(.0005,sol.t[-1],31):
            q=sol.sol(t)
            x=X*np.exp(t)
            f,g=q[:self.n],q[self.n:]
            n=self.w*f
            nf=self.wf*g
            tot=n.sum()
            diff=(tot-nf.sum())/max(tot,1e-40)
            assert min(f)>-4e-7 and min(g)>-5e-10 and max(g)<1+2e-9
            rows.append(dict(x=x,Y=tot/S,Y_sub_1e_5=n[self.ed[1:]<=1e-5].sum()/S,
                             scalar_pole_over_T=self.gap(x,f),max_fermion_occupation=max(g),
                             peak_scalar_occupation=max(f),pair_number_relative_mismatch=diff,
                             fermion_mean_p_over_T=(nf@((self.fed[1:]+self.fed[:-1])/2))/nf.sum()))
        defect=max(abs(z['pair_number_relative_mismatch']) for z in rows)
        assert defect<2e-7,(self.lam,self.b,defect)
        return dict(variables="direct",with_inverse=self.with_inverse,lambda_self_quartic=self.lam,fermion_pole_over_T=self.b,
                    scalar_bins=self.n,fermion_bins=len(self.wf),solver_calls=sol.nfev,
                    empty_daughter_row_sum_control=self.matrix_check(),
                    largest_sampled_pair_number_defect=defect,checkpoints=rows,
                    scalar_p_edges=self.ed.tolist(),fermion_p_edges=self.fed.tolist(),
                    final_scalar_n_per_visible_entropy=(self.w*sol.y[:self.n,-1]/S).tolist(),
                    final_daughter_per_spin_and_charge_occupations=sol.y[self.n:,-1].tolist())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--scalar',type=int,default=65);p.add_argument('--fermion',type=int,default=90)
    p.add_argument('--lambda',type=float,dest='lam',default=7.65175e-8)
    p.add_argument('--b',type=float,default=0.);p.add_argument('--xend',type=float,default=.05)
    p.add_argument('--qlow',type=float,default=1e-8)
    p.add_argument('--direct',action='store_true')
    p.add_argument('--no-inverse',action='store_true')
    a=p.parse_args()
    assert a.direct or not a.no_inverse
    D=DirectDriver if a.direct else Driver
    d=D(a.scalar,a.fermion,a.lam,a.b,a.qlow)
    d.with_inverse=not a.no_inverse
    print(json.dumps(d.run(a.xend),indent=2,allow_nan=False))
