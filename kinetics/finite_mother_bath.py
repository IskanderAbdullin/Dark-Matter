"""Kinetic-equilibrium parent with finite chemical relaxation in an external bath.
Tracks parent number, not its momentum distortion or the full thermal energy budget.
"""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.special import expit,roots_legendre
from radiation_fixed_bath import RadiationFixedBath,omega_from_number


class FiniteMotherBath(RadiationFixedBath):
    def __init__(self,n=200,qmin=1e-7,qmax=40.,kappa=.1347,relaxation=1.):
        super().__init__(n,qmin,qmax,kappa)
        self.relaxation=relaxation
        p,w=roots_legendre(160);p=40*(p+1);w=40*w
        self.parent_p=p;self.parent_w=4/(2*np.pi**2)*p*p*w
        self.pair_E=self.mean[:,None]+self.mean[None,:]
        self.last_state=None

    def mother_eq(self,x):
        E=np.hypot(self.parent_p,x)
        return float(self.parent_w@expit(-E))

    def mother_mu(self,x,N):
        E=np.hypot(self.parent_p,x)
        if N<=0:return -np.inf,E
        def residual(mu):return self.parent_w@expit(mu-E)-N
        mu=brentq(residual,x-750,x+80,xtol=1e-12)
        return mu,E

    def set_mother(self,x,N):
        if self.last_state==(x,N):return
        super().update(x)
        mu,E=self.mother_mu(x,N)
        f=expit(mu-E)
        self.Neq=float(self.parent_w@expit(-E));self.mu=mu
        F=expit(mu-self.pair_E)
        self.A=self.R*F;self.asum=self.A.sum(axis=1)
        if N>0:
            susceptibility=float(self.parent_w@(f*(1-f)))
            self.dA=self.R*(F*(1-F))/susceptibility
        elif N==0:
            self.dA=self.R*np.exp(np.minimum(700,x-self.pair_E))/float(self.parent_w@np.exp(x-E))
        else:
            self.dA=np.zeros_like(self.R)
        self.dasum=self.dA.sum(axis=1);self.last_state=(x,N)

    def rhs_log(self,z,y):
        x=np.exp(z);N=y[-1];self.set_mother(x,N)
        dy=x*super(RadiationFixedBath,self).rhs(x,y[:-1])
        source=self.relaxation*(self.Neq-N)
        return np.r_[dy,-self.vol@dy[:self.n]+source]

    def jac_log(self,z,y):
        x=np.exp(z);N=y[-1];self.set_mother(x,N);p,f=np.split(y[:-1],2)
        J=x*super(RadiationFixedBath,self).jac(x,y[:-1]);v=self.vol
        d=x*np.r_[(self.dA@(1-f)+p*self.dasum)/v,(self.dA@(1+p)-f*self.dasum)/(4*v)]
        lastrow=-v@J[:self.n,:]
        lastcol=-v@d[:self.n]-self.relaxation
        return np.block([[J,d[:,None]],[lastrow[None,:],np.array([[lastcol]])]])

    def evolve_mother(self,x0=.1,xend=40.,rtol=2e-6,keep_psd=False):
        init=self.mother_eq(x0);y0=np.r_[np.zeros(2*self.n),init]
        samples=sorted(set([x0]+[x for x in [.03,.05,.075,.1,.15,.2,.3,.5,1,2,3,5,10,15,20,30,40,60] if x0<x<xend]+[xend]))
        sol=solve_ivp(self.rhs_log,[np.log(x0),np.log(xend)],y0,method='BDF',jac=self.jac_log,
                      rtol=rtol,atol=1e-12,max_step=.04,t_eval=np.log(samples))
        if not sol.success:raise RuntimeError(sol.message)
        history=[]
        for x,y in zip(samples[1:],sol.y.T[1:]):
            s=self.snapshot(y[:-1]);s['x']=x;N=float(y[-1]);mu,E=self.mother_mu(x,N)
            p,f=np.split(y[:-1],2);rho1=float(self.parent_w@(E*expit(mu-E)))
            s.update(mother_number=N,mother_eq_number=self.mother_eq(x),mother_mu_over_T=float(mu) if np.isfinite(mu) else None,
                     mother_rho_over_T4=rho1,mother_plus_phi_number=N+s['number'],
                     rho_phi_over_T4=float((p*self.vol)@self.mean),rho_chi2_over_T4=float((4*f*self.vol)@self.mean),
                     omega_phi_h2_m86_g106p75=float(omega_from_number(s['number'])))
            assert N>-1e-9 and abs(s['number_mismatch'])<1e-7 and s['min_occupation']>-1e-7 and s['max_fermion']<1.00001,s
            if self.relaxation==0:assert abs(s['mother_plus_phi_number']/init-1)<1e-7,s
            history.append(s)
        out=dict(n=self.n,x0=x0,xend=xend,kappa=self.kappa,relaxation_per_H=self.relaxation,
                 initial_mother_number=init,nfev=sol.nfev,rtol=rtol,history=history,final=history[-1])
        if keep_psd:
            out['frozen_cell_spectrum']={'edges':self.edges.tolist(),'f_phi':sol.y[:self.n,-1].tolist(),
                'f_chi2':sol.y[self.n:2*self.n,-1].tolist(),'g_phi':1,'g_chi2':4}
        return out


def main():
    p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=180);p.add_argument('--x0',type=float,default=.1)
    p.add_argument('--xend',type=float,default=40.);p.add_argument('--kappa',type=float,default=.1347)
    p.add_argument('--relaxation',type=float,nargs='+',default=[0.,1.,10.,100.]);p.add_argument('--rtol',type=float,default=2e-6)
    p.add_argument('--output',required=True);p.add_argument('--keep-psd',action='store_true');a=p.parse_args();rows=[]
    for lam in a.relaxation:
        m=FiniteMotherBath(a.n,kappa=a.kappa,relaxation=lam)
        # Check the analytic Jacobian on a nonzero physical state.
        x=.2;z=np.log(x);y=np.r_[np.full(m.n,.1),np.full(m.n,.05),m.mother_eq(x)*.8]
        rng=np.random.default_rng(17);d=rng.normal(size=len(y));d/=np.linalg.norm(d)
        eps=1e-6;j=m.jac_log(z,y)@d;fd=(m.rhs_log(z,y+eps*d)-m.rhs_log(z,y-eps*d))/(2*eps)
        assert np.linalg.norm(j-fd)/max(np.linalg.norm(j),1e-30)<1e-5
        r=m.evolve_mother(a.x0,a.xend,a.rtol,a.keep_psd);rows.append(r)
        print(json.dumps({'relaxation':lam,'n':a.n,'nfev':r['nfev'],**r['final']}),flush=True)
    Path(a.output).write_text(json.dumps(rows,indent=2,allow_nan=False)+'\n')

if __name__=='__main__':main()
