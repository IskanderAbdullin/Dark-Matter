"""Radiation-era external thermal mother; massless daughter benchmark.
Constant g*, g*s, T proportional to 1/a. No mother depletion or bath energy evolution.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import expit
from scipy.constants import G, parsec, c, elementary_charge, Boltzmann, hbar
from daughter_energy_fv import DaughterEnergyFV

class RadiationFixedBath(DaughterEnergyFV):
    def __init__(self,n=200,qmin=1e-7,qmax=40.,kappa=.01):
        # Same volume means as the flat model, but the pair areas change with x.
        self.n=n; self.edges=np.geomspace(qmin,qmax,n+1); self.kappa=kappa
        l,h=self.edges[:-1],self.edges[1:]
        self.vol=(h**3-l**3)/(6*np.pi**2)
        self.mean=.75*(h**4-l**4)/(h**3-l**3)
        self.mean2=.6*(h**5-l**5)/(h**3-l**3)
        self.lq,self.hq=l[:,None],h[:,None]
        self.lr,self.hr=l[None,:],h[None,:]
        self.fparent=expit(-self.mean[:,None]-self.mean[None,:])
        self.last_x=None

    def update(self,x):
        if x==self.last_x:return
        c=x*x/4
        lq,hq,lr,hr=self.lq,self.hq,self.lr,self.hr
        # Ramp where the lower r limit is c/q; then the full r interval.
        lo=np.clip(c/hr,lq,hq)
        hi=np.maximum(lo,np.clip(c/lr,lq,hq))
        dq=hi-lo
        area=np.maximum(hr*dq-c*np.log1p(dq/lo),0.)+(hq-hi)*(hr-lr)
        # The kernel is per unit Gamma*t, then multiply by Gamma/(x H).
        self.R=(4*x/(2*np.pi**2))*(area+area.T)/2*(self.kappa*x)
        self.A=self.R*self.fparent
        self.asum=self.A.sum(axis=1);self.last_x=x

    def rhs(self,x,y):
        self.update(x);return super().rhs(x,y)

    def jac(self,x,y):
        self.update(x);return super().jac(x,y)

    def evolve(self,x0=.1,xend=30.,rtol=2e-6,keep_psd=False):
        # Integrate ln x to keep comparable relative step at early and late times.
        def rhs(z,y):
            x=np.exp(z);return x*self.rhs(x,y)
        def jac(z,y):
            x=np.exp(z);return x*self.jac(x,y)
        samples=sorted(set([x0]+[x for x in [.03,.05,.075,.1,.15,.2,.3,.5,1,2,3,5,10,15,20,30,40,60] if x0<x<xend]+[xend]))
        sol=solve_ivp(rhs,[np.log(x0),np.log(xend)],np.zeros(2*self.n),method='BDF',jac=jac,
                      rtol=rtol,atol=1e-12,max_step=.04,t_eval=np.log(samples))
        if not sol.success:raise RuntimeError(sol.message)
        history=[]
        for x,y in zip(samples[1:],sol.y.T[1:]):
            s=self.snapshot(y);s['x']=x
            p,f=np.split(y,2)
            s['rho_phi_over_T4']=float((p*self.vol)@self.mean)
            s['rho_chi2_over_T4']=float((4*f*self.vol)@self.mean)
            s['bath_energy_fraction_g106p75']=float((s['rho_phi_over_T4']+s['rho_chi2_over_T4'])/(np.pi**2/30*106.75))
            s['delta_Neff_chi2_if_stable_free_streaming']=float(s['rho_chi2_over_T4']/(7*np.pi**2/120)*(3.909/106.75)**(4/3)/(4/11)**(4/3))
            s['omega_phi_h2_m86_g106p75']=float(omega_from_number(s['number']))
            assert abs(s['number_mismatch'])<1e-7 and s['min_occupation']>-1e-7, s
            assert s['max_fermion']<1.00001,s
            history.append(s)
        result=dict(n=self.n,kappa=self.kappa,x0=x0,xend=xend,qmin=self.edges[0],qmax=self.edges[-1],
                    rtol=rtol,nfev=sol.nfev,history=history,final=history[-1])
        if keep_psd:
            result['frozen_cell_spectrum']={'edges':self.edges.tolist(),'f_phi':sol.y[:self.n,-1].tolist(),
                'f_chi2':sol.y[self.n:,-1].tolist(),'g_phi':1,'g_chi2':4}
        return result


def omega_from_number(number,m_eV=86.,gs=106.75):
    H100=1e5/(1e6*parsec);rho=3*H100**2/(8*np.pi*G)
    n_today=number*(Boltzmann*2.7255/(hbar*c))**3*3.909/gs
    return m_eV*elementary_charge/c**2*n_today/rho


def main():
    p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=200)
    p.add_argument('--kappa',type=float,nargs='+',default=[.001,.01,.1]);p.add_argument('--x0',type=float,default=.1)
    p.add_argument('--xend',type=float,default=30.);p.add_argument('--qmin',type=float,default=1e-7)
    p.add_argument('--qmax',type=float,default=40.);p.add_argument('--rtol',type=float,default=2e-6)
    p.add_argument('--output',required=True);p.add_argument('--keep-psd',action='store_true');a=p.parse_args();results=[]
    for k in a.kappa:
        model=RadiationFixedBath(a.n,a.qmin,a.qmax,k)
        # Independent check of analytic area vs quadrature in the flat FV solver.
        for x in [.1,1,5]:
            model.update(x)
            reference=DaughterEnergyFV(a.n,x,a.qmin,a.qmax,bath='balanced')
            denom=np.max(reference.R)
            assert np.max(np.abs(model.R/(k*x)-reference.R))/denom<1e-11
            eq=np.r_[1/np.expm1(model.mean),expit(-model.mean)]
            assert np.sum(np.abs(model.rhs(x,eq)[:a.n])*model.vol)<1e-11*max(model.A.sum(),1e-30)
        row=model.evolve(a.x0,a.xend,a.rtol,a.keep_psd);results.append(row)
        print(json.dumps({'kappa':k,'n':a.n,'x0':a.x0,'nfev':row['nfev'],**row['final']}),flush=True)
    Path(a.output).write_text(json.dumps(results,indent=2)+'\n')

if __name__=='__main__':main()
