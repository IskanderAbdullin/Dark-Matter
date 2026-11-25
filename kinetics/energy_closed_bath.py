"""FRW kinetic-equilibrium parent + finite chemical source + closed stress-energy.
Phenomenological relaxation; constant bath degrees of freedom; massless daughters at production.
"""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.special import expit,roots_legendre
from daughter_energy_fv import DaughterEnergyFV
from radiation_fixed_bath import RadiationFixedBath,omega_from_number


class EnergyClosedBath(RadiationFixedBath):
    def __init__(self,n=200,qmin=1e-7,qmax=40.,kappa=.1347,relaxation=10.,g_bath=106.75):
        super().__init__(n,qmin,qmax,kappa)
        self.relaxation=relaxation;self.g_bath=g_bath;self.rad0=np.pi**2/30*g_bath
        p,w=roots_legendre(160);p=40*(p+1);w=40*w
        self.parent_p=p;self.parent_w=4/(2*np.pi**2)*p*p*w
        self.pair_E=self.mean[:,None]+self.mean[None,:]
        self.ew=np.r_[self.vol*self.mean,4*self.vol*self.mean]
        self.geom_x=None
        self.energy_phi=self.mean
        self.trace_weights=np.zeros(2*n)

    def geometry(self,x):
        if x!=self.geom_x:
            super().update(x);self.R0=self.R*x;self.geom_x=x

    def collision_rhs(self,y):
        p,f=np.split(y,2)
        return np.r_[(self.A@(1-f)+p*(self.A.sum(axis=1)-self.R@f))/self.vol,
                     (self.A.T@(1+p)-f*(self.A.sum(axis=0)+self.R.T@p))/(4*self.vol)]

    def collision_jac(self,y):
        p,f=np.split(y,2);v=self.vol
        return np.block([[np.diag((self.A.sum(axis=1)-self.R@f)/v),
                          -(self.A+p[:,None]*self.R)/v[:,None]],
                         [(self.A.T-f[:,None]*self.R.T)/(4*v[:,None]),
                          -np.diag((self.A.sum(axis=0)+self.R.T@p)/(4*v))]])

    def thermo(self,x,N,ztheta):
        theta=np.exp(ztheta);E=np.hypot(self.parent_p,x);w=self.parent_w
        eq=expit(-E/theta);Neq=float(w@eq);NeqT=float(w@(eq*(1-eq)*E)/theta)
        if N>0:
            xi=brentq(lambda u:w@expit(u-E/theta)-N,x/theta-750,x/theta+80,xtol=1e-12)
            f=expit(xi-E/theta)
        else:
            xi=-np.inf;f=np.zeros_like(E)
        s=f*(1-f);A=float(w@s);B=float(w@(s*E));D=float(w@(s/E))
        rho=float(w@(f*E));trace=float(x*x*(w@(f/E)))
        if A>1e-250:
            rhoN=B/A
            rhoT=float(w@(s*(E-rhoN)**2))/theta
            force=x*x*(A-rhoN*D)/theta
        else:
            boltz=w*np.exp((x-E)/theta);rhoN=float(boltz@E/boltz.sum());rhoT=force=0.
        rb=self.rad0*theta**4
        return dict(theta=theta,xi=xi,A=A,rho=rho,trace=trace,rhoN=rhoN,rhoT=rhoT,force=force,
                    rb=rb,Neq=Neq,NeqT=NeqT,E=E,heat=4*rb+rhoT)

    def set_state(self,x,y):
        n=self.n;yd=y[:2*n];N,zt=y[-3:-1]
        self.geometry(x)
        th=self.thermo(x,N,zt);total=th['rb']+th['rho']+self.ew@yd
        Hfactor=np.sqrt(self.rad0/total)
        self.R=self.R0*Hfactor
        F=expit(th['xi']-self.pair_E/th['theta'])
        self.A=self.R*F;self.asum=self.A.sum(axis=1)
        th.update(total=total,Hfactor=Hfactor,F=F)
        return th

    def rhs_log(self,z,y):
        x=np.exp(z);n=self.n;N=y[-3];th=self.set_state(x,y)
        dy=self.collision_rhs(y[:2*n])
        dN=-self.vol@dy[:n]+self.relaxation*(th['Neq']-N)
        dtheta=(th['force']-th['rhoN']*dN-self.ew@dy)/th['heat']
        return np.r_[dy,dN,dtheta,th['trace']+self.trace_weights@y[:2*n]]

    def jac_log(self,z,y):
        x=np.exp(z);n=self.n;N=y[-3];th=self.set_state(x,y);p,f=np.split(y[:2*n],2)
        dy=self.collision_rhs(y[:2*n])
        J=self.collision_jac(y[:2*n])-.5/th['total']*np.outer(dy,self.ew)
        F=th['F'];S=F*(1-F)
        if N>0 and th['A']>1e-250:
            FN=S/th['A']
        elif N==0:
            FN=np.exp(np.minimum(700,(x-self.pair_E)/th['theta']))/float(self.parent_w@np.exp((x-th['E'])/th['theta']))
        else:
            FN=np.zeros_like(F)
        FT=S*(self.pair_E-th['rhoN'])/th['theta']
        def col(DF,total_deriv):
            dA=self.R*DF;dS=dA.sum(axis=1)
            return np.r_[(dA@(1-f)+p*dS)/self.vol,(dA.T@(1+p)-f*dA.sum(axis=0))/(4*self.vol)]-.5*dy*total_deriv/th['total']
        jN=col(FN,th['rhoN']);jT=col(FT,th['heat'])
        JD=np.column_stack((J,jN,jT,np.zeros(2*n)))
        jParent=-self.vol@JD[:n,:]
        jParent[-3]-=self.relaxation;jParent[-2]+=self.relaxation*th['NeqT']
        jTemp=(-th['rhoN']*jParent-self.ew@JD)/th['heat']
        jTrace=np.r_[self.trace_weights,0.,0.,0.]
        # Only two thermodynamic columns require finite differences of higher moments.
        for k in [-3,-2]:
            step=max(1e-10,abs(y[k])*1e-5) if k==-3 else 1e-6
            yp=y.copy();ym=y.copy();yp[k]+=step;ym[k]-=step
            if k==-3 and ym[k]<0:
                deriv=(self.rhs_log(z,yp)-self.rhs_log(z,y))/step
            else:
                deriv=(self.rhs_log(z,yp)-self.rhs_log(z,ym))/(2*step)
            jTemp[k]=deriv[-2];jTrace[k]=deriv[-1]
        return np.vstack((JD,jParent,jTemp,jTrace))

    def evolve_closed(self,x0=.1,xend=40.,rtol=2e-6,keep_psd=False):
        th0=self.thermo(x0,0.,0.);init=th0['Neq'];initial_rho=self.thermo(x0,init,0.)['rho']+self.rad0
        y0=np.r_[np.zeros(2*self.n),init,0.,0.]
        samples=sorted(set([x0]+[x for x in [.03,.05,.075,.1,.15,.2,.3,.5,1,2,3,5,10,15,20,30,40,60] if x0<x<xend]+[xend]))
        sol=solve_ivp(self.rhs_log,[np.log(x0),np.log(xend)],y0,method='BDF',jac=self.jac_log,
                      rtol=rtol,atol=1e-12,max_step=.04,t_eval=np.log(samples))
        if not sol.success:raise RuntimeError(sol.message)
        history=[];l,h=self.edges[:-1],self.edges[1:]
        for x,y in zip(samples[1:],sol.y.T[1:]):
            self.geometry(x)
            yd=y[:2*self.n];N,zt,W=y[-3:];th=self.thermo(x,N,zt);theta=th['theta'];s=self.snapshot(yd)
            p,f=np.split(yd,2);total=th['rb']+th['rho']+self.ew@yd
            cold=np.maximum(np.minimum(h,.1*theta)**3-l**3,0)/(h**3-l**3)
            rphi=float((p*self.vol)@self.energy_phi);rchi=float((4*f*self.vol)@self.mean)
            s.update(x=x,theta=theta,mother_number=float(N),mother_eq_number=th['Neq'],mother_rho_Tref4=th['rho'],
                     total_rho_Tref4=float(total),continuity_residual=float((total-initial_rho-W)/initial_rho),
                     mother_plus_phi_number=float(N+s['number']),rho_phi_Tref4=rphi,rho_chi2_Tref4=rchi,
                     mean_p_over_actual_T=s['mean_q']/theta,hot_fraction_p_over_actual_T_gt_0p1=float((p*self.vol)@(1-cold)/s['number']),
                     omega_phi_h2_m86=float(omega_from_number(s['number']/theta**3,gs=self.g_bath)),
                     delta_Neff_chi2_if_stable=float(rchi/theta**4/(7*np.pi**2/120)*(3.909/self.g_bath)**(4/3)/(4/11)**(4/3)))
            assert abs(s['number_mismatch'])<1e-7 and N>-1e-9 and s['min_occupation']>-1e-7 and s['max_fermion']<1.00001,s
            assert abs(s['continuity_residual'])<3e-5,s
            if self.relaxation==0:assert abs(s['mother_plus_phi_number']/init-1)<1e-7,s
            history.append(s)
        out=dict(n=self.n,x0=x0,xend=xend,kappa=self.kappa,relaxation_per_H=self.relaxation,g_bath=self.g_bath,
                 initial_mother_number=init,initial_total_rho_Tref4=initial_rho,nfev=sol.nfev,rtol=rtol,history=history,final=history[-1])
        if keep_psd:out['frozen_cell_spectrum']={'edges':self.edges.tolist(),'f_phi':sol.y[:self.n,-1].tolist(),
            'f_chi2':sol.y[self.n:2*self.n,-1].tolist(),'g_phi':1,'g_chi2':4,'theta_freezeout':history[-1]['theta']}
        return out


def main():
    p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=180);p.add_argument('--x0',type=float,default=.1)
    p.add_argument('--xend',type=float,default=40.);p.add_argument('--kappa',type=float,default=.1347)
    p.add_argument('--relaxation',type=float,nargs='+',default=[0.,10.,100.]);p.add_argument('--rtol',type=float,default=2e-6)
    p.add_argument('--output',required=True);p.add_argument('--keep-psd',action='store_true');a=p.parse_args();rows=[]
    for lam in a.relaxation:
        m=EnergyClosedBath(a.n,kappa=a.kappa,relaxation=lam)
        x=.2;z=np.log(x);y=np.r_[np.full(m.n,.1),np.full(m.n,.05),m.thermo(x,0.,0.)['Neq']*.8,.01,0.]
        rng=np.random.default_rng(17);d=rng.normal(size=len(y));d/=np.linalg.norm(d)
        eps=1e-6;j=m.jac_log(z,y)@d;fd=(m.rhs_log(z,y+eps*d)-m.rhs_log(z,y-eps*d))/(2*eps)
        assert np.linalg.norm(j-fd)/max(np.linalg.norm(j),1e-30)<1e-5
        r=m.evolve_closed(a.x0,a.xend,a.rtol,a.keep_psd);rows.append(r)
        print(json.dumps({'relaxation':lam,'n':a.n,'nfev':r['nfev'],**r['final']}),flush=True)
    Path(a.output).write_text(json.dumps(rows,indent=2,allow_nan=False)+'\n')

if __name__=='__main__':main()
