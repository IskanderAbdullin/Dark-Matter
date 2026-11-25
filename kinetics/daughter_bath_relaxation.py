"""Detailed-balanced effective chi2 bath interaction with finite vacuum masses.

C2=(Gamma2/H)*(m2/E2)*(f_FD(T)-f2), a phenomenological relaxation
operator representing decays and inverse creation by an equilibrated bath.
It is not a thermal field-theory calculation of a specified Lagrangian.
The exchanged energy heats/cools the same bath as the parent chemistry.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.special import expit
from scipy.integrate import solve_ivp
from two_mass_bath import TwoMassBath
from radiation_fixed_bath import omega_from_number


class DaughterBath(TwoMassBath):
    def __init__(self,kappa2=.1,**kwargs):
        super().__init__(**kwargs)
        assert self.fermion_ratio>0 and kappa2>=0
        self.kappa2=kappa2

    def bath_rates(self,x,th):
        # Rate2 scales with Hfactor; rate1 in this baseline is constant per H.
        b=x*self.fermion_ratio
        dilution=self.trace_weights[self.n:]/(self.g_fermion*self.vol*b)
        rate=self.kappa2*x*x*th['Hfactor']*dilution
        return rate,np.zeros(self.n),self.relaxation,0.,False

    def evaluate(self,z,y):
        x=np.exp(z);n=self.n;N=y[-4]
        th=self.set_state(x,y[:-1])
        yd=y[:2*n];f=yd[n:]
        rate,rateT,lam,lamT,lamH=self.bath_rates(x,th)
        th.update(bath_rate2_T=rateT,bath_rate1=lam,bath_rate1_T=lamT,bath_rate1_Hdependent=lamH)
        eq=expit(-self.energy_chi/th['theta'])
        bath=rate*(eq-f)
        dy=self.collision_rhs(yd);dy[n:]+=bath
        th['parent_bath_source']=lam*(th['Neq']-N)
        dN=-self.vol@dy[:n]+th['parent_bath_source']
        dtheta=(th['force']-th['rhoN']*dN-self.ew@dy)/th['heat']
        dS=self.g_fermion*self.vol@bath
        dW=th['trace']+self.trace_weights@yd
        out=np.r_[dy,dN,dtheta,dW,dS]
        return out,th,rate,eq,bath

    def rhs_full(self,z,y):
        return self.evaluate(z,y)[0]

    def jac_full(self,z,y):
        x=np.exp(z);n=self.n;N=y[-4];yd=y[:2*n];p,f=np.split(yd,2);v=self.vol
        dyfull,th,rate,eq,bath=self.evaluate(z,y)
        dy=dyfull[:2*n];dS=dyfull[-1]
        J=self.collision_jac(yd)
        J[n:,n:]-=np.diag(rate)
        J-=.5/th['total']*np.outer(dy,self.ew)
        F=th['F'];S=F*(1-F)
        if N>0 and th['A']>1e-250:
            FN=S/th['A']
        elif N==0:
            FN=np.exp(np.minimum(700,(x-self.pair_E)/th['theta']))/float(self.parent_w@np.exp((x-th['E'])/th['theta']))
        else:
            FN=np.zeros_like(F)
        FT=S*(self.pair_E-th['rhoN'])/th['theta']
        def col(DF,total_deriv):
            dA=self.R*DF
            return np.r_[(dA@(1-f)+p*dA.sum(axis=1))/v,
                         (dA.T@(1+p)-f*dA.sum(axis=0))/(self.g_fermion*v)]-.5*dy*total_deriv/th['total']
        jN=col(FN,th['rhoN']);jT=col(FT,th['heat'])
        temp=rate*eq*(1-eq)*self.energy_chi/th['theta']+th['bath_rate2_T']*(eq-f)
        jT[n:]+=temp
        JD=np.column_stack((J,jN,jT,np.zeros((2*n,2))))
        jParent=-v@JD[:n,:]
        jParent[-4]-=th['bath_rate1']
        jParent[-3]+=th['bath_rate1']*th['NeqT']+th['bath_rate1_T']*(th['Neq']-N)
        if th['bath_rate1_Hdependent']:
            jParent-=.5*th['parent_bath_source']/th['total']*np.r_[self.ew,th['rhoN'],th['heat'],0.,0.]
        jTemp=(-th['rhoN']*jParent-self.ew@JD)/th['heat']
        jTrace=np.r_[self.trace_weights,0.,0.,0.,0.]
        jS=np.r_[np.zeros(n),-self.g_fermion*v*rate,0.,self.g_fermion*v@temp,0.,0.]
        jS-=.5*dS/th['total']*np.r_[self.ew,th['rhoN'],th['heat'],0.,0.]
        for k in [-4,-3]:
            step=max(1e-10,abs(y[k])*1e-5) if k==-4 else 1e-6
            yp=y.copy();ym=y.copy();yp[k]+=step;ym[k]-=step
            if k==-4 and ym[k]<0:
                diff=(self.rhs_full(z,yp)-self.rhs_full(z,y))/step
            else:
                diff=(self.rhs_full(z,yp)-self.rhs_full(z,ym))/(2*step)
            jTemp[k]=diff[-3];jTrace[k]=diff[-2]
        return np.vstack((JD,jParent,jTemp,jTrace,jS))

    def test_full(self):
        x=.2;z=np.log(x);self.geometry(x);N=self.thermo(x,0.,0.)['Neq']
        y=np.r_[1/np.expm1(self.energy_phi),expit(-self.energy_chi),N,0.,0.,0.]
        dy,th,rate,eq,bath=self.evaluate(z,y)
        bal=(np.abs(dy[:self.n])@self.vol+np.abs(dy[self.n:2*self.n])@(self.g_fermion*self.vol))/max(self.A.sum(),1e-30)
        assert bal<1e-8,bal
        y=np.r_[np.full(self.n,.1),np.full(self.n,.05),N*.8,.01,0.,0.]
        rng=np.random.default_rng(17);d=rng.normal(size=len(y));d/=np.linalg.norm(d)
        eps=1e-6;j=self.jac_full(z,y)@d
        fd=(self.rhs_full(z,y+eps*d)-self.rhs_full(z,y-eps*d))/(2*eps)
        err=np.linalg.norm(j-fd)/max(np.linalg.norm(j),1e-30)
        assert err<1e-5,err
        return dict(balance=float(bal),jacobian=float(err))

    def evolve_bath(self,x0=.1,xend=80.,rtol=2e-6,fermion_atol=1e-15,keep_psd=False,parent_gev=86.,initial_parent='thermal'):
        checks=self.test_full();n=self.n
        assert initial_parent in ('thermal','empty')
        init=self.thermo(x0,0.,0.)['Neq'] if initial_parent=='thermal' else 0.
        rho0=self.thermo(x0,init,0.)['rho']+self.rad0
        y0=np.r_[np.zeros(2*n),init,0.,0.,0.]
        samples=sorted(set([x0]+[x for x in [.01,.02,.03,.05,.075,.1,.15,.2,.3,.5,1,2,3,5,10,15,20,30,40,60,80,100,120] if x0<x<xend]+[xend]))
        atol=np.r_[np.full(n,1e-12),np.full(n,fermion_atol),1e-14,1e-13,1e-12,1e-13]
        sol=solve_ivp(self.rhs_full,[np.log(x0),np.log(xend)],y0,method='BDF',
                      jac=self.jac_full,rtol=rtol,atol=atol,max_step=.035,t_eval=np.log(samples))
        if not sol.success:raise RuntimeError(sol.message+'; last sampled x='+str(np.exp(sol.t[-1])))
        history=[];l,h=self.edges[:-1],self.edges[1:]
        for x,y in zip(samples[1:],sol.y.T[1:]):
            self.geometry(x)
            yd=y[:2*n];p,f=np.split(yd,2);N,zt,W,S=y[-4:]
            th=self.thermo(x,N,zt);theta=th['theta'];s=self.snapshot(yd)
            N2=float(self.g_fermion*self.vol@f);Nphi=s['number']
            rho1=th['rho'];rphi=float((p*self.vol)@self.energy_phi);rchi=float((self.g_fermion*f*self.vol)@self.energy_chi)
            total=th['rb']+rho1+rphi+rchi
            source_err=(N2-Nphi-S)/max(Nphi,1e-12)
            cold=np.maximum(np.minimum(h,.1*theta)**3-l**3,0)/(h**3-l**3)
            s.update(x=x,theta=theta,mother_number=float(N),mother_rho_Tref4=rho1,
                     fermion_number=N2,integrated_fermion_bath_source=float(S),
                     number_source_residual=float(source_err),
                     rho_phi_Tref4=rphi,rho_chi2_Tref4=rchi,total_rho_Tref4=float(total),
                     continuity_residual=float((total-rho0-W)/rho0),
                     matter_fraction_chi2=float(rchi/total),
                     mean_p_over_actual_T=s['mean_q']/theta,
                     hot_fraction_p_over_actual_T_gt_0p1=float((p*self.vol)@(1-cold)/Nphi),
                     omega_phi_h2_m86=float(omega_from_number(Nphi/theta**3,gs=self.g_bath)),
                     omega_chi2_h2_if_no_more_reactions_NR=float(omega_from_number(N2/theta**3,m_eV=self.fermion_ratio*parent_gev*1e9,gs=self.g_bath)))
            # Number difference alone is no longer an invariant.
            s.pop('number_mismatch')
            assert abs(source_err)<1e-6,s
            assert abs(s['continuity_residual'])<5e-5,s
            assert s['max_fermion']<1.00001 and s['min_occupation']>-1e-7 and N>-1e-9,s
            history.append(s)
        out=dict(n=n,x0=x0,xend=xend,kappa=self.kappa,kappa2=self.kappa2,scalar_ratio=self.scalar_ratio,
                 fermion_ratio=self.fermion_ratio,parent_GeV=parent_gev,relaxation_per_H=self.relaxation,
                 g_bath=self.g_bath,g_fermion=self.g_fermion,rtol=rtol,fermion_atol=fermion_atol,checks=checks,nfev=sol.nfev,
                 initial_parent=initial_parent,initial_mother_number=init,initial_total_rho_Tref4=rho0,history=history,final=history[-1])
        if keep_psd:
            out['frozen_cell_spectrum']=dict(edges=self.edges.tolist(),f_phi=sol.y[:n,-1].tolist(),
                f_chi2=sol.y[n:2*n,-1].tolist(),g_phi=1,g_chi2=self.g_fermion,theta_freezeout=history[-1]['theta'])
        return out


def main():
    p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=160)
    p.add_argument('--scalar-ratio',type=float,default=1e-9);p.add_argument('--fermion-ratio',type=float,default=.7071067811865476)
    p.add_argument('--kappa',type=float,default=.05);p.add_argument('--kappa2',type=float,nargs='+',default=[.01,.03,.1])
    p.add_argument('--x0',type=float,default=.1);p.add_argument('--xend',type=float,default=80.)
    p.add_argument('--qmin',type=float,default=1e-7);p.add_argument('--parent-gev',type=float,default=86.)
    p.add_argument('--relaxation',type=float,default=100.);p.add_argument('--rtol',type=float,default=2e-6)
    p.add_argument('--fermion-atol',type=float,default=1e-15);p.add_argument('--keep-psd',action='store_true')
    p.add_argument('--output',required=True);a=p.parse_args();out=[]
    for k2 in a.kappa2:
        m=DaughterBath(n=a.n,scalar_ratio=a.scalar_ratio,fermion_ratio=a.fermion_ratio,kappa=a.kappa,
                       kappa2=k2,relaxation=a.relaxation,qmin=a.qmin)
        r=m.evolve_bath(a.x0,a.xend,a.rtol,a.fermion_atol,a.keep_psd,a.parent_gev);out.append(r)
        print(json.dumps(dict(n=a.n,kappa=a.kappa,kappa2=k2,fermion_ratio=a.fermion_ratio,
                             checks=r['checks'],**r['final'])),flush=True)
        Path(a.output).write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':main()
