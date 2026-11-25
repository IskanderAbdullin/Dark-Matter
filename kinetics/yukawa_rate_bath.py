"""Yukawa-motivated rate scaling, NOT full finite-temperature field theory.

Gamma_i^SM = y_i^2 M_i/(8 pi), Gamma_th = gamma y_i^2 T.
Use explicit vacuum inverse-decay relaxation and a smoothly suppressed
relativistic scattering estimate. The suppression is phenomenological.
All dimensionless rates are referred to H_ref(M), not the actual H.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.special import expit
from daughter_bath_relaxation import DaughterBath


class YukawaRateBath(DaughterBath):
    def __init__(self,kappa_sm1=100.,gamma=.01,scatter=True,**kwargs):
        super().__init__(**kwargs)
        assert kappa_sm1>0 and gamma>=0
        self.kappa_sm1=kappa_sm1;self.gamma=gamma;self.scatter=scatter

    def bath_rates(self,x,th):
        theta=th['theta'];Hf=th['Hfactor'];r=self.fermion_ratio;b=x*r
        eq=expit(-th['E']/theta)
        ne=th['Neq'];net=th['NeqT'];w=self.parent_w
        # Equilibrium average m1/E1 (number relaxation approximation).
        I=x*float(w@(eq/th['E']))
        IT=x*float(w@(eq*(1-eq)))/theta
        if ne>1e-200:
            dil1=I/ne;dil1T=(IT-dil1*net)/ne
        else:
            dil1=1.;dil1T=0.
        lam1=self.kappa_sm1*x*x*Hf*dil1
        lam1T=self.kappa_sm1*x*x*Hf*dil1T
        dil2=self.trace_weights[self.n:]/(4*self.vol*b)
        rate2=self.kappa2*x*x*Hf*dil2;rate2T=np.zeros(self.n)
        if self.scatter:
            # Relativistic gamma*y^2*T estimate is only valid at T>>M_i.
            # The rational switch is a deliberately explicit interpolation,
            # NOT a calculation of nonrelativistic scattering.
            U=x/theta;V=b/theta
            sc1=8*np.pi*self.gamma*self.kappa_sm1*x*theta*Hf/(1+U*U)
            sc2=8*np.pi*self.gamma*self.kappa2*x*theta*Hf/(r*(1+V*V))
            lam1+=sc1;lam1T+=sc1*(1+2*U*U/(1+U*U))
            rate2=rate2+sc2;rate2T+=sc2*(1+2*V*V/(1+V*V))
        return rate2,rate2T,lam1,lam1T,True

    def run(self,x0,xend,rtol,fermion_atol,keep_psd,parent_gev,initial_parent='thermal'):
        r=self.evolve_bath(x0,xend,rtol,fermion_atol,keep_psd,parent_gev,initial_parent)
        # The inherited constant relaxation setting is not used.
        r.pop('relaxation_per_H')
        r.update(kappa_sm1=self.kappa_sm1,gamma=self.gamma,scatter=self.scatter,
                 rate_model='vacuum time-dilated decay + phenomenological tapered relativistic scattering')
        self.geometry(x0)
        N=r['initial_mother_number']
        y=np.r_[np.zeros(2*self.n),N,0.,0.,0.]
        _,th,rate,eq,bath=self.evaluate(np.log(x0),y)
        r['initial_parent_relaxation_over_H']=float(th['bath_rate1'])
        r['initial_fermion_relaxation_over_H_q3']=float(np.interp(3.,self.mean,rate))
        return r


def rate_mapping(M,r,kappa,kappa2,kappa_sm1,gstar=106.75):
    """Illustrative tree-level mapping in specified width conventions.

    Spin-averaged Dirac pseudoscalar interaction i*g*phi*bar(chi2)*gamma5*chi1.
    For mphi negligible, Gamma_phi=g^2 M (1-r^2)(1-r)^2/(16 pi).
    If g=C*(1+r)*M/(2 f), then Gamma_phi=C^2 M^3 (1-r^2)^3/(64 pi f^2).
    SM portal width convention is y_i^2 M_i/(8 pi).
    """
    mpl=1.220890e19;H=1.66*np.sqrt(gstar)*M*M/mpl
    eta=1-r*r
    fc=np.sqrt(M**3*eta**3/(64*np.pi*kappa*H))
    g=(1+r)*M/(2*fc)
    y1=np.sqrt(8*np.pi*kappa_sm1*H/M)
    y2=np.sqrt(8*np.pi*kappa2*H/(r*M))
    return dict(parent_GeV=M,H_ref_GeV=float(H),f_over_C_GeV=float(fc),
                pseudoscalar_coupling=float(g),y1=float(y1),y2=float(y2),
                nominal_seesaw_mass1_eV=float(y1*y1*174.**2/M*1e9),
                nominal_seesaw_mass2_eV=float(y2*y2*174.**2/(r*M)*1e9))


def main():
    p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=160)
    p.add_argument('--fermion-ratio',type=float,default=.9)
    p.add_argument('--kappa',type=float,required=True);p.add_argument('--kappa2',type=float,default=.03)
    p.add_argument('--kappa-sm1',type=float,nargs='+',default=[100.,1000.,4000.])
    p.add_argument('--gamma',type=float,default=.01);p.add_argument('--no-scatter',action='store_true')
    p.add_argument('--parent-gev',type=float,default=86000.)
    p.add_argument('--x0',type=float,default=.1);p.add_argument('--xend',type=float,default=80.)
    p.add_argument('--rtol',type=float,default=2e-6);p.add_argument('--qmin',type=float,default=1e-7)
    p.add_argument('--initial-parent',choices=['thermal','empty'],default='thermal')
    p.add_argument('--keep-psd',action='store_true');p.add_argument('--output',required=True);a=p.parse_args();out=[]
    for k1 in a.kappa_sm1:
        m=YukawaRateBath(n=a.n,scalar_ratio=86./(a.parent_gev*1e9),fermion_ratio=a.fermion_ratio,
                        kappa=a.kappa,kappa2=a.kappa2,kappa_sm1=k1,gamma=a.gamma,scatter=not a.no_scatter,qmin=a.qmin)
        row=m.run(a.x0,a.xend,a.rtol,1e-15,a.keep_psd,a.parent_gev,a.initial_parent)
        row['mapping']=rate_mapping(a.parent_gev,a.fermion_ratio,a.kappa,a.kappa2,k1)
        out.append(row);s=row['final']
        print(json.dumps(dict(n=a.n,kappa=a.kappa,kappa2=a.kappa2,kappa_sm1=k1,scatter=not a.no_scatter,
                             b=a.fermion_ratio,initial_lambda1=row['initial_parent_relaxation_over_H'],
                             initial_lambda2_q3=row['initial_fermion_relaxation_over_H_q3'],
                             mapping=row['mapping'],checks=row['checks'],**s)),flush=True)
        Path(a.output).write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':main()
