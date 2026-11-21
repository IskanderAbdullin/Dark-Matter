"""Radiation-dominated FRW fixed-mother-bath 1<->2 pilot in comoving q=p/T.

Assumptions: radiation H=H1/x^2, x=M1/T proportional to scale factor,
constant relativistic entropy degrees, mother held at thermal FD equilibrium,
massless daughters, no self-consistent temperature or mother depletion.
Gamma_rest/H(T=M1) = kappa is an externally chosen dimensionless coupling.

Exact continuum daughter-energy measure:
  (g1 M1 Gamma/(2 pi^2)) dq_phi_phys dq_chi_phys
with q_phys=T*q_comoving and allowed 4 q_phi q_chi >= x^2.

Because dx/dt=H*x, the occupation ODE is:
  df/dx = (g1/(2 pi^2)) [Gamma/H(x)] * paired_AREA(x)/VOL
           * [a(1-f_chi)(1+f_phi)-(1-a) f_chi f_phi]
        = (g1 kappa x^2/(2 pi^2)) [AREA(x)/VOL] * statistics.
Liouville redshifting vanishes for q=p/T when T proportional to 1/a.

THIS IS NOT self-consistent cosmology or a physical CLASS initial PSD.
"""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import expit
from onshell_collision_2d import allowed_rectangle_areas


class FRWFixedBathPilot:
    def __init__(self,nq=80,kappa=0.01,qmin=1e-5,qmax=32.,g1=4,gchi=4):
        if nq<10 or kappa<0 or qmin<=0 or qmax<=qmin:
            raise ValueError("Bad FRW toy parameters")
        self.nq,self.kappa,self.g1,self.gchi=nq,kappa,g1,gchi
        self.edges=np.geomspace(qmin,qmax,nq+1)
        self.e=.75*(self.edges[1:]**4-self.edges[:-1]**4)/(self.edges[1:]**3-self.edges[:-1]**3)
        self.vphi=(self.edges[1:]**3-self.edges[:-1]**3)/(6*math.pi**2)
        self.vchi=gchi*self.vphi
        self.E=self.e[:,None]+self.e[None,:]
        self.a=expit(-self.E)
        self.pair_product=4*self.e[:,None]*self.e[None,:]

    def flux(self,x,b,f):
        if x<=0:raise ValueError("x=M1/T must be positive")
        area=allowed_rectangle_areas(self.edges,M1_T=x)
        on_shell=self.pair_product>=x*x
        K=(self.g1*self.kappa*x*x/(2*math.pi**2))*area*on_shell
        A=K*self.a*(1+b[:,None])*(1-f[None,:])
        B=K*(1-self.a)*b[:,None]*f[None,:]
        return A-B

    def rhs(self,x,y):
        b,f=y[:self.nq],y[self.nq:]
        J=self.flux(x,b,f)
        return np.r_[J.sum(axis=1)/self.vphi,J.sum(axis=0)/self.vchi]

    def snapshot(self,x,y):
        b,f=y[:self.nq],y[self.nq:]
        mass=b*self.vphi
        nphi=float(mass.sum())
        nchi=float(np.dot(f,self.vchi))
        d=self.rhs(x,y)
        return {
          "x_M1_over_T":float(x),
          "Gamma_over_H":float(self.kappa*x*x),
          "nphi_over_T3":nphi,
          "nchi2_over_T3":nchi,
          "number_difference_over_T3":nphi-nchi,
          "mean_phi_q_over_T":float(np.dot(mass,self.e)/nphi) if nphi>0 else None,
          "phi_fraction_q_below_0p1":float(mass[self.e<.1].sum()/nphi) if nphi>0 else None,
          "phi_occupation_max_grid_dependent":float(np.max(b)),
          "fermion_occupation_max":float(np.max(f)),
          "min_occupation":float(np.min(y)),
          "nphi_source_d_by_dx":float(np.dot(self.vphi,d[:self.nq])),
          "number_rate_difference_d_by_dx":float(np.dot(self.vphi,d[:self.nq])-
                                                   np.dot(self.vchi,d[self.nq:]))
        }

    def integrate(self,x_start=.05,x_mid=20.,x_end=30.,
                  rtol=1e-5,atol=1e-10,max_step=.05):
        if not (0<x_start<x_mid<x_end):
            raise ValueError("Require 0 < x_start < x_mid < x_end")
        init=np.zeros(2*self.nq)
        sol=solve_ivp(self.rhs,(x_start,x_end),init,method="BDF",
                      t_eval=[x_mid,x_end],rtol=rtol,atol=atol,max_step=max_step)
        if not sol.success:raise RuntimeError(sol.message)
        mids=self.snapshot(x_mid,sol.y[:,0])
        final=self.snapshot(x_end,sol.y[:,1])
        if abs(final["number_difference_over_T3"])>2e-7:
            raise AssertionError("Daughter number invariant failed")
        if final["min_occupation"]<-2e-8 or final["fermion_occupation_max"]>1.000001:
            raise AssertionError("Statistical occupation bounds failed")
        drift=abs(final["nphi_over_T3"]/mids["nphi_over_T3"]-1) if mids["nphi_over_T3"]>0 else None
        return {
          "integration":{"x_start":x_start,"x_mid":x_mid,"x_end":x_end,
                         "rtol":rtol,"atol":atol,"max_step":max_step,
                         "nfev":sol.nfev,"success":sol.success},
          "mid":mids,"final":final,
          "relative_number_change_from_mid_to_end":drift,
          "frozen_number_to_1e-4_over_late_window":drift is not None and drift<1e-4
        }


def run(quick=False):
    cases=([(40,.001),(40,.01)] if quick else
           [(40,.001),(80,.0001),(80,.001),(80,.01),
            (80,.05),(80,.1),(80,.2),(120,.01),(120,.1),(120,.2)])
    data=[]
    for nq,kappa in cases:
        model=FRWFixedBathPilot(nq=nq,kappa=kappa)
        result=model.integrate()
        data.append({"nq":nq,"Gamma_over_H_at_x1":kappa,"result":result})
    return {
      "name":"Radiation-dominated fixed-bath FRW kinetic PILOT; not self-consistent cosmology",
      "x_definition":"x=M1/T, q=p/T comoving, H(x)=H(T=M1)/x^2",
      "kappa_definition":"kappa=Gamma_rest/H(T=M1)",
      "mother_status":"externally maintained thermal Fermi-Dirac bath",
      "gstar_s_constant_assumption":True,
      "daughter_masses_zero_during_production":True,
      "source_energy_budget_not_closed":True,
      "not_valid_CLASS_initial_spectrum":True,
      "warning":"The very cold instantaneous fixed-T burst is not necessarily a frozen FRW spectrum. Investigate mother depletion, thermal history, IR dynamics and energy budget before interpreting omega_DM.",
      "cases":data
    }


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--quick",action="store_true")
    p.add_argument("--output",default="outputs/frw_fixed_bath_pilot.json")
    opts=p.parse_args()
    ans=run(opts.quick)
    path=Path(opts.output);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(ans,indent=2),encoding="utf-8")
    print(json.dumps(ans,indent=2))
