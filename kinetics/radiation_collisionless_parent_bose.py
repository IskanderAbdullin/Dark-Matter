"""Independent exact initial Bose-gain kernel on a prescribed radiation FRW.

This is an analytically reduced OPEN approximation to BEST's full CLOSED
chi1 <-> chi2 + massive phi kinetic system. It retains collisionless
redshift of the already-produced massive mother's frozen initial FD
COMOVING distribution, rather than incorrectly imposing FD at each T.

RD a(t)=sqrt(t/t0), initial mother f1(q1) =
  1/[1+exp(sqrt(q1²+M1²))], q1=a p1, T(t0)=a(t0)=1.
With negligible chi2, fixed mother and a nonzero boson mphi, the local
nonlinear Bose occupation solves:
  df_phi(q,t)/dt = A(q,t)*(1+f_phi(q,t))
  f_phi(q,t) = [1+fseed(q)]*exp(Integral A dt)-1.
A is independently obtained by exact Lorentz two-body decay geometry:
  A = (g1/gphi) Gamma M1²/(2 p* pphi Ephi)
      Integral_[E1min,E1max] f1(a sqrt(E1²-M1²)) dE1.
The q->0 limit is finite and has the explicit mother boost:
  A0 = (g1/gphi) Gamma (M1/mphi)^3
       f1(a sqrt((M1 Ephi*/mphi)^2-M1²)).
This is NOT the instantaneous FD expression at the decreasing
radiation temperature T=1/a.

No daughter chi2 blocking/inverse decay; no mother depletion, no SM
Q/gstarS, no full BEC zero-mode, no self-consistent Friedmann equation.
Never treat this benchmark as a CLASS-ready frozen 86-eV DM spectrum.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.special import expit

from massive_daughter_decay_source import rest_kinematics


def scale_factor(t,t0=1e-4):
    if not math.isfinite(t0) or t0<=0 or not math.isfinite(t) or t<t0:
        raise ValueError("Require finite radiation cosmic time t>=t0>0")
    return math.sqrt(t/t0)


def collisionless_parent_fd_at_physical_energy(E,a,M):
    if not all(math.isfinite(x) for x in (E,a,M)) or a<=0 or M<=0 or E<M:
        raise ValueError("Require physical on-shell mother energy E>=M and a>0")
    q1=a*math.sqrt(max((E-M)*(E+M),0.))
    return float(expit(-math.hypot(q1,M)))


def gain_coefficient(q,t,t0=1e-4,M=.1,mphi=.001,mchi=0.,
                     g1=4,gphi=1,rest_width=1.):
    """Exactly reduce quantum vacuum mother gain to 1D E1 quadrature.

    q is COMOVING daughter momentum and p=q/a(t).
    The boson factor (1+fphi) is not included in this A coefficient.
    """
    if (not math.isfinite(q) or q<0 or
            not all(math.isfinite(x) for x in (g1,gphi,rest_width)) or
            g1<=0 or gphi<=0 or rest_width<=0 or mphi<=0):
        raise ValueError("Finite nonnegative q, positive phi mass/degeneracy/width required")
    kin=rest_kinematics(M,mphi,mchi)
    a=scale_factor(t,t0)
    epstar=kin["rest_phi_energy"]
    pstar=kin["rest_daughter_momentum"]
    if q==0:
        E1=M*epstar/mphi
        return float((g1/gphi)*rest_width*(M/mphi)**3*
            collisionless_parent_fd_at_physical_energy(E1,a,M))
    p=q/a
    Ep=math.hypot(p,mphi)
    denominator=Ep*epstar+p*pstar
    Emin=M*(epstar*epstar+p*p)/denominator
    # The width is computed WITHOUT subtracting close endpoints.
    span=2*M*p*pstar/(mphi*mphi)
    integrand=lambda u:collisionless_parent_fd_at_physical_energy(
        Emin+span*u,a,M)
    if span<1e-5:
        parent_moment=span*integrand(.5)
    else:
        parent_moment=span*quad(integrand,0.,1.,epsabs=1e-13,
                                 epsrel=1e-11,limit=150)[0]
    return float((g1/gphi)*rest_width*M*M/
                 (2*pstar*p*Ep)*parent_moment)


def frozen_parent_rd_occupation(q,t=1.4e-4,t0=1e-4,
                                M=.1,mphi=.001,mchi=0.,
                                g1=4,gphi=1,rest_width=1.,seed=1e-12):
    if (not math.isfinite(seed) or seed<0 or
            not math.isfinite(t) or t<t0):
        raise ValueError("Require nonnegative finite seed and time t>=t0")
    A_integral=quad(lambda s:gain_coefficient(
        q,s,t0,M,mphi,mchi,g1,gphi,rest_width),
        t0,t,epsabs=1e-10,epsrel=1e-10,limit=150)[0]
    fseed=seed*math.exp(-q)
    return float((1+fseed)*math.expm1(A_integral)+fseed)


def frozen_parent_rd_ir_number(q_cut=1e-4,t=1.4e-4,t0=1e-4,
                               M=.1,mphi=.001,mchi=0.,
                               g1=4,gphi=1,rest_width=1.,seed=1e-12):
    if not math.isfinite(q_cut) or q_cut<=0:
        raise ValueError("Positive finite qcut required")
    integral=quad(
        lambda q:q*q*frozen_parent_rd_occupation(
            q,t,t0,M,mphi,mchi,g1,gphi,rest_width,seed),
        0,q_cut,epsabs=1e-19,epsrel=2e-9,limit=160)[0]
    return float(gphi/(2*math.pi**2)*integral)


def evaluate(t0=1e-4,end_t=1.4e-4,qmin=1e-4):
    result={
        "physical_model":"thermal FD mother frozen in COMOVING coordinates, collisionless RD a=sqrt(t/t0), massive phi; fixed chi1 reservoir and empty chi2 approximation",
        "cosmic_time_initial":t0,"cosmic_time_final":end_t,
        "scale_factor_final":scale_factor(end_t,t0),
        "qmin_comoving":qmin,
        "M1":.1,"mphi":.001,"mchi2":0.,
        "g1":4,"gphi":1,"rest_width":1.,
        "initial_vacuum_A_qmin":gain_coefficient(qmin,t0,t0),
        "redshifted_final_A_qmin":gain_coefficient(qmin,end_t,t0),
        "fphi_at_qzero":frozen_parent_rd_occupation(0.,end_t,t0),
        "fphi_at_qmin":frozen_parent_rd_occupation(qmin,end_t,t0),
        "below_qmin_comoving_number":frozen_parent_rd_ir_number(qmin,end_t,t0),
        "time_integrated_Bose_gain_not_instantaneous_FD_at_T_over_a":True,
        "short_time_frozen_parent_approximation_only":True,
        "macroscopic_condensate_or_real_86eV_FRW_relic_or_CLASS_validated":False
    }
    return result


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",default="outputs/frozen_collisionless_rd_phi.json")
    args=parser.parse_args()
    r=evaluate()
    path=Path(args.output);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(r,indent=2),encoding="utf-8")
    print(json.dumps(r,indent=2))
