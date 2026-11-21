"""Collisionless massive parent during SM entropy release: exact FRW test.

Physical SM entropy conservation: a³*gstarS(T)*T³=constant.
For x=M1/T and s(x)=[gstarS(x)/gstarS(x_dec)]^(1/3),
the invariant free-streaming momentum is q_c=p/[T*s(x)], whereas
p/T=s(x)*q_c. This equals q=p/T only if gstarS is constant.

For conserved f(q_c), n/T³ ∝ s³ and physical n*a³ is constant.
With alpha=d ln gstarS/d ln x,
  d ln a / d ln x = 1-alpha/3,
  x d(rho/T⁴)/dx = (1+alpha)*rho/T⁴+(alpha-3)*P/T⁴,
  d ln rho_physical / d ln a = -3(1+P/rho).

The gstarS(x) below is an ARTIFICIAL SM entropy-step test function.
Not a physical particle-theory/temperature table or connected 1<->2
production scenario.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from scipy.integrate import quad
from scipy.special import expit


class SmoothEntropyStep:
    def __init__(self,g_early=55.,g_late=3.909,
                 x_transition=1.,power=4.):
        if not (g_early>=g_late>0 and x_transition>0 and power>0):
            raise ValueError("positive entropy steps: early>=late>0")
        self.g_early=g_early
        self.g_late=g_late
        self.x_transition=x_transition
        self.power=power

    def value(self,x):
        if x<=0:raise ValueError("x=M1/T must be positive")
        u=(x/self.x_transition)**self.power
        return self.g_late+(self.g_early-self.g_late)/(1+u)

    def alpha_d_ln_gstar_d_ln_x(self,x):
        if x<=0:raise ValueError("positive x required")
        u=(x/self.x_transition)**self.power
        return (-(self.g_early-self.g_late)*self.power*u/
                ((1+u)**2*self.value(x)))


def frozen_entropy_mother(x,x_dec=.05,entropy=None,g1=4):
    """Independent 1D continuum FD source and entropy-dilution moments."""
    if x<=0 or x_dec<=0 or g1<=0:
        raise ValueError("positive x,x_dec,g1 required")
    entropy=entropy if entropy is not None else SmoothEntropyStep()
    g_dec=entropy.value(x_dec)
    g_x=entropy.value(x)
    alpha=entropy.alpha_d_ln_gstar_d_ln_x(x)
    s=(g_x/g_dec)**(1./3)
    pref=g1/(2*math.pi**2)
    def occupation(y):
        return float(expit(-math.hypot(y,x_dec)))
    def moment(integrand):
        return float(quad(lambda y:integrand(y)*occupation(y),
                          0,math.inf,epsabs=1e-12,epsrel=1e-11,
                          limit=300)[0])
    n0=pref*moment(lambda y:y*y)
    e=lambda y:math.hypot(s*y,x)
    rho=pref*s**3*moment(lambda y:y*y*e(y))
    pressure=pref*s**5/3*moment(lambda y:y**4/e(y))
    n=pref*s**3*moment(lambda y:y*y)
    # The derivative includes changes to both the SM entropy scaling
    # s(x) and the massive on-shell dispersion E=sqrt(s² y²+x²).
    x_drho=pref*s**3*moment(
        lambda y:y*y*(x*x+(alpha/3)*s*s*y*y)/e(y)) + alpha*rho
    from_pressure=(1+alpha)*rho+(alpha-3)*pressure
    a_ratio=(x/x_dec)/s
    expansion_dln_a_dln_x=1-alpha/3
    physical_slope=(-4+x_drho/rho)/expansion_dln_a_dln_x
    return {
        "x_M1_over_T":x,"decoupling_x":x_dec,
        "gstarS_current":g_x,"gstarS_decoupling":g_dec,
        "entropy_dln_gstar_dln_x":alpha,
        "invariant_q_definition":"q=p/[T*(gstarS(x)/gstarS(x_dec))^(1/3)]",
        "p_over_T_to_comoving_q_scale":s,
        "a_over_a_decoupling":a_ratio,
        "d_ln_a_d_ln_x":expansion_dln_a_dln_x,
        "initial_mother_n_over_Tdec3":n0,
        "current_mother_n_over_T3":n,
        "gstarS_number_ratio":g_x/g_dec,
        "comoving_number_charge":n*(x_dec/x)**3*a_ratio**3,
        "rho_over_T4":rho,"pressure_over_T4":pressure,
        "equation_of_state_P_over_rho":pressure/rho,
        "x_drho_over_T4_dx_from_integral":x_drho,
        "x_drho_over_T4_dx_from_pressure":from_pressure,
        "entropy_aware_energy_identity_residual":x_drho-from_pressure,
        "physical_d_ln_rho_d_ln_a":physical_slope,
        "continuity_expected_slope":-3*(1+pressure/rho),
        "collisionless":True,
        "gstar_history_only_toy":True,
    }


def run():
    entropy=SmoothEntropyStep()
    x0=.05
    xs=(.05,.2,.5,1.,2.,5.,10.,20.)
    cases=[frozen_entropy_mother(x,x_dec=x0,entropy=entropy) for x in xs]
    original=cases[0]["initial_mother_n_over_Tdec3"]
    for case in cases:
        assert abs(case["comoving_number_charge"]/original-1)<1e-11,case
        assert abs(case["entropy_aware_energy_identity_residual"])<1e-10,case
        assert abs(case["physical_d_ln_rho_d_ln_a"]-
                   case["continuity_expected_slope"])<1e-10,case
    return {
      "scope":"Collisionless massive mother + synthetic entropy release, no 1<->2 collisions",
      "entropy_model":"Smooth 55 -> 3.909 effective SM entropy degrees across M1/T~1; merely a test function, not a physical SM table",
      "initial_x":x0,
      "cases":cases,
      "scientific_result":"Constant p/T fails across SM entropy release; conserved collisionless variable is p/[T*gstarS^(1/3)]. Entropy-aware background exactly conserves n*a³ and satisfies rho+P FRW continuity.",
      "not_a_physical_DM_relic_prediction":True
    }


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--output",default="outputs/frw_entropy_release_mother.json")
    o=p.parse_args()
    result=run()
    path=Path(o.output)
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(result,indent=2),encoding="utf-8")
    print(json.dumps(result,indent=2))
