"""Exact collisionless massive chi1 FRW background in comoving q=p/T.

For constant g_star_S, T~a^-1, x=M1/T. A freely propagating particle
has constant q=p/T; after decoupling at x0 the mother distribution is
 f1(q,x)=f1(q,x0)=1/(1+exp[sqrt(q^2+x0^2)]).
It must NOT be continuously reset to FD(sqrt(q^2+x^2)) without an
explicit external rethermalization/number-changing collision term.

Number n/T^3 is constant. Dimensionless rho(T)/T^4 and P(T)/T^4
obey x*d(rho_bar)/dx = rho_bar-3*P_bar, equivalent to physical FRW
 d rho/d ln a + 3 (rho+P)=0.
This is an EXACT collisionless background analytic/numerical benchmark
for a future closed FRW daughter+mother kinetic solver, not a full
collision-producing cosmological dark-matter scenario.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from scipy.integrate import quad
from scipy.special import expit


def frozen_mother_background(x,x_dec=.05,g1=4):
    if x<=0 or x_dec<=0 or g1<=0:
        raise ValueError("Require positive mass/temperature ratios and degeneracy")
    const=g1/(2*math.pi**2)
    def f(q):
        return float(expit(-math.hypot(q,x_dec)))
    def moment(func):
        return float(quad(lambda q:func(q)*f(q),0,math.inf,
                          epsabs=1e-12,epsrel=1e-11,limit=300)[0])
    n=const*moment(lambda q:q*q)
    rho=const*moment(lambda q:q*q*math.hypot(q,x))
    pressure=const/3*moment(lambda q:q**4/math.hypot(q,x))
    rho_derivative=const*moment(lambda q:q*q*x/math.hypot(q,x))
    # The instantaneous thermal FD formula is a DIFFERENT physics model:
    # a source externally replenishes/removes/rethermalizes mother states.
    instantaneous=const*quad(
        lambda q:q*q*float(expit(-math.hypot(q,x))),
        0,math.inf,epsabs=1e-12,epsrel=1e-11,limit=300)[0]
    continuity=(x*rho_derivative-rho+3*pressure)
    eos=pressure/rho
    return {
        "x_M1_over_T":float(x),
        "initial_decoupling_x":float(x_dec),
        "comoving_mother_number_over_T3":n,
        "dimensionless_rho_over_T4":rho,
        "dimensionless_pressure_over_T4":pressure,
        "equation_of_state_P_over_rho":eos,
        "d_dimensionless_rho_dx":rho_derivative,
        "FRW_collisionless_continuity_residual":continuity,
        "physical_d_ln_rho_by_d_ln_a":-4+x*rho_derivative/rho,
        "FRW_expected_d_ln_rho_by_d_ln_a":-3*(1+eos),
        "instantaneous_thermal_mother_n_over_T3":instantaneous,
        "instantaneous_thermal_to_frozen_number_ratio":instantaneous/n,
        "not_a_mother_rethermalization_model":True,
    }


def run():
    x0=.05
    xs=(.05,.1,.5,1.,2.,5.,10.,20.)
    cases=[frozen_mother_background(x,x_dec=x0) for x in xs]
    numbers=[r["comoving_mother_number_over_T3"] for r in cases]
    n0=numbers[0]
    assert max(abs(n/n0-1) for n in numbers)<1e-12
    assert max(abs(r["FRW_collisionless_continuity_residual"])
               for r in cases)<1e-11
    assert all(a["instantaneous_thermal_to_frozen_number_ratio"]>
               b["instantaneous_thermal_to_frozen_number_ratio"]
               for a,b in zip(cases,cases[1:]))
    return {
        "scope":"Collisionless FRW Liouville background for freely streaming initially thermal massive mother",
        "initial_decoupling_M1_over_T":x0,
        "assumption":"gstarS fixed, T proportional 1/a, massive constant M1, no collisions",
        "cases":cases,
        "important_result":"Resetting massive parent FD to instantaneous temperature is NOT collisionless evolution and implies external physics. Frozen q=p/T distribution conserves n/T³ and exactly obeys FRW energy continuity.",
        "warning":"Does not include 1<->2 reactions, stimulated Bose spectrum, SM energy transfer or CLASS. A later coupled closed FRW model needs an energy-consistent Liouville operator plus collisions."
    }


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--output",default="outputs/frw_collisionless_mother.json")
    o=p.parse_args()
    r=run()
    target=Path(o.output)
    target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(json.dumps(r,indent=2),encoding="utf-8")
    print(json.dumps(r,indent=2))
