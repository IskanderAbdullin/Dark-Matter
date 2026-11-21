"""Analytic vacuum-source infrared bound for the on-shell 1<->2 model.

For mass M1/T=m, massless phi/chi2, thermal FD mother, EMPTY daughters,
and c=m^2/4, the exact differential vacuum-production source is

 dR_vac/dq_phi = g1*m/(2*pi^2) * log(1+exp[-q_phi-c/q_phi]).

It follows from the on-shell domain 4q_phi q_chi2>=m^2, integrating
q_chi2 from c/q_phi to infinity. For any epsilon>0:

 R(q_phi<epsilon)/R(total)
 <= epsilon*E2(c/epsilon)/I(total)
 <= epsilon*exp(-c/epsilon)/I(total),

where E2 is the order-two exponential integral and
I(total)=integral_0^inf log(1+exp[-q-c/q]) dq.

This is a rigorous bound ONLY on the empty-daughter vacuum source.
It is NOT a bound on the stimulated NONLINEAR Bose population, where
the 1+b factor can change the infrared flux dramatically.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.special import expn

from onshell_collision_2d import vacuum_reference


def _log_fd_integrand(q,c):
    return 0. if q<=0 else float(np.logaddexp(0.,-q-c/q))


def thermal_vacuum_ir_fraction(epsilon,mass=.1,g1=4):
    if epsilon<=0 or mass<=0 or g1<=0:
        raise ValueError("positive q cutoff, mass and degeneracy required")
    c=mass*mass/4.
    pref=g1*mass/(2*math.pi**2)
    r_full=vacuum_reference(M1_T=mass,g1=g1)
    I_full=r_full/pref
    sharp=epsilon*float(expn(2,c/epsilon))
    # The integral may be ~1e-116: an absolute tolerance of 1e-20
    # allows QUADPACK to stop long before a rigorous inequality can
    # be checked. Scale epsabs to the known analytic upper bound.
    epsabs=min(1e-20,max(1e-300,sharp*1e-11))
    I_low=quad(lambda q:_log_fd_integrand(q,c),0.,epsilon,
               epsabs=epsabs,epsrel=1e-10,limit=300)[0]
    elementary=epsilon*math.exp(-c/epsilon)
    rel=I_low/I_full
    assert 0<=rel<=sharp/I_full*(1+1e-10)+1e-200
    assert sharp<=elementary*(1+1e-10)+1e-200
    return {
       "IR_threshold_q_over_T":float(epsilon),
       "M1_over_T":float(mass),
       "on_shell_c_over_T2":c,
       "total_vacuum_source_per_T3_Gamma":r_full,
       "actual_vacuum_fraction_below_threshold":float(rel),
       "sharper_exponential_integral_bound":float(sharp/I_full),
       "elementary_exponential_bound":float(elementary/I_full),
       "source_per_q_definition":"g1*m/(2pi²)*log(1+exp[-q-m²/(4q)])",
       "applicable_to_stimulated_nonlinear_source":False,
       "caveat":"Only initially empty daughters and fixed thermal parent bath; not an upper bound for Bose-stimulated occupation or late frozen DM."
    }


def run():
    q=(1e-5,1e-4,3e-4,5e-4,1e-3,2e-3,1e-2)
    cases=[thermal_vacuum_ir_fraction(x) for x in q]
    assert all(a["actual_vacuum_fraction_below_threshold"]<
               b["actual_vacuum_fraction_below_threshold"]
               for a,b in zip(cases,cases[1:]))
    return {
       "study":"rigorous EMPTY-daughter vacuum IR-tail source bound",
       "M1_over_T":.1,
       "source_model":"thermally fixed parent FD, massless daughters, 1<->2 on-shell",
       "exact_total_source":cases[0]["total_vacuum_source_per_T3_Gamma"],
       "cases":cases,
       "limits":"The observed nonlinear cold stimulated burst may violate these vacuum small-q tail fractions without violating the vacuum bound, due to Bose gain 1+f and inverse reactions.",
       "not_a_cosmological_freezeout":True,
    }


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--output",default="outputs/vacuum_ir_tail_bound.json")
    opts=p.parse_args()
    result=run()
    out=Path(opts.output);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,indent=2),encoding="utf-8")
    print(json.dumps(result,indent=2))
