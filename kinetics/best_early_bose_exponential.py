"""Analytical local Bose stimulation in a FROZEN parent bath, massive phi.

The exact initial vacuum differential source R'(q) defines
A(q) = [2pi²/(gphi*q²)] R'(q), with the finite-mass q->0 limit
derived independently in massive_daughter_decay_source.py.

IF mother chi1 remains fixed and chi2 occupation stays negligible,
the full local Bose factor evolves:
   df_phi(q)/dt = A(q) * [1+f_phi(q)].
Consequently
   f(q,t) = [1+f_seed(q)] * expm1(A(q)*t) + f_seed(q).

This is a CLOSED-FORM SHORT-TIME APPROXIMATION, not the actual closed
3-species quantum Boltzmann equation: once chi2 is occupied, inverse
decays and mother depletion matter. Finite q<qmin number is
 n = gphi/(2pi²) integral_0^qmin q² f(q,t)dq.

The above below-qmin finite-particle-domain number is NOT automatically
a macroscopic condensate. Independently compare this approximation
against actual BEST qmin occupation and massive zero-mode N0; inspect
the difference rather than equating the two.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import numpy as np
from scipy.integrate import quad

from massive_daughter_decay_source import (
    massive_daughter_vacuum_spectrum,
    zero_momentum_bose_occupation_source,
)


def initial_gain_coefficient(q,M=.1,mphi=.001,mchi=0.,
                             g1=4,gphi=1,rest_width=1.):
    Q=np.asarray(q,dtype=float)
    if np.any(~np.isfinite(Q)) or np.any(Q<0) or gphi<=0:
        raise ValueError("Nonnegative finite momentum and Bose degeneracy required")
    gain=np.zeros_like(Q)
    present=Q>0
    if np.any(present):
        vals=Q[present]
        source=massive_daughter_vacuum_spectrum(
            vals,M=M,mphi=mphi,mchi=mchi,g1=g1,rest_width=rest_width)
        gain[present]=2*np.pi**2*source/(gphi*vals*vals)
    if np.any(~present):
        gain[~present]=zero_momentum_bose_occupation_source(
            M=M,mphi=mphi,mchi=mchi,g1=g1,gphi=gphi,rest_width=rest_width)
    if Q.ndim==0:return float(gain)
    return gain


def frozen_bath_phi_occupation(q,t,M=.1,mphi=.001,mchi=0.,
                               g1=4,gphi=1,rest_width=1.,seed=1e-12):
    Q=np.asarray(q,dtype=float)
    if t<0 or seed<0 or not np.isfinite(t) or not np.isfinite(seed):
        raise ValueError("Positive time and nonnegative finite seed required")
    A=initial_gain_coefficient(Q,M,mphi,mchi,g1,gphi,rest_width)
    f0=seed*np.exp(-Q)
    predicted=(1+f0)*np.expm1(A*t)+f0
    if Q.ndim==0:return float(predicted)
    return predicted


def frozen_bath_ir_number(qcut=.0001,t=4e-5,M=.1,mphi=.001,mchi=0.,
                          g1=4,gphi=1,rest_width=1.,seed=1e-12):
    if qcut<=0 or not math.isfinite(qcut) or t<0:
        raise ValueError("Positive qcut and nonnegative t required")
    N=quad(lambda q:q*q*frozen_bath_phi_occupation(
        q,t,M,mphi,mchi,g1,gphi,rest_width,seed),
        0,qcut,epsabs=1e-19,epsrel=1e-10,limit=350)[0]
    return float(gphi/(2*np.pi**2)*N)


def evaluate(t=4e-5,qmin=1e-4):
    gain0=initial_gain_coefficient(0.)
    gainq=initial_gain_coefficient(qmin)
    predictedq=frozen_bath_phi_occupation(qmin,t)
    predictedN0=frozen_bath_ir_number(qmin,t)
    return {
        "mass_over_T":.1,"mphi_over_T":.001,
        "tGamma":t,"qmin":qmin,
        "independent_massive_vacuum_A0":gain0,
        "independent_massive_vacuum_A_qmin":gainq,
        "linear_initial_source_f_qmin":gainq*t,
        "nonlinear_frozen_bath_phi_f_qmin":predictedq,
        "nonlinear_frozen_bath_phi_q_below_qmin_number":predictedN0,
        "initial_bose_seed_shape":"1e-12 exp(-q)",
        "assumptions":"chi1 thermal reservoir fixed; chi2 remains empty; phi evolves only via Bose 1+f; no inverse decay",
        "has_no_explicit_condensate_equation":True,
        "not_self_consistent_FRW_or_CLASS":True
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",default="outputs/frozen_bath_massive_bose_early.json")
    p.add_argument("--t",type=float,default=4e-5)
    args=p.parse_args()
    r=evaluate(args.t)
    out=Path(args.output);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(r,indent=2),encoding="utf-8")
    print(json.dumps(r,indent=2))


if __name__=="__main__":
    main()
