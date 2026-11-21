"""Charge-constrained e-/e+ finite-volume heat susceptibility.

References the ideal paired EM equilibrium EOS but keeps net electron charge
density FIXED during an instantaneous heat-transfer collision step.
No e+e- <-> gamma gamma reaction rates or expansion/eta(T) are implied.
"""
from __future__ import annotations

import argparse
import json
import math

import numpy as np
from scipy.optimize import brentq
from scipy.special import expit

from kinetics.em_pair_thermodynamics import (
    M_E_GEV, N_PHOTON_OVER_T3, PHOTON_RHO_OVER_T4,
    quadrature, charge_neutral_mu,
)


def _charge_from_y(x, y, q, weights):
    eps=np.sqrt(q*q+x*x)
    difference=np.sinh(y)/(np.cosh(eps)+np.cosh(y))
    return float(np.dot(weights*q*q,difference)/math.pi**2)


def solve_mu_fixed_charge(T, charge_density, *, me=M_E_GEV, n=256):
    """Return y=mu_e-/T for a constant *physical volume* net charge nQ."""
    if not (math.isfinite(T) and math.isfinite(me) and T>0 and me>0
            and .02<=T/me<=100 and math.isfinite(charge_density)
            and charge_density>=0):
        raise ValueError("finite nonnegative charge and supported T/me required")
    x=me/T
    q,w=quadrature(n)
    target=charge_density/T**3
    return brentq(lambda y:_charge_from_y(x,y,q,w)-target,
                  0.,x+100.,xtol=1e-13)


def pair_energy(T, y, *, me=M_E_GEV, n=256):
    """e-/e+ rho [GeV^4] at fixed dimensionless chemical potential y."""
    q,w=quadrature(n)
    eps=np.sqrt(q*q+(me/T)**2)
    return T**4/math.pi**2*float(np.dot(
        w*q*q,eps*(expit(y-eps)+expit(-y-eps))))


def charge_constrained_heat(T, *, me=M_E_GEV, eta=6e-10, n=256):
    """Positive Schur-complement thermal capacity at fixed NET charge.

    Evaluates C_Q = rho_T|y - (rho_y * n_T|y)/n_y at fixed volume.
    eta is used ONLY to label reference nQ(T_initial), and is NOT held
    constant during a temperature derivative (because n_gamma~T^3).
    """
    if not (math.isfinite(T) and math.isfinite(me) and T>0 and me>0
            and .02<=T/me<=100 and math.isfinite(eta) and eta>0):
        raise ValueError("supported T/me and positive illustrative eta required")
    y=charge_neutral_mu(T,m=me,eta=eta,n=n)["mu_e_minus_over_T"]
    x=me/T
    q,w=quadrature(n)
    eps=np.sqrt(q*q+x*x)
    fminus=expit(y-eps)
    fplus=expit(-y-eps)
    a=fminus*(1-fminus)
    b=fplus*(1-fplus)
    common=1/math.pi**2
    B=common*float(np.dot(w*q*q,eps**2*(a+b)))
    D=common*float(np.dot(w*q*q,eps*(a-b)))
    S=common*float(np.dot(w*q*q,a+b))
    Cfixed=B-D*D/S
    if Cfixed<0:
        raise ArithmeticError("violated covariance Cauchy-Schwarz bound")
    return {
       "T_over_me":T/me, "eta_for_initial_charge":eta,
       "mu_e_minus_over_T":y,
       "initial_charge_density_GeV3":eta*N_PHOTON_OVER_T3*T**3,
       "electron_pair_C_at_fixed_y_over_T3":B,
       "electron_pair_C_at_fixed_net_charge_over_T3":Cfixed,
       "electron_capacity_ratio_CQ_to_Cy":Cfixed/B,
       "EM_C_at_fixed_net_charge_over_T3":Cfixed+4*PHOTON_RHO_OVER_T4,
       "thermodynamic_charges":"net electron minus positron; photons have zero electric charge",
       "not_annihilation_or_cosmological_evolution":True
    }


def independent_finite_difference(T, *, me=M_E_GEV, eta=6e-10,
                                  h=1e-4, n=256):
    """Independent FD energy derivative at fixed physical volume nQ."""
    if not 0<h<.001 or T*(1-h)/me<.02 or T*(1+h)/me>100:
        raise ValueError("stay inside supported range")
    initial_charge=eta*N_PHOTON_OVER_T3*T**3
    Tp,Tm=T*(1+h),T*(1-h)
    yp=solve_mu_fixed_charge(Tp,initial_charge,me=me,n=n)
    ym=solve_mu_fixed_charge(Tm,initial_charge,me=me,n=n)
    Cp=(pair_energy(Tp,yp,me=me,n=n)-pair_energy(Tm,ym,me=me,n=n))/(Tp-Tm)/T**3
    return Cp


if __name__ == "__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--T-over-me",type=float,default=.025)
    p.add_argument("--eta",type=float,default=6e-10)
    a=p.parse_args()
    T=a.T_over_me*M_E_GEV
    c=charge_constrained_heat(T,eta=a.eta)
    c["independent_charge_fixed_fd_C_over_T3"]=independent_finite_difference(T,eta=a.eta)
    print(json.dumps(c,indent=2))
