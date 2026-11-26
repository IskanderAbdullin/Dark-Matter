"""Finite equilibrium e-/e+ + photon thermodynamic *reference*, not a plasma solver.

Ideal massive Dirac pairs in chemical equilibrium (mu_e+=-mu_e-), photons
with zero chemical potential. The mu=0 heat capacity assumes FAST pair
chemistry; no e+e- <-> gamma gamma rate, neutrino coupling, baryons,
finite-T QED or cosmological evolution is modeled. Units GeV.
"""
from __future__ import annotations
import argparse
import json
import math
from functools import lru_cache

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import brentq
from scipy.special import expit, kv, zeta

M_E_GEV = 0.00051099895
PHOTON_RHO_OVER_T4 = math.pi**2/15
N_PHOTON_OVER_T3 = 2*zeta(3,1)/math.pi**2


@lru_cache(maxsize=16)
def quadrature(n):
    if not isinstance(n,int) or n<16:
        raise ValueError("Gauss order must be integer >= 16")
    x,w=leggauss(n)
    return 24*(x+1),24*w  # dimensionless p/T: 0..48


def equilibrium(T, *, m=M_E_GEV, n=256):
    """mu_e±=mu_gamma=0; rho, P, C at fixed m, mu=0, in thermal units.

    C_em/T^3 = d(rho_e+rho_gamma)/dT /T^3, assuming fast e+e-
    chemical equilibration (an assumption, not a derived collision rate).
    """
    if not (math.isfinite(T) and math.isfinite(m) and T>0 and m>0
            and .02<=T/m<=100):
        raise ValueError("supported thermodynamic range .02<=T/m<=100")
    x=m/T
    q,w=quadrature(n)
    eps=np.sqrt(q*q+x*x)
    f=expit(-eps)
    common=4/(2*math.pi**2)
    rho=common*np.dot(w*q*q,eps*f)
    P=common/3*np.dot(w*q**4,f/eps)
    Ce=common*np.dot(w*q*q,eps*eps*f*(1-f))
    pair_n=common*np.dot(w*q*q,f)
    return {
      "T_over_me":T/m, "me_GeV":m,
      "rho_pair_over_T4":float(rho), "P_pair_over_T4":float(P),
      "n_pair_over_T3":float(pair_n),
      "C_pair_at_fixed_mu0_over_T3":float(Ce),
      "C_gamma_plus_pair_over_T3":float(Ce+4*PHOTON_RHO_OVER_T4),
      "gstar_rho_gamma_plus_pair":float(30*(rho+PHOTON_RHO_OVER_T4)/math.pi**2),
      "gstar_s_gamma_plus_pair":float(45*(rho+P+4*PHOTON_RHO_OVER_T4/3)/(2*math.pi**2)),
      "fast_chemical_equilibrium_assumed":True,
      "not_charge_neutral_mu0_at_late_times":True
    }


def bessel_pair(x, *, terms=400):
    """Independent fugacity/Bessel series for pressure and energy density."""
    if not (math.isfinite(x) and x>0 and isinstance(terms,int) and terms>=1):
        raise ValueError
    pr=rho=0.
    for n in range(1,terms+1):
        sg=1 if n%2 else -1
        pr+=sg*kv(2,n*x)/n**2
        rho+=sg*(3*kv(2,n*x)/n**2+x*kv(1,n*x)/n)
    factor=4*x*x/(2*math.pi**2)
    return factor*pr, factor*rho


def charge_neutral_mu(T, *, m=M_E_GEV, eta=6e-10, n=256):
    """Solve illustrative n_e- - n_e+ = eta*n_gamma, at fixed T and eta.

    Assumes massless equilibrium photons, homogeneous positive background
    charge and fast pair chemical equilibration. NOT post-annihilation
    baryon kinetic evolution or an inferred cosmological eta history.
    """
    if not (math.isfinite(eta) and eta>0 and math.isfinite(T)
            and .02<=T/m<=100):
        raise ValueError("finite positive eta, supported T/m")
    q,w=quadrature(n)
    x=m/T
    eps=np.sqrt(q*q+x*x)
    target=eta*N_PHOTON_OVER_T3
    def charge(y):
        # Stable small-y FD difference, avoids f_minus-f_plus cancellation.
        delta=np.sinh(y)/(np.cosh(eps)+np.cosh(y))
        return float(np.dot(w*q*q,delta)/math.pi**2)
    y=brentq(lambda y:charge(y)-target,0.,x+100,xtol=1e-13)
    ne_minus=float(np.dot(w*q*q,expit(y-eps))/math.pi**2)
    ne_plus=float(np.dot(w*q*q,expit(-y-eps))/math.pi**2)
    return {"T_over_me":T/m,"eta_assumed":eta,
            "mu_e_minus_over_T":y,"mu_e_plus_over_T":-y,
            "net_electron_charge_over_T3":charge(y),
            "target_charge_over_T3":target,
            "n_e_minus_over_T3":ne_minus,"n_e_plus_over_T3":ne_plus,
            "positron_to_electron_density_ratio":ne_plus/ne_minus,
            "not_a_dynamical_annihilation_solution":True}


if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--T-over-me",type=float,default=1.)
    ap.add_argument("--eta",type=float,default=6e-10)
    op=ap.parse_args()
    T=op.T_over_me*M_E_GEV
    print(json.dumps({"mu0_thermo":equilibrium(T),
                      "neutral_charge":charge_neutral_mu(T,eta=op.eta)},indent=2))
