"""Tree-level MICROPHYSICS PILOT, not a thermal kinetic-equilibration rate.

Natural units hbar=c=k_B=1. Masses and temperature in GeV, y dimensionless,
G in GeV^-2, widths / frequencies in GeV. The electron/positron bath
assumes massless e, zero chemical potential and ignores final-state
Pauli blocking, neutrino bath channels, thermal screening and DM recoil.
It is NOT a physically calibrated replacement for gamma_el in BEST+SM.

Lagrangian benchmark: y phi (bar chi2 chi1 + h.c.); contact operator
G (bar chi1 gamma_mu chi1)(bar e gamma^mu e), below EWSB.
Z2: chi1,phi odd, chi2 and SM even; phi is lightest odd particle.
"""
from __future__ import annotations

import json
import math
from scipy.integrate import quad

PI = math.pi
MPL_GEV = 1.2209e19  # non-reduced Planck mass, radiation-era approximation


def _finite_nonnegative(**values):
    if not all(math.isfinite(float(v)) and float(v) >= 0 for v in values.values()):
        raise ValueError("Expected finite nonnegative physical parameters")


def _decay_masses(m1, m2, mphi):
    _finite_nonnegative(m1=m1, m2=m2, mphi=mphi)
    if m1 <= 0:
        raise ValueError("Parent mass must be positive")


def kallen_sqrt_decay(m1, m2, mphi):
    """sqrt[lambda(m1^2,m2^2,mphi^2)]; zero at/below threshold."""
    _decay_masses(m1, m2, mphi)
    if m1 <= m2 + mphi:
        return 0.0
    return math.sqrt((m1 - m2 - mphi) * (m1 + m2 + mphi) *
                     (m1 - m2 + mphi) * (m1 + m2 - mphi))


def scalar_yukawa_spin_averaged_m2(m1, m2, mphi, y):
    """One particle channel, averaged over its two parent spins."""
    _decay_masses(m1, m2, mphi)
    _finite_nonnegative(y=y)
    return y*y*((m1 + m2)**2 - mphi*mphi)


def scalar_yukawa_width(m1, m2, mphi, y):
    """Tree-level partial rest width chi1 -> chi2 + real phi in GeV.

    Gamma = [y^2 ((m1+m2)^2-mphi^2)] sqrt(lambda)/(16 pi m1^3).
    No medium effects, thermal masses, backreaction or other channels.
    """
    amp = scalar_yukawa_spin_averaged_m2(m1, m2, mphi, y)
    return amp * kallen_sqrt_decay(m1, m2, mphi)/(16*PI*m1**3)


def best_spin_summed_m2(m1, m2, mphi, y, parent_dof=4):
    """BEST sum convention from existing decay benchmarks.

    BEST protocol uses g_parent * 16 pi m1^3 Gamma/sqrt(lambda).
    For a Dirac particle+antiparticle (g_parent=4) this equals
    4*y^2*((m1+m2)^2-mphi^2); other g need physical channel review.
    Return the matrix element even when decay kinematics is closed.
    """
    _finite_nonnegative(parent_dof=parent_dof)
    if parent_dof <= 0:
        raise ValueError("parent_dof must be positive")
    return parent_dof*scalar_yukawa_spin_averaged_m2(m1,m2,mphi,y)


def width_from_best_spin_sum(m1,m2,mphi,amplitude,parent_dof=4):
    """Inverse width only in the same BEST normalization convention."""
    _decay_masses(m1,m2,mphi)
    _finite_nonnegative(amplitude=amplitude,parent_dof=parent_dof)
    if parent_dof <= 0:
        raise ValueError("parent_dof must be positive")
    return (amplitude/parent_dof*kallen_sqrt_decay(m1,m2,mphi)/
            (16*PI*m1**3))


def vector_contact_m2(s,t,mchi,G):
    """Spin-averaged |M|^2 for chi e -> chi e; e massless, vector contact.

    u=2*mchi^2-s-t. Two incoming spin-1/2 particles are averaged.
    Physical t: -(s-mchi^2)^2/s <= t <= 0.
    """
    _finite_nonnegative(s=s,mchi=mchi,G=G)
    if mchi <= 0 or s <= mchi*mchi:
        raise ValueError("Physical s > mchi^2 > 0 required")
    lo=-(s-mchi*mchi)**2/s
    if not math.isfinite(t) or t < lo-1e-12*abs(lo) or t > 1e-12*abs(lo):
        raise ValueError("t outside elastic scattering range")
    u=2*mchi*mchi-s-t
    return 2*G*G*((s-mchi*mchi)**2+(u-mchi*mchi)**2+2*mchi*mchi*t)


def vector_contact_dsdt(s,t,mchi,G):
    return vector_contact_m2(s,t,mchi,G)/(16*PI*(s-mchi*mchi)**2)


def sigma_total_massless_e(s,mchi,G):
    """Tree-level unblocked total cross section in GeV^-2."""
    _finite_nonnegative(s=s,mchi=mchi,G=G)
    if mchi <= 0 or s <= mchi*mchi:
        raise ValueError("Physical s > mchi^2 > 0 required")
    D=s-mchi*mchi
    return G*G*D*D*(D*D+3*s*s)/(24*PI*s**3)


def sigma_transport_massless_e(s,mchi,G):
    """Integral (1-cos theta_cm) d sigma, massless electron, in GeV^-2.

    NOT the kinetic/energy relaxation rate for a moving mother particle.
    """
    _finite_nonnegative(s=s,mchi=mchi,G=G)
    if mchi <= 0 or s <= mchi*mchi:
        raise ValueError("Physical s > mchi^2 > 0 required")
    D=s-mchi*mchi
    return G*G*D*D*(3*D*D+4*s*s)/(48*PI*s**3)


def bath_transport_frequency_chi_at_rest(mchi,T,G,g_bath=4):
    """Unblocked chi-at-rest transport collision FREQUENCY in GeV.

    nu = g/(2*pi^2) int dk k^2 f_FD(k/T) sigma_T(mchi^2+2*mchi*k).
    Counts only massless e-/e+ (g=4); no final Pauli, full bath
    collision kernel, drag/diffusion, neutrino or hard SM channels.
    It MUST NOT be passed as gamma_el to BEST's BGK without validation.
    """
    _finite_nonnegative(mchi=mchi,T=T,G=G,g_bath=g_bath)
    if mchi <= 0 or T <= 0 or g_bath <= 0:
        raise ValueError("mchi,T,g_bath must be positive")
    def integrand(x):
        # Positive, stable Fermi-Dirac; x>=0 on quadrature interval.
        f=math.exp(-x)/(1+math.exp(-x))
        s=mchi*mchi+2*mchi*T*x
        return x*x*f*sigma_transport_massless_e(s,mchi,G)
    val,err=quad(integrand,0,45,epsabs=1e-13,epsrel=1e-10,limit=200)
    return g_bath*T**3/(2*PI**2)*val


def bath_transport_frequency_heavy_limit(T,G,g_bath=4):
    """M/T -> infinity reference, nu ~15*g*zeta(5)*G^2*T^5/(4*pi^3)."""
    _finite_nonnegative(T=T,G=G,g_bath=g_bath)
    if T <= 0 or g_bath <= 0:
        raise ValueError("T and g_bath must be positive")
    return 15*g_bath*1.0369277551433699*G*G*T**5/(4*PI**3)


def radiation_H(T,gstar_rho):
    """H [GeV] = 1.66 sqrt(gstar_rho) T^2/M_Pl in radiation domination."""
    _finite_nonnegative(T=T,gstar_rho=gstar_rho)
    if T <= 0 or gstar_rho <= 0:
        raise ValueError("T and gstar_rho must be positive")
    return 1.66*math.sqrt(gstar_rho)*T*T/MPL_GEV


def exploratory_snapshot(*,m1=0.1,m2=0.,mphi=86e-9,T=0.1,
                         yukawa=1e-9,eft_scale=1.,coefficient=1.,
                         gstar_rho=55.):
    """Dimensionful illustrative scan; NO physical model fit or BEC claim."""
    _finite_nonnegative(eft_scale=eft_scale,coefficient=coefficient)
    if eft_scale <= 0: raise ValueError("EFT scale must be positive")
    width=scalar_yukawa_width(m1,m2,mphi,yukawa)
    G=coefficient/eft_scale**2
    nu=bath_transport_frequency_chi_at_rest(m1,T,G)
    H=radiation_H(T,gstar_rho)
    return {
      "units":"GeV for masses,T,width,nu,H; G in GeV^-2",
      "parameters":{"m1":m1,"m2":m2,"mphi":mphi,"T":T,
                    "yukawa":yukawa,"Lambda":eft_scale,"C":coefficient,
                    "gstar_rho":gstar_rho},
      "tree_vacuum_Gamma_rest":width,
      "tree_BEST_spin_summed_m2":best_spin_summed_m2(m1,m2,mphi,yukawa),
      "unblocked_e_only_rest_transport_collision_frequency":nu,
      "radiation_H_approx":H,
      "vacuum_Gamma_rest_over_H":width/H,
      "rest_transport_frequency_over_H_NOT_gamma_el_over_H":nu/H,
      "eft_warning":"Contact EFT valid only for momenta well below mediator scale; "
                    "electroweak completion and neutrino channels omitted.",
      "no_thermal_Gamma_or_bgk_calibration":True,
      "no_frozen_relic_or_BEC_or_CLASS":True,
    }


if __name__=="__main__":
    print(json.dumps(exploratory_snapshot(),indent=2))
