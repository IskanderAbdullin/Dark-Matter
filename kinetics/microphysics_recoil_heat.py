"""Microscopically derived electron -> chi1 RECOIL HEAT reference, not BGK gamma.

For a chi1 at REST scattering on a massless electron/positron FD bath,
tree vector four-fermion contact, no incoming chi thermal distribution,
no outgoing chi Pauli, no screening or mediator propagator. Exact
two-body recoil k'=k+t/(2 M); final e blocking included.

Q_chi_per_chi = (g/(2pi^2)) int dk k^2 f_e(k)
  int_tmin^0 dt (d sigma/dt) [1-f_e(k')] [-t/(2M)].

Natural units: Q [GeV^2] (energy gained per unit natural time).
A positive rest-particle heating power cannot be used as a full
linearized drag/diffusion operator for arbitrary f_chi(p).
"""
from __future__ import annotations
import math
from scipy.integrate import quad
from microphysics_portal import vector_contact_dsdt


def _validate(M,T,G,g):
    if (not all(math.isfinite(x) for x in (M,T,G,g)) or
            M<=0 or T<=0 or G<0 or g<=0):
        raise ValueError("Finite M,T,g>0, G>=0 required")


def fermi_dimensionless(x):
    if x>=0:
        a=math.exp(-x)
        return a/(1+a)
    return 1/(1+math.exp(x))


def energy_gain_per_incident_energy(k,M,G,T,blocking=True):
    """int d sigma/dt * Delta E_chi (1-f_final), in GeV^-1.

    k initial e energy (GeV), M chi mass, t<0, k_final=k+t/(2M).
    """
    _validate(M,T,G,1.)
    if not math.isfinite(k) or k<0:
        raise ValueError("finite k>=0 required")
    if k==0 or G==0:
        return 0.
    s=M*M+2*M*k
    tmin=-(s-M*M)**2/s
    def at_z(z):
        t=tmin*z
        e_final=k+t/(2*M)
        if e_final<0:
            raise ArithmeticError("negative recoil electron energy")
        blocking_factor=(1-fermi_dimensionless(e_final/T)) if blocking else 1.
        return (-tmin)*vector_contact_dsdt(s,t,M,G)*(-t)/(2*M)*blocking_factor
    result,_=quad(at_z,0,1,epsabs=0,epsrel=1e-10,limit=100)
    return result


def recoil_heating_rate_chi_at_rest(M,T,G,g_bath=4,blocking=True):
    """Rest-chi e+e- heating Q per chi, GeV^2; not a thermal gamma_el."""
    _validate(M,T,G,g_bath)
    if G==0:
        return 0.
    def at_x(x):
        if x<=0:
            return 0.
        return (x*x*fermi_dimensionless(x)*
                energy_gain_per_incident_energy(x*T,M,G,T,blocking))
    integral,_=quad(at_x,0,36,epsabs=0,epsrel=1e-9,limit=120)
    return g_bath*T**3/(2*math.pi**2)*integral


def heavy_chi_recoil_asymptote(M,T,G,g_bath=4):
    """M/T -> infinity with final-state electron Pauli blocking included.

    Analytic low-energy d sigma_T ~ G^2 k^2/(3pi);
    int dx x^6 f(1-f) = 31 pi^6/42.
    Thus Q ~ 31*g*pi^3*G^2*T^7/(252*M).
    """
    _validate(M,T,G,g_bath)
    return 31*g_bath*math.pi**3*G*G*T**7/(252*M)
