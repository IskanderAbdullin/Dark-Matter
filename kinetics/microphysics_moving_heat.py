"""Exact unblocked energy moment of vector-contact chi+massless-e 2->2.

For a moving chi and isotropic initial FD electron bath, the outgoing
CM scattering azimuth averages exactly. The resulting lab energy kick
per differential angular cross section is proportional to (1-cos theta*),
hence the FULL outgoing angular energy moment reduces analytically to
the physical sigma_transport(s) times an exact boost factor.

Q(p) [GeV^2/chi] includes correct incoming Moller flux and recoil;
e final-state Pauli, chi final-state blocking, inverse thermal population
and plasma backreaction are NOT included. Not a gamma_el or full Boltzmann
kernel, even though this is a useful physical energy-moment benchmark.
"""
from __future__ import annotations
import math
from scipy.integrate import quad
from microphysics_portal import sigma_transport_massless_e
from microphysics_recoil_heat import recoil_heating_rate_chi_at_rest


def _inputs(M,p,T,G,g):
    if (not all(math.isfinite(z) for z in (M,p,T,G,g)) or
            M<=0 or p<0 or T<=0 or G<0 or g<=0):
        raise ValueError("Finite positive M,T,g and nonnegative p,G required")


def energy_kick_scattering_azimuth_average(M,p,k,mu):
    """Boost-exact coefficient: outgoing-angle integral is coeff*sigma_T.

    Echi=sqrt(M^2+p^2), s=M^2+2*k*(Echi-p*mu).
    coeff = -(Echi-gamma_cm*Echi_cm*) =
          k*[ k*(Echi-p*mu)+Echi*p*mu-p^2 ]/s.
    This is a signed energy in GeV.
    """
    if (not all(math.isfinite(z) for z in (M,p,k,mu)) or
            M<=0 or p<0 or k<=0 or not -1<=mu<=1):
        raise ValueError("M,k>0, p>=0 and mu in [-1,1] required")
    E=math.hypot(M,p)
    s=M*M+2*k*(E-p*mu)
    return k*(k*(E-p*mu)+E*p*mu-p*p)/s


def incoming_s_and_moller_flux(M,p,k,mu):
    if (not all(math.isfinite(z) for z in (M,p,k,mu)) or
            M<=0 or p<0 or k<=0 or not -1<=mu<=1):
        raise ValueError("Unphysical incoming 2->2 momenta")
    E=math.hypot(M,p)
    flux=1-p*mu/E
    s=M*M+2*k*E*flux
    return s,flux


def moving_chi_unblocked_energy_moment(M,p,T,G,g_bath=4):
    """Net lab chi energy transfer per time for fixed incoming chi momentum.

    g/(2 pi^2) int dk k^2 fFD(k/T) * (1/2) int dmu v_M
                                         sigma_T(s) E_kick_avg.
    """
    _inputs(M,p,T,G,g_bath)
    if G==0:return 0.
    if p==0:
        return recoil_heating_rate_chi_at_rest(
            M,T,G,g_bath,blocking=False)
    def at_x(x):
        if x<=0:return 0.
        k=x*T
        def at_mu(mu):
            s,flux=incoming_s_and_moller_flux(M,p,k,mu)
            kick=energy_kick_scattering_azimuth_average(M,p,k,mu)
            return flux*sigma_transport_massless_e(s,M,G)*kick
        angular,_=quad(at_mu,-1,1,epsabs=0,epsrel=1e-9,limit=140)
        f=math.exp(-x)/(1+math.exp(-x))
        return x*x*f*angular/2
    result,_=quad(at_x,0,40,epsabs=0,epsrel=1e-9,limit=140)
    return g_bath*T**3/(2*math.pi**2)*result
