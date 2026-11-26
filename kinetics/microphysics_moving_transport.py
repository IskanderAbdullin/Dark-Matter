"""Momentum-dependent UNBLOCKED transport collision frequency for moving chi.

Exact incoming Moller flux of massless e in an isotropic FD bath, with
s=M^2+2k(E_chi-p*cos theta), v_M=1-(p/E_chi)*cos(theta).
Average CM (1-cos theta_CM) sigma; NOT a physical kinetic relaxation
eigenvalue, NOT final-state Pauli-corrected, NOT an energy-transfer kernel.
No full SM neutrino bath, UV completion, thermal mass or BEST coupling.

Natural units M,p,T in GeV, G in GeV^-2, result in GeV.
"""
from __future__ import annotations
import math
from scipy.integrate import quad
from microphysics_portal import (sigma_transport_massless_e,
                                 bath_transport_frequency_chi_at_rest)


def moving_chi_transport_frequency(M,p,T,G,g_bath=4):
    vals=(M,p,T,G,g_bath)
    if (not all(math.isfinite(v) for v in vals) or
            M<=0 or p<0 or T<=0 or G<0 or g_bath<=0):
        raise ValueError("finite M,T,g>0 and p,G>=0 required")
    if G==0:return 0.
    if p==0:
        return bath_transport_frequency_chi_at_rest(M,T,G,g_bath)
    energy=math.hypot(p,M)
    beta=p/energy
    def angular_average(k):
        def at_mu(mu):
            flux=1-beta*mu
            s=M*M+2*k*energy*flux
            return flux*sigma_transport_massless_e(s,M,G)
        ans,_=quad(at_mu,-1.,1.,epsabs=0.,epsrel=1e-10,limit=100)
        return ans/2
    def at_x(x):
        if x<=0:
            return 0.
        f=math.exp(-x)/(1+math.exp(-x))
        return x*x*f*angular_average(T*x)
    value,_=quad(at_x,0,40,epsabs=0,epsrel=1e-9,limit=140)
    return g_bath*T**3/(2*math.pi**2)*value
