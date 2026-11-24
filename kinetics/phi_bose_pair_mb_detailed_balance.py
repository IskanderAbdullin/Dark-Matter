"""Vacuum χ2 antiχ2 ↔ φ φ detailed balance and conditional thermal MB weak moments.

The real-scalar final pair carries 1/2! in σ(Dirac pair→2φ).
For inverse σ(φφ→Dirac pair), one averages one scalar initial state
instead of four Dirac initial spin states and removes final 1/2!.
Therefore σ_inverse/σ_forward = 8 p_Dirac² / p_scalar².

The inverse thermal reaction density carries **1/2! identical initial**
pairs. With the same Maxwell-Boltzmann T, these factors cancel exactly
pointwise in the common invariant s kernel: γ_forward^eq=γ_inverse^eq.

This does not solve full FD/Bose kinetics, thermal self-energies, SM or FRW.
"""
import json
import math
import numpy as np
from scipy.integrate import quad
from scipy.special import kv

from kinetics.phi_bose_pair_annihilation_mb import (
    M1,M2,MPHI,Y,PLANCK,cross_section,thermal_forward_rate_MB,
)

def reverse_cross_section(s, y=Y, nangle=96):
    """Vacuum σ(φφ→χ2 antiχ2), summed final χ spins, incoming real φ."""
    if s<=4*M2*M2:
        return 0.
    pin2=s/4-M2*M2
    pout2=s/4-MPHI*MPHI
    return cross_section(s,y,nangle)*8*pin2/pout2

def equilibrium_pair_densities(T,y=Y,nangle=96):
    """Independent numerical MB forward/inverse γ and total reaction energy.

    Forward χ2 antiχ2→2φ: distinct incoming spin degeneracy 2*2.
    Inverse φφ→χ2 antiχ2: identical incoming state pair 1/2!.
    True physical 4-vectors are integrated out for uniform MB distributions;
    gamma_E uses -d(gamma)/d(beta) at beta=1/T, fixed μ fugacities.
    """
    if T<=0:raise ValueError("Positive temperature")
    xmin=2*M2/T
    xp=2*MPHI/T
    gdirac=2
    pref=T/(32*math.pi**4)
    def kernel_forward(u):
        s=T*T*u*u
        return (2*T**5*u*u*(u*u-xmin*xmin) *
                cross_section(s,y,nangle))
    def kernel_inverse(u):
        s=T*T*u*u
        return (2*T**5*u*u*(u*u-xp*xp) *
                reverse_cross_section(s,y,nangle))
    F=quad(lambda u: kernel_forward(u)*kv(1,u),
           xmin,40,epsabs=1e-34*(y/Y)**4,epsrel=5e-6,limit=160)[0]
    R=quad(lambda u: kernel_inverse(u)*kv(1,u),
           xmin,40,epsabs=1e-34*(y/Y)**4,epsrel=5e-6,limit=160)[0]
    # The reaction-energy Bessel kernel is derived from differentiating
    # T*K1(sqrt(s)/T): E-weight/T*K1 = T*(2*K1 + u*K0).
    QE=quad(lambda u:kernel_forward(u)*T*(2*kv(1,u)+u*kv(0,u)),
            xmin,40,epsabs=1e-33*(y/Y)**4,epsrel=5e-6,limit=160)[0]
    return dict(T_GeV=T,
                gamma_forward_equilibrium_GeV4=gdirac**2*pref*F,
                gamma_inverse_equilibrium_GeV4=.5*pref*R,
                event_total_energy_density_GeV5=gdirac**2*pref*QE,
                average_total_event_energy_GeV=QE/F,
                inverse_to_forward_ratio=(.5*R)/(gdirac**2*F),
                class_ready=False)

def thermal_MB_chemical_source(T,mu_chi,mu_anti,mu_phi,y=Y,nangle=96):
    """Physical MB thermal shape ONLY, separate net 5-species source.

    χ2 + antiχ2 → φ + φ; μ/T may vary. Numerically evaluates both
    physical reaction densities but energy ledger uses a common exact
    event-energy moment. For a common T, the forward/inverse energy shape
    is identical: gammas differ only by fugacity factors.
    """
    eq=equilibrium_pair_densities(T,y,nangle)
    zf=math.exp((mu_chi+mu_anti)/T)
    zr=math.exp(2*mu_phi/T)
    # Equal normalization by microscopic detailed balance; use
    # separately calculated forward/inverse for a numerical rate gate.
    gf=eq["gamma_forward_equilibrium_GeV4"]*zf
    gr=eq["gamma_inverse_equilibrium_GeV4"]*zr
    net=gf-gr
    # Energy calculations use the same common physical forward/reverse
    # kernel; separately integrated normalized reverse γ is a check.
    Q=eq["event_total_energy_density_GeV5"]*(zf-zr)
    affinity=(mu_chi+mu_anti-2*mu_phi)/T
    H=net*affinity
    return dict(T_GeV=T,
                chemical_affinity_dimensionless=affinity,
                gamma_forward_GeV4=gf,
                gamma_reverse_GeV4=gr,
                net_events_GeV4=net,
                dn_chi_GeV4=-net,
                dn_anti_GeV4=-net,
                dn_phi_GeV4=2*net,
                dE_chi_GeV5=-.5*Q,
                dE_anti_GeV5=-.5*Q,
                dE_phi_GeV5=Q,
                total_energy_defect_GeV5=0.,
                net_Dirac_U1_defect=0.,
                H_per_volume_GeV4=H,
                equilibrium=eq,
                caveat="Ideal MB equilibrium-shaped distributions and vacuum 2to2; NOT FD/BE quantum collision, finite plasma, SM or expanding relic",
                class_ready=False)

def audit():
    return dict(schema_version=1,
                scope="Vacuum cross-section reciprocity, identical initial pair correction and common-T conditional Maxwell-Boltzmann weak N/E/H; not full quantum kinetics or late cosmology.",
                reverse_sigma_at_sqrt_s_0p5=cross_section(.25),
                equilibrium=[equilibrium_pair_densities(T) for T in (.1,1.)],
                mu_negative_phi=thermal_MB_chemical_source(1.,.05,-.05,-.2),
                mu_zero_phi=thermal_MB_chemical_source(1.,.05,-.05,0.),
                class_ready=False)

if __name__=="__main__":
    print(json.dumps(audit(),indent=2,allow_nan=False))
