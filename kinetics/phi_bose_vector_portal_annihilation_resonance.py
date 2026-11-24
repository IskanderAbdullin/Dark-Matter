"""ON-SHELL dark-vector resonance in χ+antiχ <-> e-+e+ chemical portal.

L ⊃ gD A'μ [barχ1 γμ χ1 +barχ2 γμ χ2]
        + ε e A'μ bar e γμ e
        + y φ barχ2 χ1 +h.c.
Real φ and χ1/χ2 have compatible qχ1=qχ2 under dark U(1), Z2 as specified.

Tree spin-averaged χχbar-> e-e+ s-channel finite-width vector:
   σ(s)= (gD ε e)^2/(12π) * (s+2mχ²)/
          [βχ ((s-mV²)²+mV² ΓV²)],  βχ=sqrt(1-4mχ²/s).

ΓV is a declared LOWEST-ORDER MINIMAL χ1+χ2+e partial-width sum.
Kinetic mixing to muons/quarks etc also contributes physically,
but is intentionally excluded from this partial-channel comparison.

Thermal MB forward reaction density is integrated with the exact
Lorentz-invariant flux and a TAN transformation across the resonance:
s=mV² + mV ΓV tanθ, ds/[(s-mV²)²+mV² ΓV²]=dθ/(mV ΓV).

The exact MB reaction-density 'energy-weighted' quantity follows from
-∂γ/∂(1/T): derivative of [T K1(√s/T)] =
 [2T² K1(√s/T)+T√s K0(√s/T)].

Direct e-e+->χχbar has σ_reverse=βχ² σ_forward by detailed balance
and exactly same γ(T) after including physical initial degeneracies,
threshold and flux. Thus net rate at fixed separate thermal MB baths,
chemical µ/T and total Qdark=-QSM can be formed consistently.
NOT quantum FD/BE Pauli, full actual SM, finite thermal width,
chemical-relaxation eigenvalue, FRW, 86eV BEC or validated relic.
"""
import math
import json
import numpy as np
from scipy.integrate import quad
from scipy.special import kve

from kinetics.phi_bose_vector_portal_finite_mediator import (
    ALPHA_EM,MCHI1,MCHI2,vector_portal_setup
)
from kinetics.phi_bose_yukawa_hubble_gate import hubble_rd


def vector_to_Dirac_width(mV,mf,g):
    if min(mV,mf)<=0 or g<0:raise ValueError("Physical positive masses and coupling")
    if mV<=2*mf:return 0.
    z=(mf/mV)**2
    return g*g*mV/(12*math.pi)*(1+2*z)*math.sqrt(1-4*z)


def minimal_vector_width(mV=10.,gD=.1,epsilon=1e-5):
    """Tree inclusive chi1/chi2 and MASSLESS e widths, not full SM width."""
    v=vector_portal_setup(gD,epsilon,mV)
    channels={
        "chi1+antichi1":vector_to_Dirac_width(mV,MCHI1,gD),
        "chi2+antichi2":vector_to_Dirac_width(mV,MCHI2,gD),
        "e-+e+":gD*0.+(epsilon*v["electron_charge"])**2*mV/(12*math.pi),
    }
    return dict(width_GeV=sum(channels.values()),channels_GeV=channels,
                missing_open_muon_hadronic_photon_Z_mixing=True)


def annihilation_sigma(s,mchi,mV=10.,gD=.1,epsilon=1e-5,GammaV=None):
    """χχbar->e+e- total unpolarized σ [GeV^-2], width in resonance."""
    if not (s>4*mchi*mchi):return 0.
    if GammaV is None:GammaV=minimal_vector_width(mV,gD,epsilon)["width_GeV"]
    if GammaV<=0:raise ValueError("Finite vector width required")
    z=4*mchi*mchi/s
    coupling=vector_portal_setup(gD,epsilon,mV)["vector_coupling_product"]
    D=(s-mV*mV)**2+(mV*GammaV)**2
    return coupling**2/(12*math.pi)*(s+2*mchi*mchi)/(math.sqrt(1-z)*D)


def inverse_ee_sigma(s,mchi,mV=10.,gD=.1,epsilon=1e-5,GammaV=None):
    """Inverse e+e-->χχbar with correct physical 2->2 relative p² flux."""
    if s<=4*mchi*mchi:return 0.
    return (1-4*mchi*mchi/s)*annihilation_sigma(
        s,mchi,mV,gD,epsilon,GammaV)


def s_channel_thermodynamic_MB(T,mchi=.02,mV=10.,gD=.1,epsilon=1e-5,
                               GammaV=None,gstar=61.75,
                               max_thermal_energy=50.):
    """MB forward χχbar reaction density γ, energy transfer Eγ and Γ/H.

    Completely independent of a constant 'chemical equilibration rate'.
    The resonance is mapped exactly to theta and separately partitioned
    at thermal momentum scales, NOT integrated by naive uniform √s Gauss.

    The constant-width Breit-Wigner approximation, Bessel-weighted ideal
    MB plasma and frozen mediator mass are explicitly restricted model
    assumptions. γ is FORWARD only; actual net requires reverse.
    """
    if min(T,mchi,mV,gstar)<=0 or max_thermal_energy<=8:
        raise ValueError("Physical parameters positive, tail >8T")
    if GammaV is None:
        GammaV=minimal_vector_width(mV,gD,epsilon)["width_GeV"]
    if not (GammaV>0 and mV>2*mchi):
        raise ValueError("Open vector resonance mV>2mchi and ΓV>0 required")
    gprod=vector_portal_setup(gD,epsilon,mV)["vector_coupling_product"]
    w=mV*GammaV
    smin=4*mchi*mchi
    Emax=max(max_thermal_energy*T,mV+30*T)
    smax=Emax*Emax
    def theta(E):
        return math.atan((E*E-mV*mV)/w)
    lo=theta(2*mchi)
    hi=theta(Emax)
    # Explicitly resolve BOTH the potentially narrow thermally populated
    # non-resonant left θ shoulder and the mV-on-shell pole.
    Ebreaks=[2*mchi,mV,Emax]
    Ebreaks += [u*T for u in
                 (.2,.5,1.,1.5,2.,3.,4.,5.,6.,8.,10.,15.,20.,30.,40.,50.)]
    Ebreaks += [mV+g*GammaV for g in
                 (-100.,-30.,-10.,-3.,-1.,1.,3.,10.,30.,100.)]
    breaks=sorted(set([lo,hi]+[theta(E) for E in Ebreaks if 2*mchi<E<Emax]))
    if min(np.diff(breaks))<=0:
        raise ArithmeticError("Breit-Wigner integration physical segment collapsed")
    def integral_pair(the):
        s=mV*mV+w*math.tan(the)
        if s<=smin:return 0.,0.
        E=math.sqrt(s)
        beta=math.sqrt(1-smin/s)
        common=(gprod*gprod/(12*math.pi))*(s+2*mchi*mchi)/beta
        common*= (s-smin)*E/w
        x=E/T
        k1=kve(1,x)*math.exp(-x)
        k0=kve(0,x)*math.exp(-x)
        reaction=common*T*k1
        # Differentiation of MB Lorentz invariant density wrt beta=1/T
        heat=common*(2*T*T*k1+T*E*k0)
        return reaction,heat
    integrals=[]
    for comp in (0,1):
        integrals.append(sum(quad(
            lambda th:integral_pair(th)[comp],a,b,epsabs=0.,
            epsrel=3e-9,limit=160)[0]
            for a,b in zip(breaks[:-1],breaks[1:])))
    pref=4/(32*math.pi**4)
    gamma=pref*integrals[0]
    qloss=pref*integrals[1]
    x=mchi/T
    n=2*mchi*mchi*T*kve(2,x)*math.exp(-x)/(2*math.pi**2)
    H=hubble_rd(T,gstar)
    return dict(T_GeV=T,mchi_GeV=mchi,mV_GeV=mV,
                GammaV_minimal_channels_GeV=GammaV,
                gamma_forward_MB_GeV4=gamma,
                energy_forward_MB_GeV5=qloss,
                mean_lab_total_energy_per_event_GeV=qloss/gamma,
                n_chi_MB_GeV3=n,
                forward_per_chi_rate_GeV=gamma/n,
                H_RD_GeV=H,
                forward_over_H_NOT_NET=gamma/n/H,
                num_theta_segments=len(breaks)-1,
                resonant_s_mapping="s=mV²+mVΓ tan θ; exact ds/BW=dθ/(mVΓ)",
                no_thermal_loop_or_full_SM_or_relic=True)



def narrow_width_vector_onshell_MB(T,mchi=.02,mV=10.,gD=.1,epsilon=1e-5,
                                   GammaV=None):
    """INDEPENDENT narrow-width sequential χχbar↔A'→e+e− reference.

    Exact tree two-body partial widths and a spin-3 massive real vector
    MB equilibrium number nA', with Lorentz-averaged decay K1/K2:
      γ_NWA=nA' Γ(A'→χχbar) Br(A'→e+e−) K1/K2.
    This is a physical real-intermediate-state (RIS) contribution. If A'
    is later explicitly evolved as a separate kinetic population, one MUST
    subtract RIS from a simultaneous full 2->2 Breit-Wigner kernel to
    avoid double counting. NWA is invalid when off-shell dominates at T≪mV.
    """
    if min(T,mchi,mV)<=0 or mV<=2*mchi:
        raise ValueError("Open vector resonance, positive T required")
    info=minimal_vector_width(mV,gD,epsilon)
    if GammaV is None:GammaV=info["width_GeV"]
    if GammaV<=0:raise ValueError("positive total vector width required")
    partial_chi=vector_to_Dirac_width(mV,mchi,gD)
    partial_e=info["channels_GeV"]["e-+e+"]
    z=mV/T
    nV=3*mV*mV*T*kve(2,z)*math.exp(-z)/(2*math.pi**2)
    avg_dilation=kve(1,z)/kve(2,z)
    gamma=nV*partial_chi*(partial_e/GammaV)*avg_dilation
    energy_per_event=2*T+mV*kve(0,z)/kve(1,z)
    return dict(T_GeV=T,mchi_GeV=mchi,mV_GeV=mV,
                effective_total_GammaV_GeV=GammaV,
                GammaV_to_chichi_GeV=partial_chi,
                GammaV_to_electron_GeV=partial_e,
                Br_electron=partial_e/GammaV,
                spin_three_vector_MB_density_GeV3=nV,
                decay_time_dilation_K1_over_K2=avg_dilation,
                real_intermediate_onshell_forward_GeV4=gamma,
                real_intermediate_energy_forward_GeV5=gamma*energy_per_event,
                mean_lab_energy_per_event_GeV=energy_per_event,
                RIS_double_count_if_explicit_vector_species=True,
                NWA_valid_only_if_on_shell_peak_dominates=True)

def net_MB_two_temperature_ledger(Tchi,Te,mchi=.02,
                                   mu_chi_sum=0.,mu_e_sum=0.,**kwargs):
    """Full forward-minus-reverse MB reaction and exact quantum-free entropy.

    Chemical potentials are sums of particle+antiparticle µ for each pair,
    not a single arbitrary BGK number reservoir. True U1 and electron
    charge each conserved: χ and antiχ lost/created in pairs.
    """
    if Tchi<=0 or Te<=0:raise ValueError("Both bath temperatures positive")
    d=s_channel_thermodynamic_MB(Tchi,mchi,**kwargs)
    e=s_channel_thermodynamic_MB(Te,mchi,**kwargs)
    fchi=math.exp(mu_chi_sum/Tchi)
    fe=math.exp(mu_e_sum/Te)
    R=d["gamma_forward_MB_GeV4"]*fchi-e["gamma_forward_MB_GeV4"]*fe
    Q=d["energy_forward_MB_GeV5"]*fchi-e["energy_forward_MB_GeV5"]*fe
    H=Q*(1/Te-1/Tchi)+R*(mu_chi_sum/Tchi-mu_e_sum/Te)
    return dict(mchi_GeV=mchi,Tchi_GeV=Tchi,Te_GeV=Te,
                R_forward_minus_reverse_GeV4=R,
                Q_dark_energy_loss_GeV5=Q,
                Q_dark_GeV5=-Q,Q_e_GeV5=Q,
                U1_net_defect=0.,electric_charge_net_defect=0.,
                H_ideal_MB_GeV4=H,
                microscopic_H_nonnegative_within_numeric_precision=H>=-1e-25,
                forward_model=d,inverse_detailed_balance_model=e,
                detailed_balance_one_event="F/R=exp(E_total*(1/Te-1/Tchi)+Σmu_chi/Tchi-Σmu_e/Te)",
                limits="ideal MB fixed separate baths, no Pauli blocking or BE factors, no in-medium vector width/SM/FRW; not true chemical relaxation or physical DM relic")


def audit():
    mass=MCHI2; v=[]
    for T,gstar in ((.1,17.25),(.5,61.75),(1.,61.75),(2.,75.)):
        v.append(s_channel_thermodynamic_MB(T,mass,gstar=gstar))
    ledger=net_MB_two_temperature_ledger(.8,1.,mass)
    return dict(schema_version=1,
                scope="conditional finite-width dark vector χ+antiχ↔e+e− MB microreversible chemical/energy source, not full quantum SM or FRW/Ω/BEC",
                partial_width=minimal_vector_width(),
                chi2_thermal_rows=v,nonequilibrium_two_T=ledger,
                class_ready=False)


if __name__=="__main__":
    print(json.dumps(audit(),indent=2))
