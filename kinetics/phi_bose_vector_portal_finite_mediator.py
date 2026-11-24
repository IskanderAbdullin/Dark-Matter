"""Finite-mass vector-portal CANDIDATE consistent with real-phi Yukawa charges.

Renormalizable off-diagonal real scalar Yukawa y phi bar(chi2)chi1+h.c.
under vectorlike U(1)_D requires q(chi1)=q(chi2) when phi is REAL
and gauge neutral. A possible additional Z2 makes chi1,phi odd,
chi2 even. Below electroweak breaking, an illustrative heavy kinetically
mixed A' couples the electromagnetic e current with epsilon*e and BOTH
Dirac chi vector currents with gD, assuming equal unit charges.

Tree physical t-channel e+chi_i -> e+chi_i with finite m_A' and
massless e, unpolarized spins. This is a CANDIDATE portal extension
(not unique and not proven phenomenologically allowed).

The MB <sigma_transport v_Moller> per-particle collision-frequency proxy
is NOT the actual linearized relativistic energy or momentum relaxation,
which requires full FD final-state Pauli, SM e/nu/gamma, in-medium
screening/thermal masses, an actual mediator UV completion and FRW.
No physical scalar DM Omega/BEC/CLASS result follows.
"""
import json
import math
from scipy.integrate import quad
from scipy.special import kve

from kinetics.microphysics_portal import (
    vector_contact_m2,vector_contact_dsdt,
    sigma_transport_massless_e,radiation_H,
)
from kinetics.gauge_charge_yukawa_audit import analyze as gauge_charge_analyze

ALPHA_EM=1/137.035999084
MCHI1=0.1
MCHI2=0.02
GE_EPLUS_EMINUS=4  # e-/e+ spin 2 each
GCHI_ONE_DIRAC_CHARGE=2  # one chi and its antiparticle separately


def vector_portal_setup(gD=0.1,epsilon=1e-5,mV=10.,alpha=ALPHA_EM):
    """Explicit electroweak-broken U(1)_D-EM mixing EFT benchmark only."""
    if not (mV>0 and alpha>0 and gD>=0 and epsilon>=0):
        raise ValueError("positive mediator mass, alpha and nonnegative couplings required")
    charge=gauge_charge_analyze(1,1,0,real_phi=True)
    assert charge["renormalizable_yukawa_gauge_invariant"]
    assert charge["chi2_diagonal_vector_current_nonzero"]
    electron=math.sqrt(4*math.pi*alpha)
    pair=gD*epsilon*electron
    return dict(gD=gD,epsilon=epsilon,mV_GeV=mV,
                electron_charge=electron,vector_coupling_product=pair,
                G_contact_GeV_minus2=pair/(mV*mV),
                q_chi1=1,q_chi2=1,
                dark_Z2="chi1,real_phi odd; chi2, SM and vector mediator even",
                scope="candidate vector portal, not full experimentally validated SM model")


def t_bounds_massless_e(s,mchi):
    if not (mchi>0 and s>mchi*mchi):
        raise ValueError("Physical chi e elastic s>mchi² required")
    tmin=-(s-mchi*mchi)**2/s
    return tmin,0.


def amplitude_squared(s,t,mchi,gD=0.1,epsilon=1e-5,mV=10.):
    """Unpolarized Dirac-vector 4-spin averaged |M|² with full spacelike A'."""
    a=t_bounds_massless_e(s,mchi)
    if t<a[0]*(1+1e-12) or t>1e-12:
        raise ValueError("t outside physical elastic range")
    Gt=vector_portal_setup(gD,epsilon,mV)["vector_coupling_product"]/(mV*mV-t)
    return vector_contact_m2(s,t,mchi,Gt)


def differential_ds_dt(s,t,mchi,gD=0.1,epsilon=1e-5,mV=10.):
    """Positive unpolarized tree dσ/dt (GeV^-4) for massless electron."""
    Gt=vector_portal_setup(gD,epsilon,mV)["vector_coupling_product"]/(mV*mV-t)
    return vector_contact_dsdt(s,t,mchi,Gt)


def finite_vector_cross_section(s,mchi,gD=0.1,epsilon=1e-5,mV=10.,
                                transport=True):
    """Actual total or COM momentum-transfer-weighted vacuum σ, GeV^-2."""
    if s<=mchi*mchi:return 0.
    tmin,_=t_bounds_massless_e(s,mchi)
    p2=(s-mchi*mchi)**2/(4*s)
    def integrand(x):
        t=tmin*x
        factor=-t/(2*p2) if transport else 1.
        return (-tmin)*differential_ds_dt(s,t,mchi,gD,epsilon,mV)*factor
    return quad(integrand,0,1,epsabs=0,epsrel=3e-10,limit=150)[0]


def MB_transport_collision_proxy(T,mchi=MCHI2,mV=10.,gD=.1,
                                 epsilon=1e-5,gstar=61.75,gchi=2,ge=4):
    """Transport-weighted 2-particle MB reaction density / chi density.

    For unequal MASSLESS-e and massive-chi, MB invariant reaction density:
        γT = gχ ge T/(32π⁴) ∫_{mχ²}∞ ds σT(s)
                                      (s−mχ²)² /sqrt(s) K1(sqrt(s)/T).
    This is NOT drag gamma_el (spectrum relaxation eigenvalue).
    """
    if not (T>0 and mchi>0 and mV>0 and gstar>0 and gchi>0 and ge>0):
        raise ValueError("physical positive temperature/mass/degrees required")
    x=mchi/T
    if x>=80:
        raise ValueError("Extremely dilute MB parents require separate log-tail analysis")
    def integrand(u):
        s=(T*u)**2
        return (finite_vector_cross_section(s,mchi,gD,epsilon,mV)*
                (s-mchi*mchi)**2*kve(1,u)*math.exp(-u))
    moment=quad(integrand,x,max(50.,x+50.),epsabs=0,
                epsrel=2e-6,limit=240)[0]
    gamma=gchi*ge*T*T/(16*math.pi**4)*moment
    nchi=gchi*mchi*mchi*T*kve(2,x)*math.exp(-x)/(2*math.pi**2)
    freq=gamma/nchi
    H=radiation_H(T,gstar)
    ref=vector_portal_setup(gD,epsilon,mV)
    return dict(T_GeV=T,mchi_GeV=mchi,gstar=gstar,
                physical_vector_mediator_GeV=mV,gD=gD,epsilon=epsilon,
                G_contact_GeV_minus2=ref["G_contact_GeV_minus2"],
                MB_chi_density_GeV3=nchi,
                MB_transport_weighted_forward_scattering_density_GeV4=gamma,
                MB_transport_weighted_per_chi_collision_frequency_GeV=freq,
                radiation_H_GeV=H,
                MB_transport_proxy_over_H_NOT_drag=freq/H,
                does_not_calibrate_BGK=True,
                dark_and_SM_heat_or_chemistry_not_solved=True)


def ultrarelativistic_MB_contact_limit(T,Gcontact,ge=4):
    """Analytical mchi/T->0, mV/T->infty MB *transport collision rate*.

    In massless limit σT(s)=7 G² s/(48π).
    MB direct moment <Eχ Ee (1-cosθ)^2>=9T²*(4/3).
    Thus <σT v>=7 G² T²/(2π) and n_e_MB=ge T³/π².
    """
    if not (T>0 and Gcontact>=0 and ge>0):
        raise ValueError("Positive physical inputs")
    return 7*ge*Gcontact*Gcontact*T**5/(2*math.pi**3)


def audit():
    setups=[]
    for mchi in (MCHI1,MCHI2):
        for T in (.1,.5,1.,2.):
            setups.append(MB_transport_collision_proxy(T,mchi))
    cross=[]
    for sqrt_s in (.1,1.,5.,10.):
        s=sqrt_s*sqrt_s
        for mchi in (MCHI1,MCHI2):
            if s<=mchi*mchi:continue
            p=vector_portal_setup()
            ct=sigma_transport_massless_e(s,mchi,p["G_contact_GeV_minus2"])
            finite=finite_vector_cross_section(s,mchi,transport=True)
            cross.append(dict(mchi=mchi,sqrt_s=sqrt_s,
                              sigma_T_finite_GeV_minus2=finite,
                              sigma_T_contact_GeV_minus2=ct,
                              ratio_finite_to_contact=finite/ct))
    return dict(schema_version=1,scope="CANDIDATE same-charge vectorlike chi1 chi2 U1D kinetically mixed A' electron-only broken-EW finite propagator, massless e and MB unblocked transport proxy only, not gamma_el, not full SM FRW or bosonic DM",
                gauge_charges=gauge_charge_analyze(1,1,0,real_phi=True),
                setup=vector_portal_setup(),cross_section_rows=cross,
                thermal_rows=setups,class_ready=False)


if __name__=="__main__":print(json.dumps(audit(),indent=2))
