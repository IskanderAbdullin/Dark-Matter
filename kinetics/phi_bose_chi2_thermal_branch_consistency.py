"""Thermal-branch consistency gate for minimal stable χ2 thermal relic.

The prior symmetric 20-MeV χ2+antiχ2 fixed-g*=10.75 MB ODE
implicitly PRESCRIBES a common thermal bath χ2/φ/SM down to x=m/T~10.
But the explicit gauge-compatible 10GeV vector A' portal to e± has
a very small physically normalized transport-weighted χ2/e collision
frequency by then.

We calculate two INDEPENDENT instantaneous forward thermal rates
on the same ideal MB background:
   Γpair/H = n_antiχ2^(eq) <σv>_exact_(χ2anti2→φφ) / H,
   νT,e/H = [γ_transport,MB(χ2 e→χ2 e)/nχ2]/H,
where the first exact χ2-virtual-Yukawa σ(s) has full Dirac t/u
and the second exact finite-mass vector t-channel has its own Bessel
integration. Both were independently checked in their source tests.

νT,e is NOT an eigenvalue of the nonlinear Boltzmann operator or
relativistic energy-relaxation gamma. A tiny νT/H is a strong
necessary-equilibration warning, but not a proof of failed dark
χ2-φ internal thermalization. Electron mass, annihilation, muons,
hadrons, full FD blocking, FRW thermal spectra, dark bath Tphi,
in-medium widths and SM entropy evolution are omitted.

Thus we may not use the illustrative thermal ODE overclosure as a
literal prediction of the declared A' benchmark, only as a scoped
negative thermal branch.
"""
import json
import math
from scipy.special import kve

from kinetics.phi_bose_stable_chi2_relic_gate import (
    MCHI2,Y_REF,exact_MB_pair_sigma_v,exact_scalar_pwave_b,
    vector_ee_swave_a,MPL_NONREDUCED,
)
from kinetics.phi_bose_vector_portal_finite_mediator import (
    MB_transport_collision_proxy,
)


def instantaneous_thermal_consistency(x,gstar=10.75):
    if not (3<=x<=20 and gstar>0):
        raise ValueError("Physical x in tested thermal window 3..20")
    T=MCHI2/x
    # Stable exp(-x) times scaled K2 is safe for x<=20.
    n_eq=2*MCHI2*MCHI2*T*kve(2,x)*math.exp(-x)/(2*math.pi**2)
    H=1.66*math.sqrt(gstar)*T*T/MPL_NONREDUCED
    sv=exact_MB_pair_sigma_v(x)
    v=MB_transport_collision_proxy(
        T,mchi=MCHI2,mV=10.,gD=.1,epsilon=1e-5,
        gstar=gstar,gchi=2,ge=4)
    p_wave=6*exact_scalar_pwave_b(Y_REF)/x
    a=vector_ee_swave_a()
    return dict(x=x,T_GeV=T,
                n_chi2_one_charge_equilibrium_GeV3=n_eq,
                H_GeV=H,
                chi2_to_phi_phi_full_GG_sigma_v_GeV_minus2=sv,
                exact_GG_to_pwave_6b_over_x=sv/p_wave,
                chi2_pair_forward_per_chi_over_H=n_eq*sv/H,
                chi2_pair_offshell_vector_e_per_chi_over_H=n_eq*a/H,
                finite_mediator_e_transport_proxy_over_H_NOT_drag=
                    v["MB_transport_proxy_over_H_NOT_drag"],
                ratio_elastic_proxy_to_dark_pair_forward=
                    v["MB_transport_proxy_over_H_NOT_drag"]/(n_eq*sv/H),
                conditional_consistency_warning="e± 10GeV A' transport proxy far below H near χ2 chemical freezeout; assumed common SM/χ/φ thermal history is not derived",
                plasma_electron_assumed_massless=True,
                dark_phi_thermal_contact_not_proven=True,
                no_physical_gamma_el_or_full_SM_FWR=True,
                no_actual_relativistic_relic_prediction=True)



def conditional_electron_transport_crossing(gstar_rho=61.75):
    """Solve ν_T,e(T)/H(T)=1 at constant supplied external g*rho.

    A positive transport frequency >H is NOT sufficient to assert
    physical relativistic kinetic equilibrium (relaxation eigenmode
    can be much slower). A crossing far ABOVE χ2 chemical freezeout
    signals the particular shared-SM thermal toy assumption needs
    actual full nonthermal and dark-bath treatment.
    """
    from scipy.optimize import brentq
    def ratio(T):
        return MB_transport_collision_proxy(
            T,mchi=MCHI2,mV=10.,gD=.1,epsilon=1e-5,
            gstar=gstar_rho,gchi=2,ge=4)[
            "MB_transport_proxy_over_H_NOT_drag"]
    Tcross=float(brentq(lambda T:ratio(T)-1.,.1,2.,xtol=1e-8))
    return dict(fixed_gstar_rho=gstar_rho,
                conditional_transport_proxy_crossing_T_GeV=Tcross,
                ratio_below_cross_at_T_0p1GeV=ratio(.1),
                ratio_above_cross_at_T_2GeV=ratio(2.),
                mchi2_over_crossing_T=MCHI2/Tcross,
                no_actual_physical_kinetic_decoupling_temperature=True,
                BBN_and_neutrino_electron_mass_not_included=True)


def epsilon_needed_for_illustrative_electron_transport_proxy_equal_H(
        x=8.,epsilon_reference=1e-5):
    """Mathematical rescaling εcrit=ε0/sqrt(ν_T/H) at target x.

    The original mV10/gD.1 vector source is perturbative ∝ε² in
    this prescribed vacuum/electron-only proxy. This is NOT a
    phenomenologically permitted ε, excludes no real plasma models.
    """
    if epsilon_reference<=0:raise ValueError("Positive ε reference")
    r=instantaneous_thermal_consistency(x)[
        "finite_mediator_e_transport_proxy_over_H_NOT_drag"]
    return dict(x=x,T_GeV=MCHI2/x,
                original_epsilon=epsilon_reference,
                original_transport_proxy_over_H=r,
                epsilon_for_proxy_over_H_1=epsilon_reference/math.sqrt(r),
                experimentally_allowed_unknown=True,
                true_heat_relaxation_not_computed=True)

def audit():
    rows=[instantaneous_thermal_consistency(x) for x in (3,5,8,10,12,15,20)]
    return dict(schema_version=1,scope="strictly conditional shared-MB bath and massless e vector transport proxy vs independent quantum-Yukawa exact GG pair forward at x3..20, not full selfconsistent cosmology",
                rows=rows,class_ready=False)


if __name__=="__main__":
    print(json.dumps(audit(),indent=2))
