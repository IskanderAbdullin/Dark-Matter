"""Mandatory visible channels of a kinetically mixed 10-GeV dark vector.

The candidate A' couples to electromagnetic current ε e J_EM and BOTH
χ1/χ2 vector currents; a single e+e− line is not the full visible width.

For charged lepton ℓ:
 Γ(A'→ℓℓ) = α ε² mV/3 (1+2mℓ²/mV²) sqrt(1−4mℓ²/mV²).
For hadrons, physical inclusive data must supply R_had(mV²) =
σ(e+e−→hadrons)/σ(e+e−→µ+µ−) for a data-defined convention.
Then Γhad=Γmu,pointlike(mV) * R_had at the chosen mass, only when the
R definition and muon reference are consistent. Do NOT invent R from
free b quarks near BB threshold or from low-energy quark masses.
This function defaults to R_had=None and reports unresolved hadrons,
not a fake complete SM width.

Near an isolated NARROW A' pole, factorization through an explicit
massive spin-1 intermediate state yields
 γ(χχbar→all charged visible) ≈ nA Γχ Br_visible K1/K2.
The ratio to the previously computed e-only NWA is
 (Γ_visible/Γ_e) * (Γ_min/Γ_total).
The old e-only Γ_min already includes both DARK Dirac channels.
This ratio is a conditional tree/NWA result, not a full thermal quantum
collision source, and no actual current experimental R or ΓV(T) is used.

A gauged 10-GeV kinetically mixed photon can have thresholds, Υ
resonances, Z mixing, running α and QCD corrections; no claim of a
validated experimentally allowed benchmark, BEC or 86-eV relic.
"""
import json
import math
from kinetics.phi_bose_vector_portal_finite_mediator import (
    ALPHA_EM,vector_portal_setup,
)
from kinetics.phi_bose_vector_portal_annihilation_resonance import (
    minimal_vector_width,narrow_width_vector_onshell_MB,
)

ME_GEV=.00051099895
MMU_GEV=.1056583755
MTAU_GEV=1.77686


def charged_lepton_width(mV,epsilon,mass):
    if min(mV,mass)<=0 or epsilon<0:raise ValueError("Positive physical masses")
    if mV<=2*mass:return 0.
    z=(mass/mV)**2
    return ALPHA_EM*epsilon*epsilon*mV/3*(1+2*z)*math.sqrt(1-4*z)


def visible_width_gate(mV=10.,gD=.1,epsilon=1e-5,R_had=None):
    if mV<=0 or gD<0 or epsilon<0:raise ValueError("Valid mediator")
    if R_had is not None and (not math.isfinite(R_had) or R_had<0):
        raise ValueError("Rhad must be a physically supplied nonnegative ratio")
    widths={
        "electron":charged_lepton_width(mV,epsilon,ME_GEV),
        "muon":charged_lepton_width(mV,epsilon,MMU_GEV),
        "tau":charged_lepton_width(mV,epsilon,MTAU_GEV),
    }
    # Definition requires that physical R is quoted against the massless
    # pointlike muon σ, not a reconstructed free-quark pole threshold.
    ref_massless_mu=ALPHA_EM*epsilon*epsilon*mV/3.
    hadronic=None if R_had is None else ref_massless_mu*R_had
    lepton_sum=sum(widths.values())
    if hadronic is not None:
        widths["hadronic_data_driven_R"]=hadronic
    min_dark= minimal_vector_width(mV,gD,epsilon)
    total_partial=sum(widths.values())+sum(
        x for k,x in min_dark["channels_GeV"].items() if k!="e-+e+")
    baseline=min_dark["width_GeV"]
    fraction=lepton_sum/widths["electron"] if widths["electron"] else None
    if hadronic is not None and widths["electron"]:
        f_all=(lepton_sum+hadronic)/widths["electron"]
    else:
        f_all=None
    return dict(mV_GeV=mV,gD=gD,epsilon=epsilon,
                gamma_charged_leptons_GeV=widths,
                gamma_leptons_only_GeV=lepton_sum,
                gamma_hadrons_data_driven_GeV=hadronic,
                R_had_user_supplied=R_had,
                dark_plus_visible_partial_GammaV_GeV=total_partial,
                dark_plus_eonly_baseline_GammaV_GeV=baseline,
                lepton_to_electron_width_ratio=fraction,
                all_visible_to_electron_width_ratio=f_all,
                total_width_relative_shift=total_partial/baseline-1,
                missing_observed_R_had_and_running_coupling=R_had is None,
                no_full_sm_electroweak_or_thermal_width=True)


def narrow_resonant_all_charged_MB(T=1.,mchi=.02,R_had=None,
                                  mV=10.,gD=.1,epsilon=1e-5):
    """NWA rescaling from independently derived real-vector three spin
    channels. A supplied R_had is a PARAMETER, not a measured value.
    """
    x=visible_width_gate(mV,gD,epsilon,R_had)
    baseline=narrow_width_vector_onshell_MB(T,mchi,mV,gD,epsilon)
    w=x["dark_plus_visible_partial_GammaV_GeV"]
    w0=x["dark_plus_eonly_baseline_GammaV_GeV"]
    gamma_e=baseline["real_intermediate_onshell_forward_GeV4"]
    correction=w0/w
    electron=gamma_e*correction
    leptons=gamma_e*x["lepton_to_electron_width_ratio"]*correction
    all_visible=(None if R_had is None else
                 gamma_e*x["all_visible_to_electron_width_ratio"]*correction)
    return dict(T_GeV=T,mchi_GeV=mchi,epsilon=epsilon,mV_GeV=mV,
                forward_ee_NWA_GeV4=electron,
                forward_all_lepton_pairs_NWA_GeV4=leptons,
                forward_lepton_plus_hadron_NWA_GeV4=all_visible,
                changed_total_width_correction=correction,
                lepton_only_relative_to_prior_eonly=leptons/gamma_e,
                visible_widths=x,
                physically_measured_R_had_used=False,
                no_true_SM_plasma_or_dark_relic=True)


def audit():
    return dict(schema_version=1,
                scope="mandatory e,mu,tau tree A' widths and user-supplied hadronic R, exact NWA ratio for visible channels, no observed R/final thermal plasma",
                default_without_hadron_data=visible_width_gate(),
                illustrative_only_Rhad_3p5=visible_width_gate(R_had=3.5),
                resonance_T1_leptons=narrow_resonant_all_charged_MB(),
                resonance_T1_hadronic_example=narrow_resonant_all_charged_MB(R_had=3.5),
                class_ready=False)


if __name__=="__main__":print(json.dumps(audit(),indent=2))
