"""Exact asymmetric stable-fermion abundance floor in real-phi Yukawa U1.

Lagrangian (χ1,antiχ1,χ2,antiχ2,φ), with real neutral φ, allows
χ1↔χ2+φ and χ_i+antiχ_i↔2φ; a common vectorlike A' coupling gives
χ_i+e↔χ_i+e and χ_i+antiχ_i↔SM+antiSM. ALL conserve Q_D
  Q_D = Nχ1 − Nantiχ1 + Nχ2 − Nantiχ2.
Real φ, A' and SM have Q_D=0.

Even if all charge-symmetric fermion pairs perfectly annihilate, the
surviving total charged-fermion number obeys N_charged>=|Q_D|. Since
mχ1>=mχ2, the late physical energy density obeys
  ρ_dark,stable_fermions >= mχ2 * |Q_D|/a³.
In terms of today's *actual* SM entropy S_SM,0 (which need not have
been comovingly conserved in the early coupled plasma), define
  eta_D,0=Q_D,comoving/S_SM,0.
Then Ω_stablefermion h² >= mχ2(s_SM,0/(ρcrit/h²)) |eta_D,0|.
No thermal MB assumption, no Bose condensate/φ PSD, no Hubble model,
no chemical equilibrium or elastic Γ required: this is an algebraic
Noether-charge bound, under the stated field content and no U1 breaking.

It is NOT a lower bound on the neutral φ component. If q1=q2=0 and
additional charge-breaking terms are allowed, explicitly rederive model.
"""
from __future__ import annotations
import math
import json

MCHI1_GEV=.1
MCHI2_GEV=.02
PHI_GEV=86e-9
S0_CM3=2891.2
RHOC_OVER_H2_GEV_CM3=1.05375e-5
RELIC_OBS_OMEGA_C_H2=.120


def exactly_conserved_dark_charge(N):
    """N=(χ1,antiχ1,χ2,antiχ2,realφ), physical comoving counts."""
    if len(N)!=5 or not all(math.isfinite(float(x)) and x>=0 for x in N):
        raise ValueError("five nonnegative finite species counts")
    return N[0]-N[1]+N[2]-N[3]


def minimum_charged_fermion_omega_h2(eta_D_today,mchi2=MCHI2_GEV):
    """Strict bound for conserved net DARK U1 per today's SM entropy."""
    if not (math.isfinite(eta_D_today) and mchi2>0):
        raise ValueError("finite asymmetry and positive lighter charged mass")
    return abs(eta_D_today)*mchi2*S0_CM3/RHOC_OVER_H2_GEV_CM3


def maximum_asymmetry_given_fraction(fraction=1.,mchi2=MCHI2_GEV,
                                     omega_C_h2=RELIC_OBS_OMEGA_C_H2):
    """Necessary (not sufficient) allowable |Qcomoving|/S_SM,0."""
    if not (0<fraction<=1 and mchi2>0 and omega_C_h2>0):
        raise ValueError("physical fraction and positive masses")
    return fraction*omega_C_h2*RHOC_OVER_H2_GEV_CM3/(mchi2*S0_CM3)


def verify_exact_event_charge_conservation():
    """Physical integer stoichiometric changes for all dark+SM channels."""
    q=(1,-1,1,-1,0)
    # Δn in (χ1,antiχ1,χ2,antiχ2,φ):
    events={
        "chi1_decay":(-1,0,+1,0,+1),
        "anti_chi1_decay":(0,-1,0,+1,+1),
        "chi1_antichi1_to_phiphi":(-1,-1,0,0,+2),
        "chi2_antichi2_to_phiphi":(0,0,-1,-1,+2),
        "chi1_antichi1_to_SM":(-1,-1,0,0,0),
        "chi2_antichi2_to_SM":(0,0,-1,-1,0),
        "vector_dark_SM_elastic":(0,0,0,0,0),
    }
    charges={name:sum(x*y for x,y in zip(q,dn)) for name,dn in events.items()}
    if any(charges.values()):
        raise ArithmeticError("Broken exact dark Dirac U1 by claimed event")
    return charges


def audit():
    return dict(schema_version=1,
                scope="strict algebraic dark U1 asymmetry floor for neutral real phi and common-charge vector portal; no thermal history, no full FRW/Ωphi",
                conserved_integer_dark_U1_events=verify_exact_event_charge_conservation(),
                stable_mass_lower_GeV=MCHI2_GEV,
                omega_charged_min_per_abs_eta_today=minimum_charged_fermion_omega_h2(1.),
                abs_charge_per_SM_entropy_max_for_100pct_DM=
                      maximum_asymmetry_given_fraction(1.),
                abs_charge_per_SM_entropy_max_for_1pct_DM=
                      maximum_asymmetry_given_fraction(.01),
                explicit_example_eta_D_1e_minus10_omega=
                      minimum_charged_fermion_omega_h2(1e-10),
                exact_bound_independent_of_thermal_MB_chemical_equilibrium=True,
                caveat="Today's SM entropy yield used for normalization; if early SM entropy was heated, early yield differs. No neutral phi lower bound, no physical measured asymmetry for dark sector, no full model exclusion",
                class_ready=False)


if __name__=="__main__":
    print(json.dumps(audit(),indent=2))
