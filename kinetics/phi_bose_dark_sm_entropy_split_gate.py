"""Conditional exact two-reservoir comoving entropy temperature-ratio identity.

When a dark sector containing χ1 Dirac+anti (g=4), χ2 Dirac+anti
(g=4) and one real φ (g=1) once shares Tdec with the SM, but then
their two comoving entropies become SEPARATELY conserved, the Tdark/TSM
ratio is determined only by THEIR OWN separately evolving effective
entropy degrees g*S(T):

  xi = [gSM,S(Tsm)/gSM,S(Tdec) *
        gD,S(Tdec)/gD,S(Tdark)]^(1/3).

The actual dark sector's gD,S changes as massive χ1/χ2 become NR,
decay, pair-annihilate or acquire chemical potentials; if freezeout
is nonthermal there is NO single gD,S(T) description.

This is an exact algebraic adiabatic limit, with illustrative
gSM,dec=61.75 (QCD crossover approximate), gSM,MeV=10.75.
For three initial relativistic dark fields, gD,dec=1+7/8*8=8.
Two endpoint scenarios: no dark entropy release (gD remains 8)
gives xi=(10.75/61.75)^(1/3)~0.558; complete transfer into
one thermal real scalar (gDlate=1) gives xi~1.117.
The latter is INCOMPATIBLE with an appreciable relic χ2, so cannot
be used as a physical computed temperature after actual freezeout.
The difference is a QUANTITATIVE uncertainty warning, not a claimed
bound for arbitrary nonequilibrium spectra, reheating, SM entropy
injection, or actual g*(T).

No CLASS, relic abundance, BEC or galaxy-halo evidence follows.
"""
from __future__ import annotations
import json
import math


def relativistic_dark_gstar_s():
    """One real scalar + two distinct four-state Dirac fermions."""
    return 1. + 7/8*4 + 7/8*4


def decoupled_entropy_temperature_ratio(gsm_before=61.75,gsm_after=10.75,
                                        gd_before=8.,gd_after=8.):
    """Exact a³ s(g,T) reservoir identity for separate adiabatic sectors."""
    if not all(math.isfinite(x) and x>0 for x in
               (gsm_before,gsm_after,gd_before,gd_after)):
        raise ValueError("Positive finite entropy DOF required")
    return ((gsm_after/gsm_before)*(gd_before/gd_after))**(1/3)


def benchmark():
    gd=relativistic_dark_gstar_s()
    no_dark_heat=decoupled_entropy_temperature_ratio(gd_before=gd,
                                                      gd_after=gd)
    all_into_phi=decoupled_entropy_temperature_ratio(gd_before=gd,
                                                     gd_after=1.)
    return dict(schema_version=1,
                model="one real phi plus two Dirac χ species with anti states; 1+7/8*(4+4)=8 initial thermal entropy DOF",
                initial_relativistic_dark_gstar_s=gd,
                prescribed_external_SM_gstar_s_before=61.75,
                prescribed_external_SM_gstar_s_MeV=10.75,
                no_dark_internal_entropy_release_xi=no_dark_heat,
                complete_hypothetical_dark_reheating_into_phi_xi=all_into_phi,
                ratio_of_two_cases=all_into_phi/no_dark_heat,
                SM_QCD_gstar_time_dependence_not_computed=True,
                dark_Bose_Fermi_chemical_entropy_not_computed=True,
                actual_chi2_relic_compatible_with_gD_after_1=False,
                no_common_dark_SM_temperature_validated=True,
                class_ready=False)


if __name__=="__main__":
    print(json.dumps(benchmark(),indent=2))
