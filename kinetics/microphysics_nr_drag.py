"""Controlled NONRELATIVISTIC heavy chi kinetic drag from microscopic recoil.

In a dilute nondegenerate chi gas with M/T >> 1, small momentum
transfers allow a Fokker-Planck/Langevin description
 dp_i = -eta p_i dt + sqrt(kappa) dW_i
where <dp_i^2>/dt=kappa, Maxwell equilibrium variance M*T and
the Einstein relation eta=kappa/(2*M*T).  At p=0, energy gain
Q_rest=3*kappa/(2*M), so eta=Q_rest/(3*T).

Our Q_rest includes final-state electron Pauli and full vector-contact
two-body recoil. This is leading small-transfer matching, NOT an exact
relativistic BGK gamma_el(p,T); requires M/T>>1, valid contact EFT,
negligible chi Pauli blocking and a near-thermal bath. Do not apply to
early BEST M/T=0.1 (relativistic parent).
"""
from __future__ import annotations
import math
from microphysics_recoil_heat import (
    recoil_heating_rate_chi_at_rest,heavy_chi_recoil_asymptote)


def heavy_nr_fokker_planck_coefficients(M,T,G,g_bath=4,
                                      min_mass_over_temperature=100.):
    vals=(M,T,G,g_bath,min_mass_over_temperature)
    if (not all(math.isfinite(v) for v in vals) or M<=0 or T<=0 or
            G<0 or g_bath<=0 or min_mass_over_temperature<50):
        raise ValueError("Finite M,T,g>0, G>=0 and conservative validity threshold required")
    if M/T<min_mass_over_temperature:
        raise ValueError("Relativistic/transition regime: heavy-chi Fokker-Planck invalid")
    Q=recoil_heating_rate_chi_at_rest(M,T,G,g_bath,blocking=True)
    kappa=2*M*Q/3
    eta=Q/(3*T)
    Q_asym=heavy_chi_recoil_asymptote(M,T,G,g_bath)
    eta_asym=Q_asym/(3*T)
    return {
      "M_over_T":M/T,
      "Q_rest_GeV2_per_chi":Q,
      "momentum_diffusion_kappa_GeV3":kappa,
      "leading_drag_eta_GeV":eta,
      "leading_drag_heavy_limit_GeV":eta_asym,
      "heavy_limit_relative_deviation":eta/eta_asym-1 if eta_asym>0 else 0.,
      "fluctuation_dissipation_Einstein_relation":True,
      "validity":"leading heavy nondegenerate chi, M/T>>1, near-equilibrium "
                 "electron bath, contact EFT, neglect chi blocking; gamma "
                 "is NOT valid at the relativistic BEST benchmark M/T=0.1",
      "full_thermal_collision_operator_validated":False,
      "CLASS_ready":False
    }
