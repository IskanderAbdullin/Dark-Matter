"""Conditional SU(2)L completion of an electron-vector chi current portal.

L_eff=(chibar gamma chi)[G_L(Lbar_e gamma L_e)+G_R(eRbar gamma eR)].
L_e=(nu_eL,eL). In a massless common-temperature bath, one chiral
lepton and antilepton carry two helicity/charge states total.

For an unpolarized Dirac chi vector current the parity-odd chiral
lepton trace vanishes contracted with the symmetric dark spin tensor;
one chiral lepton state has the same angular |M|² kernel as a pure
vector Dirac charged-lepton helicity state.
Conditional tree-level bath weights:
eL+anti 2 GL²; eR+anti 2 GR²; nuL+antinuR 2 GL².

These weights assume the same thermal distributions for all species;
they are NOT UV completion, EW thermal masses or a neutrino collision
operator with distinct species temperatures, flavour effects or bounds.
"""
from __future__ import annotations
import math


def same_temperature_massless_lepton_weights(GL,GR):
    if not all(math.isfinite(float(x)) for x in (GL,GR)):
        raise ValueError("finite real chiral couplings required")
    left=2*GL*GL
    right=2*GR*GR
    neutrino=2*GL*GL
    electron=left+right
    full=electron+neutrino
    return {
      "G_e_vector":(GL+GR)/2,
      "G_e_axial":(GR-GL)/2,
      "electron_left_weight":left,
      "electron_right_weight":right,
      "neutrino_left_weight":neutrino,
      "electron_only_weight":electron,
      "leptons_completed_weight":full,
      "ratio_full_over_electron_only":full/electron if electron>0 else None,
      "requires_same_temperature_massless_bath":True,
      "does_not_include_UV_completion_or_SM_rates":True,
    }
