"""Numerical counterexample: exact quantum detailed balance with wrong Bose cell number.

Every 1<->2 Yukawa microevent has F=R at common T, yet collocated
BE(mphi) times the first positive radial energy cell is not the actual
thermal scalar number for mphi << first grid spacing. Not a cosmological
production or physical BEC result.
"""
from __future__ import annotations
import json
from kinetics.four_species_yukawa_shared_lepton import configuration,initial,ledger,rhs
from kinetics.phi_ir_first_cell_artifact import thermal_first_cell_and_total


def compare(bins=64, mphi=86e-9, temperature=1., upper_cut=6.):
    cfg=configuration(bins=bins,mphi=mphi,elastic_cut=upper_cut,
                      elastic_G=0.,elastic_nr=3,elastic_nmu=3,
                      elastic_nc=6,elastic_nphi=8,
                      yukawa_nr=8,yukawa_nangle=8)
    state=initial(cfg,T1=temperature,T2=temperature,
                  Tlep=temperature,Tphi=temperature)
    physical=thermal_first_cell_and_total(
        mphi,upper=float(cfg['yukawa'].ephi[-1]),bins=bins,T=temperature)
    charges=ledger(cfg,state)
    dx,info=rhs(cfg,state)
    return {
        'mphi':mphi,'temperature':temperature,'bins':bins,
        'common_energy_upper':float(cfg['yukawa'].ephi[-1]),
        'grid_phi_number':charges['Nphi'],
        'direct_excited_thermal_phi_number':physical['full_exact_thermal_total'],
        'false_full_thermal_number_ratio':charges['Nphi']/physical['full_exact_thermal_total'],
        'false_first_cell_ratio':physical['first_cell_false_enhancement_ratio'],
        'absolute_yukawa_rate':abs(info['Yukawa_Rphi']),
        'absolute_entropy_production':abs(info['yukawa_entropy']),
        'absolute_energy_balance':abs(info['energy_residual']),
        'max_absolute_charge_defect':max(map(abs,info['charge_residuals'])),
        'discrete_detailed_balance_not_correct_thermal_number':True,
        'not_a_BEC_or_relic':True
    }


def audit():
    cases=[]
    for m in (.001,86e-9):
        for bins in (16,64,256):
            cases.append(compare(bins,m))
    return {'schema_version':1,'toy_Gev_thermal_equality_only':True,'cases':cases,
            'no_real_scalar_condensate_or_Gamma_el':True,'class_ready':False}


if __name__=='__main__':
    print(json.dumps(audit(),indent=2))
