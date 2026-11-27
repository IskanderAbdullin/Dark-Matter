"""Conditional prior-anchor count edits and late dilution feasibility ledger."""
import argparse
import gzip
import json
from pathlib import Path
import numpy as np

from archived_collisionless_streaming_reference import Background
from conditional_relic_normalization_reference import conversion
from three_species_decay_onsager_reference import m, mp, M
from warm_prompt_reheating_expanding_slow_reference import visible
from warm_prompt_reheating_cap_free_thresholds_reference import (
    coefficients, slab_passage_lower,
)


def passage_for_q_interval(lo, hi, T, n0, rho, eta, norm, target):
    """Turn the band of possible initial counts into one higher loss envelope."""
    cut=target/(norm*lo)
    if cut>=1:
        return 0.
    energy=lambda q:rho-eta*(1-q)*m*n0
    co=coefficients(T,n0*hi,energy(hi),.0049)
    # r(q) = eta + (r(1)-eta)/q, decreasing since r(1)>eta.
    co['r']=energy(lo)/(m*n0*lo)
    # All other coefficients are nonnegative constants times positive powers
    # of n or rho, including the optional per-pair vector ceiling.
    return slab_passage_lower(cut,co,1024)


def read_gzip(path):
    with gzip.open(path,'rt') as file:
        return json.load(file)


def run(eos,sample,isolated,budget):
    a=read_gzip(sample)['sampled_anchors'][-1]
    ref=json.loads(isolated.read_text())
    ene=ref['anchors'][-1]
    vis=json.loads(budget.read_text())['anchors'][-1]
    T,n0,rho=a['T_visible_GeV'],a['n0'],a['rho_cap']
    assert T==ene['visible_T']==vis['T_visible_GeV']
    assert abs(rho/vis['remaining_energy_ceiling']-1)<2e-14
    assert rho>m*n0
    conv=conversion(eos)['Omega_h2_per_GeV_times_yield']
    bg=Background(eos)
    _,s=visible(bg,T)
    norm=m*conv*n0/s
    H=ene['photon_H_floor_at_anchor']
    target_total=.120
    target_half=.060
    y_total_sc=target_half/(mp*conv)
    q0=ene['initial_weighted_yield']
    max_recycle=lambda r:ref['bounded_heavy_coefficient']*r*r/(M*M*H*s)
    y_preexisting=q0+vis['direct_weighted_yield_upper']+max_recycle(rho)
    assert abs(y_preexisting*mp*conv/vis['total_scalar_omega_upper_including_leptons']-1)<1e-14
    change=[]
    for eta in (0.,.5,1.):
        for target in (target_total,target_half):
            controls=[]
            for N in (256,1024):
                q=np.linspace(.001,1.,N+1)
                bad=[]
                for lo,hi in zip(q[:-1],q[1:],strict=True):
                    p=passage_for_q_interval(lo,hi,T,n0,rho,eta,norm,target)
                    if p<=1.:
                        bad.append(hi)
                cut=max(bad) if bad else q[0]
                assert cut<1.
                controls.append(float(cut))
            assert controls[-1]<=controls[0]+1e-8
            change.append(dict(anchored_stable_target=target,energy_export_fraction_per_removed_stable_rest_mass=eta,
                               remaining_initial_stable_fraction_threshold_for_obstruction=controls[-1],
                               necessary_earlier_stable_fraction_removed_to_escape_obstruction=1-controls[-1],
                               fixed_count_band_checks=controls))
    dilution=[]
    for v in (target_total,target_half):
        dil=a['omega_h2_lower']/v
        sc=y_preexisting/dil*mp*conv
        dilution.append(dict(stable_goal=v,minimum_pure_late_entropy_multiplier=dil,
                             scalar_omega_upper_after_dilution=sc,
                             missing_scalar_omega_to_reference_half=max(0.,target_half-sc)))
    return dict(scope="Only algebraic anchored moment modifications; not a realizable thermal or "
                      "production history. Assumes same post-anchor cap-free channel and ordinary "
                      "cosmology except the separately labelled after-reaction dilution scenario. "
                      "All quoted changes are necessary to escape a comparison, not sufficient.",
                anchor_visible_T_GeV=T,modified_count_scenarios=change,
                scalar_goal_yield=y_total_sc,scalar_max_before_any_new_channel_yield=y_preexisting,
                scalar_deficit_to_half_without_dilution=max(0.,y_total_sc-y_preexisting),
                pure_dilution_scenarios=dilution)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('eos','sample','isolated','budget'):
        p.add_argument(name,type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.sample,a.isolated,a.budget),indent=2,allow_nan=False))
