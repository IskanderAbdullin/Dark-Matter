"""Entropy-per-particle gate for condensation in an isolated fixed-count ideal scalar gas."""
import argparse
import gzip
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.special import zeta
from three_species_decay_onsager_reference import mp


def thermodynamic_entropy(r,phase):
    # p/T = u; x=(E-m)/T is evaluated without subtracting nearly equal numbers.
    def base(u,typ):
        x=u*u/(np.hypot(u,r)+r)
        Z=x-np.log(phase)
        f=np.exp(-Z)/(-np.expm1(-Z))
        if typ=='n':
            c=f
        elif typ=='rho':
            c=np.hypot(u,r)*f
        elif typ=='p':
            c=u*u/(3*np.hypot(u,r))*f
        elif typ=='entropy':
            c=-np.log(-np.expm1(-Z))+Z*f
        return u*u*c
    out={k:quad(lambda u:base(u,k),0.,200.+40*np.sqrt(max(1.,2*r)),
                epsabs=1e-12,epsrel=5e-10)[0] for k in ('n','rho','p','entropy')}
    S=out['entropy']/out['n']
    check=(out['rho']+out['p'])/out['n']-r-np.log(phase)
    assert abs(check/S-1.)<2e-9
    return S


def run(history,saturated):
    path=Path(history)
    with gzip.open(path,'rt') as f:
        ans=json.load(f)['arms'][-1]['trajectory']
    zero=(2*np.pi**4/45)/zeta(3)
    nr=2.5*zeta(2.5)/zeta(1.5)
    assert abs(thermodynamic_entropy(0.,1.)/zero-1)<2e-10
    grid=[]
    for r in (0.,.01,.1,1.,3.,10.,100.,10000.):
        s=thermodynamic_entropy(r,1.)
        if grid:
            assert grid[-1]['s_per_particle']>s>nr
        grid.append(dict(m_over_T=r,s_per_particle=s))
    assert abs(grid[-1]['s_per_particle']/nr-1)<2e-4
    rows=[]
    for r in ans:
        T=r['light_scalar_T'][1]
        eta=r['phase_factors'][1]
        S=thermodynamic_entropy(mp/T,eta)
        rows.append(dict(visible_T_GeV=r['T_visible'],s_per_particle=S,
                         critical_s_per_excitation_at_same_x=thermodynamic_entropy(mp/T,1.),
                         entropy_fraction_to_remove_at_fixed_count_before_saturation=1-zero/S,
                         necessary_count_multiplier_at_fixed_entropy_before_saturation=S/zero))
    assert min(r['s_per_particle'] for r in rows)>zero
    c=float(next(a['minimum_condensate_fraction_if_excited_component_is_normal_thermal']
            for a in json.loads(saturated.read_text())['cold_population_condensate_examples']
            if a['q_cut_over_Tgamma0']==.1 and a['target_cold_fraction']==.5))
    return dict(scope="Ideal scalar thermodynamics with conserved scalar number, no heat or entropy "
                      "exchange or phase-dependent new physics after each hypothesized anchor. "
                      "Boltzmann entropy is not assumed to decrease; number or entropy transfers "
                      "could evade this comparison. The modeled thermal states are illustrative.",
                saturated_UL_critical_entropy_per_excitation=zero,
                saturated_NR_critical_entropy_per_excitation=nr,
                saturated_tests=grid,modeled_anchor_rows=rows,
                optimistic_maximum_total_s_per_particle_for_point1_half_cold_exemplar=(1-c)*zero,
                entropy_loss_if_half_target_inherits_80keV_entropy_per_particle=1-(1-c)*zero/rows[-1]['s_per_particle'])


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('history',type=Path)
    p.add_argument('saturated',type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.history,a.saturated),indent=2,allow_nan=False))
