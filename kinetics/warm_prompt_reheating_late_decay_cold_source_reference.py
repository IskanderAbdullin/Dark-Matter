"""Cold vacuum-daughter fraction upper for prompt heavy decays after the cold anchor."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad

from archived_collisionless_streaming_reference import Background
from conditional_relic_normalization_reference import conversion
from three_species_decay_onsager_reference import m, mp, M
from warm_prompt_reheating_expanding_slow_reference import visible


def cold_fraction(y,e,p,cut):
    low=mp*np.cosh(y-np.arccosh(e/mp))
    width=2*p*np.sinh(y)
    return np.clip((cut-low)/width,0.,1.) if y else float(e<cut)


def run(eos,iso,ext):
    i=json.loads(iso.read_text())
    joule=json.loads(ext.read_text())['anchors'][-1]
    st=i['anchors'][-1]
    T=st['visible_T']
    assert T==joule['T_visible_GeV']
    _,s=visible(Background(eos),T)
    bg=Background(eos)
    a=bg.scale(T)
    e=(M*M-m*m+mp*mp)/(2*M)
    p=np.sqrt((e-mp)*(e+mp))
    ystar=np.arccosh(e/mp)
    rho_cap=joule['remaining_energy_ceiling']
    # Grant every particle initially representable at the heavy rest mass.
    # Bound any generation of heavy particles by the prior gross input.
    nh=st['initial_dark_energy_density']/(M*s)
    nvh=joule['additional_comoving_pair_energy_upper']/(M*s)
    nrh=joule['regeneration_weighted_yield_upper']
    all_heavy=nh+nvh+nrh
    target=.06/(mp*conversion(eos)['Omega_h2_per_GeV_times_yield'])
    assert all_heavy>0
    rows=[]
    for q in (.01,.1,1.,10.):
        pc=q*bg.T0/a
        cut=np.hypot(mp,pc)
        ymin=ystar-np.arccosh(cut/mp)
        assert ymin>0
        min_parent=M*np.cosh(ymin)
        frac=pc*pc/(2*p*(p+np.sqrt((e-cut)*(e+cut))))
        yopt=np.arccosh(e/cut)
        assert abs(cold_fraction(yopt,e,p,cut)/frac-1)<5e-14
        ys=np.linspace(ymin,2*ystar,300)
        sample=np.array([cold_fraction(y,e,p,cut) for y in ys])
        assert max(sample)<frac*(1+2e-14)
        for y in (ymin+.1,ystar,ystar+.5):
            yc=(cut/(np.cosh(y)+0.)-e)/(p*np.tanh(y))
            true=quad(lambda c:.5,-1.,max(-1.,min(yc,1.)),epsabs=1e-17)[0]
            assert abs(cold_fraction(y,e,p,cut)-true)<2e-15
        rows.append(dict(q_cut_over_Tgamma0=q,lab_p_cut_GeV=pc,
                         min_parent_energy_GeV=min_parent,
                         one_decay_cold_probability_upper=frac,
                         sampled_max_one_decay_probability=max(sample),
                         integrated_cold_scalar_yield_upper=all_heavy*frac,
                         ratio_to_reference_half_abundance=all_heavy*frac/target))
    return dict(scope="One collisionless scalar per unpolarized vacuum heavy decay after "
                      "the modeled 80 keV anchor; no angular Bose enhancement, final-state "
                      "blocking, late cooling, inverse prompt decays, additional heavy sources or coherent field. "
                      "The bound assumes the previously allowed lepton energy and dark recycling.",
                T_visible_GeV=T, energy_budget_GeV4=rho_cap,
                rest_frame_scalar_E_GeV=e,admitted_heavy_yield_upper=all_heavy,
                granted_initial_heavy_yield=nh, regeneration_yield_upper=nrh,
                direct_visible_heavy_yield_upper=nvh,examples=rows)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('eos','isolated','external'):
        p.add_argument(name,type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.isolated,a.external),indent=2,allow_nan=False))
