"""Compare elastic subcycling to refining the complete prepared split."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_joined_temporal_refinement_reference import conservation


def small_ledger(z):
    n=z['steps']
    k=z['elastic_substeps_per_macro_step']
    macro=z['elastic']
    if k!=4 or len(macro)!=n or z['elastic_substep_energy_budget']!=5e-15:
        raise ValueError('Unexpected subcycle count or budget.')
    all_steps=[]
    for m in macro:
        steps=m['elastic_substeps']
        if len(steps)!=k:
            raise ValueError('Incomplete local elastic receipt.')
        np.testing.assert_allclose(fsum(x['h'] for x in steps),z['duration']/n,rtol=1e-15,atol=0.)
        np.testing.assert_equal(m['rho_drift'],fsum(x['rho_drift'] for x in steps))
        np.testing.assert_equal(m['entropy_gain'],fsum(x['entropy_gain'] for x in steps))
        all_steps.extend(steps)
    peak=max(abs(x['rho_drift']) for x in all_steps)
    if peak>=z['elastic_substep_energy_budget']:
        raise ValueError('Peak recorded small step exceeded its gate.')
    result=conservation(z)
    result['peak_abs_macro_signed_elastic_defect']=result.pop('largest_step_energy_defect')
    result['elastic_substep_energy_budget']=z['elastic_substep_energy_budget']
    return dict(result,peak_substep_energy_defect=peak,number_of_elastic_steps=len(all_steps))


def compare(base, subcycled, finer, exact):
    controls=('projection','duration','grid','aperture_grid','hard_edge_sets','pair','method',
              'strength','same_soft_strength','compton_strength','conversion_strength',
              'screening_update','response_smear')
    signed={}
    all_ledgers={}
    for profile,arm in (('base16',base),('subcycled16',subcycled),('fine64',finer)):
        expected_steps=64 if profile=='fine64' else 16
        for r,z in zip(base,arm,strict=True):
            if any(r[k]!=z[k] for k in controls):
                raise ValueError('Mismatched physical input and controls.')
            if z['screening_update']!='frozen' or z['steps']!=expected_steps:
                raise ValueError('Mismatched response or outer clock.')
            for k in ('cold_bins','hard_moments','instantaneous_source_by_q'):
                np.testing.assert_array_equal(r['snapshots'][0][k],z['snapshots'][0][k])
        a,b=(np.asarray(z['state_history']) for z in arm)
        delta=((b-b[0])-(a-a[0]))[::expected_steps//16]
        signed[profile]=delta
        all_ledgers[profile]=[(small_ledger if profile=='subcycled16' else conservation)(x) for x in arm]
    if exact[1].get('prepared_charged_swap')!='local-daughter':
        raise ValueError('Expected aperture-local daughter chart.')
    for x in exact:
        if x['duration']!=base[0]['duration'] or x['report_panels']!=16:
            raise ValueError('Mismatched simultaneous-control duration or sample grid.')
    a,b=(np.asarray(z['cold_count_history']) for z in exact)
    signed['unsplit']=(b-b[0])-(a-a[0])
    times=exact[0]['history_exposure']
    np.testing.assert_array_equal(times,np.linspace(0,base[0]['duration'],17))
    mean={k:float(v[-1]/times[-1]) for k,v in signed.items()}
    common_width=base[0]['duration']/64
    for z in finer:
        if len(z['elastic'])!=64:
            raise ValueError('Expected every fine elastic receipt.')
        for q in z['elastic']:
            np.testing.assert_allclose(q['h'],common_width,rtol=1e-15,atol=0.)
    same_width_difference=signed['subcycled16']-signed['fine64']
    # Conditional fixed-inner-width cancellation of a leading outer H**2 term.
    outer_extrapolated=(16*signed['fine64']-signed['subcycled16'])/15
    outer=dict(inner_step=common_width,mean_difference_16x4_minus_64x1=float(
                   same_width_difference[-1]/times[-1]),
               maximum_absolute_cumulative_difference=float(np.max(abs(same_width_difference))),
               extrapolated_mean=float(outer_extrapolated[-1]/times[-1]),
               extrapolated_mean_error_vs_unsplit=float((outer_extrapolated[-1]-signed['unsplit'][-1])/times[-1]),
               extrapolated_cumulative_effect=outer_extrapolated.tolist(),
               assumption='Optional smooth outer quadratic diagnostic at fixed inner width, not an order proof.')
    full_change=mean['fine64']-mean['base16']
    return dict(final_mean_shift=mean,exposure=times,
                cumulative_effect={k:x.tolist() for k,x in signed.items()},
                elastic_only_share_of_total_observed_refinement=(
                    (mean['subcycled16']-mean['base16'])/full_change if full_change else None),
                common_inner_width=outer,conservation=all_ledgers,
                scope='Finite held graph; inner subcycle changes elastic discretization at a fixed outer split width.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('base','subcycled','finer','exact'):
        p.add_argument('--'+k,nargs=2,required=True,help='Native and aperture-local outputs, in order.')
    a=p.parse_args()
    print(json.dumps(compare(*([json.loads(Path(p).read_text()) for p in paths]
                              for paths in vars(a).values())),indent=2,allow_nan=False))
