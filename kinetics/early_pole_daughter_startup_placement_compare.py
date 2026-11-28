"""Gate equal-work first-only and second-only refinements against an offline prediction."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_joined_temporal_refinement_reference import conservation


def compare(coarse,first,second,fine,allocation):
    reference=fine[0]
    h=reference['duration']/64
    settings=('projection','duration','grid','aperture_grid','hard_edge_sets','strength','pair','method',
              'same_soft_strength','compton_strength','conversion_strength','screening_update','response_smear')
    controls={}
    for tag,arms,n,block in (('coarse',coarse,16,None),('first_only',first,19,1),
                            ('second_only',second,19,2),('fine',fine,64,None)):
        scores=[]
        for i,(v,anchor) in enumerate(zip(arms,fine,strict=True)):
            if (v['steps']!=n or v['screening_update']!='frozen' or any(v[k]!=reference[k] for k in settings)
                    or v.get('prepared_charged_swap','native')!=('native','local-daughter')[i]
                    or v.get('outer_order','slow-elastic-slow')!='slow-elastic-slow'):
                raise ValueError('Different original held comparison setup.')
            start=v['snapshots'][0]
            for k in ('hard_moments','cold_bins','instantaneous_source_by_q'):
                np.testing.assert_array_equal(start[k],anchor['snapshots'][0][k])
            z=np.asarray(v['state_history'])
            if (z.shape!=(n+1,) or not np.all(np.isfinite(z)&(z>0.))
                    or z[0]!=anchor['state_history'][0] or z[-1]-z[0]!=v['scalar_gain']):
                raise ValueError('Missing admissible matched history or stored scalar gain.')
            stats=conservation(v)
            if block is not None:
                expected=[w for j in range(1,17) for w in ([1]*4 if j==block else [4])]
                if v['macro_step_weights']!=expected or v['elastic_substep_energy_budget']!=5e-15:
                    raise ValueError('Wrong refinement weights or elastic budget.')
                np.testing.assert_allclose(v['macro_exposure'],np.r_[0,np.cumsum(expected)]*h,rtol=2e-15,atol=0.)
                for q,m in zip(v['elastic'],expected,strict=True):
                    zs=q['elastic_substeps']
                    if len(zs)!=m or abs(fsum(x['rho_drift'] for x in zs)-q['rho_drift'])>1e-28:
                        raise ValueError('Bad elastic aggregate.')
                    np.testing.assert_allclose([k['h'] for k in zs],h,atol=0.,rtol=1e-15)
                peak=max(abs(x['rho_drift']) for q in v['elastic'] for x in q['elastic_substeps'])
                if peak>=5e-15:
                    raise ValueError('A fixed-width elastic leaf failed.')
                if block==1:
                    np.testing.assert_array_equal(z[:5],anchor['state_history'][:5])
                else:
                    np.testing.assert_array_equal(z[:2],coarse[i]['state_history'][:2])
                stats['max_abs_leaf_energy_defect']=peak
            scores.append(dict(stats,slow_rhs_calls=v['slow_rhs_calls']))
        controls[tag]=dict(mean_shift=(arms[1]['scalar_gain']-arms[0]['scalar_gain'])/reference['duration'],
                           outer_intervals=n,arm_ledgers=scores)
    means={tag:row['mean_shift'] for tag,row in controls.items()}
    span=means['fine']-means['coarse']
    if span!=allocation['slow_coupled_full_gap'] or span<=0.:
        raise ValueError('A distinct full-slow reference gap was used.')
    expected=allocation['second_only_predeclared_prediction']
    if expected['refined_blocks']!=[2]:
        raise ValueError('Wrong predeclared placement.')
    if means['first_only']!=allocation['previous_hybrid_controls']['1']['observed_mean']:
        raise ValueError('First-only raw outputs do not match the archived measurement.')
    f=(means['second_only']-means['coarse'])/span
    error=means['second_only']-expected['paired_mean']
    criteria=dict(positive_second_correction=means['second_only']>means['coarse'],
                  recovered_fraction_in_predeclared_range=.065<=f<=.085,
                  mean_error_within_predeclared_bound=abs(error)<=1e-7)
    return dict(scope='Only second versus only first of 16 original macro blocks; fixed 64 elastic leaves and full slow field.',
                controlled_profiles=controls,
                first_only_gap_recovered=(means['first_only']-means['coarse'])/span,
                second_only_gap_recovered=f,second_prediction_error=error,
                meets_predeclared_criteria=all(criteria.values()),criteria=criteria)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('coarse','first','second','fine'):
        p.add_argument('--'+k,nargs=2,required=True,help='Native then aperture-local.')
    p.add_argument('--allocation',required=True)
    a=p.parse_args()
    r=lambda path:json.loads(Path(path).read_text())
    print(json.dumps(compare(*([r(f) for f in getattr(a,k)] for k in ('coarse','first','second','fine')),
                             r(a.allocation)),indent=2,allow_nan=False))
