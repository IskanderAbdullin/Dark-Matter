"""Audit opposite symmetric orderings at the same first-window elastic resolution."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

from early_pole_joined_temporal_refinement_reference import conservation


def compare(profiles):
    controls=('projection','duration','grid','aperture_grid','hard_edge_sets',
              'strength','same_soft_strength','compton_strength','conversion_strength',
              'screening_update','method','pair','response_smear')
    ref=profiles['fine64']
    duration=ref[0]['duration']
    t=np.linspace(0.,duration,17)
    out={}
    for profile,n in (('elastic_centered',16),('slow_centered',16),('fine64',64)):
        arms=profiles[profile]
        settings=[]
        for a,b in zip(arms,ref,strict=True):
            if any(a[k]!=b[k] for k in controls) or a['steps']!=n or a['screening_update']!='frozen':
                raise ValueError('Different original held control.')
            for k in ('hard_moments','cold_bins','instantaneous_source_by_q'):
                np.testing.assert_array_equal(a['snapshots'][0][k],b['snapshots'][0][k])
            if a.get('outer_order','slow-elastic-slow')!=(
                    'elastic-slow-elastic' if profile=='slow_centered' else 'slow-elastic-slow'):
                raise ValueError('Wrong outer order metadata.')
            conform=conservation(a)
            if n==16:
                if a['elastic_substep_energy_budget']!=5e-15 or len(a['elastic'])!=16:
                    raise ValueError('Unexpected subcycle count or energy gate.')
                for v in a['elastic']:
                    steps=v['elastic_substeps']
                    if len(steps)!=4 or any(abs(z['rho_drift'])>=5e-15 for z in steps):
                        raise ValueError('Missing or invalid inner elastic ledger.')
                    np.testing.assert_allclose([z['h'] for z in steps],duration/64,rtol=1e-15,atol=0.)
                    for k in ('rho_drift','entropy_gain'):
                        observed=v[k]
                        target=math.fsum(z[k] for z in steps)
                        if not math.isclose(observed,target,rel_tol=2e-15,abs_tol=1e-28):
                            raise ValueError('Aggregate ledger does not equal paired leaf sum.')
                peak=max(abs(z['rho_drift']) for v in a['elastic'] for z in v['elastic_substeps'])
                conform['peak_abs_macro_signed_elastic_defect']=conform.pop('largest_step_energy_defect')
                conform['peak_abs_individual_elastic_defect']=peak
            else:
                if len(a['elastic'])!=64 or a['elastic_energy_budget']!=5e-15:
                    raise ValueError('Missing fine elastic ledger.')
                np.testing.assert_allclose([z['h'] for z in a['elastic']],duration/64,rtol=1e-15,atol=0.)
            settings.append(dict(conform,slow_rhs_calls=a['slow_rhs_calls']))
        x,y=(np.asarray(a['state_history']) for a in arms)
        if x.shape!=y.shape or len(x)!=n+1 or x[0]!=y[0]:
            raise ValueError('Mismatched cold history.')
        d=((y-y[0])-(x-x[0]))[::n//16]
        out[profile]=dict(final_mean_shift=float(d[-1]/duration),first_block_paired_count=float(d[1]),
                          paired_count_history=d.tolist(),conservation=settings)
    a,b,c=(out[k]['final_mean_shift'] for k in ('elastic_centered','slow_centered','fine64'))
    return dict(cases=out,exposure=t.tolist(),
                elastic_centered_minus_fine=a-c,slow_centered_minus_fine=b-c,
                mixed_orders_average_minus_fine=(a+b)/2-c,
                scope='Matched same inner widths, different outer order only; finite diagnostic without order claim.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('elastic_centered','slow_centered','fine64'):
        p.add_argument('--'+k,nargs=2,required=True,help='Native then local output.')
    a=p.parse_args()
    print(json.dumps(compare({k:[json.loads(Path(f).read_text()) for f in fs]
                              for k,fs in vars(a).items()}),indent=2,allow_nan=False))
