"""Audit opposite symmetric orderings at the same first-window elastic resolution."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

from early_pole_joined_temporal_refinement_reference import conservation


def compare(profiles, coarse_macro_steps=16):
    if coarse_macro_steps not in (16,32):
        raise ValueError('Only the checked 16x4 and 32x2 held layouts are supported.')
    controls=('projection','duration','grid','aperture_grid','hard_edge_sets',
              'strength','same_soft_strength','compton_strength','conversion_strength',
              'screening_update','method','pair','response_smear')
    ref=profiles['fine64']
    duration=ref[0]['duration']
    t=np.linspace(0.,duration,17)
    out={}
    for profile,n in (('elastic_centered',coarse_macro_steps),('slow_centered',coarse_macro_steps),('fine64',64)):
        arms=profiles[profile]
        settings=[]
        for a,b,label in zip(arms,ref,('native','local-daughter'),strict=True):
            if a.get('prepared_charged_swap','native')!=label:
                raise ValueError('Native and aperture-local arms must retain their labels.')
            for k,v in (('scalar_bins',8),('aperture_subdivision',1),('witness_mix',0.),
                        ('response_radial_order',None)):
                if a.get(k,v)!=b.get(k,v):
                    raise ValueError('A secondary time/geometry control differs.')
            if any(k in a for k in ('macro_step_weights','macro_exposure')):
                raise ValueError('The dyadic comparison was defined on uniform grids.')
            hist=np.asarray(a['state_history'])
            if hist.shape!=(n+1,) or not np.all(np.isfinite(hist)&(hist>0.)):
                raise ValueError('Require full admissible cold history.')
            if hist[0]!=b['state_history'][0] or hist[-1]-hist[0]!=a['scalar_gain']:
                raise ValueError('Unexpected change in seed or stored gain.')
            if any(a[k]!=b[k] for k in controls) or a['steps']!=n or a['screening_update']!='frozen':
                raise ValueError('Different original held control.')
            for k in ('hard_moments','cold_bins','instantaneous_source_by_q'):
                np.testing.assert_array_equal(a['snapshots'][0][k],b['snapshots'][0][k])
            if a.get('outer_order','slow-elastic-slow')!=(
                    'elastic-slow-elastic' if profile=='slow_centered' else 'slow-elastic-slow'):
                raise ValueError('Wrong outer order metadata.')
            conform=conservation(a)
            if n<64:
                if (a['elastic_substep_energy_budget']!=5e-15 or len(a['elastic'])!=n
                        or a.get('elastic_substeps_per_macro_step')!=64//n):
                    raise ValueError('Unexpected subcycle count or energy gate.')
                for v in a['elastic']:
                    steps=v['elastic_substeps']
                    if len(steps)!=64//n or any(abs(z['rho_drift'])>=5e-15 for z in steps):
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
    p.add_argument('--coarse-macro-steps',type=int,choices=(16,32),default=16)
    a=p.parse_args()
    print(json.dumps(compare({k:[json.loads(Path(f).read_text()) for f in fs]
                              for k,fs in vars(a).items() if k!='coarse_macro_steps'},
                             a.coarse_macro_steps),indent=2,allow_nan=False))
