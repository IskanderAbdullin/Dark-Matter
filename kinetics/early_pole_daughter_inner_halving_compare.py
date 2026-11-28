"""Compare one and two elastic substeps at the same 64 held outer steps."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_joined_temporal_refinement_reference import conservation


def compare(coarse,fine,chart):
    controls=('duration','projection','steps','grid','aperture_grid','hard_edge_sets',
              'strength','same_soft_strength','compton_strength','conversion_strength',
              'screening_update','pair','method','response_smear')
    ledgers=[]
    for a,b in zip(coarse,fine,strict=True):
        if any(a[k]!=b[k] for k in controls) or a['steps']!=64 or a['screening_update']!='frozen':
            raise ValueError('Require the same 64 outer split intervals and held response.')
        for key in ('hard_moments','cold_bins','instantaneous_source_by_q'):
            np.testing.assert_array_equal(a['snapshots'][0][key],b['snapshots'][0][key])
        if len(a['elastic'])!=64 or len(b['elastic'])!=64 or b['elastic_substeps_per_macro_step']!=2:
            raise ValueError('Incomplete inner comparison.')
        gate=b['elastic_substep_energy_budget']
        if gate!=a['elastic_energy_budget'] or gate!=5e-15:
            raise ValueError('Mismatched per-elastic-step energy gate.')
        terms=[]
        for q,r in zip(a['elastic'],b['elastic'],strict=True):
            k=r['elastic_substeps']
            if len(k)!=2 or max(abs(x['rho_drift']) for x in k)>=gate:
                raise ValueError('Invalid inner step receipt.')
            for x in k:
                np.testing.assert_allclose(2*x['h'],q['h'],rtol=1e-15,atol=0.)
            np.testing.assert_equal(fsum(x['rho_drift'] for x in k),r['rho_drift'])
            terms.extend(k)
        c=conservation(b)
        c['peak_abs_macro_signed_elastic_defect']=c.pop('largest_step_energy_defect')
        c['peak_abs_individual_elastic_defect']=max(abs(x['rho_drift']) for x in terms)
        ledgers.append([conservation(a),c])
    if chart[1].get('prepared_charged_swap')!='local-daughter':
        raise ValueError('Expected native and local data in order.')
    t=coarse[0]['duration']
    at=np.linspace(0,t,17)
    for q in chart:
        np.testing.assert_array_equal(q['history_exposure'],at)
        if not (abs(q['energy_defect'])<3e-12 and abs(q['entropy_defect'])<3e-11):
            raise ValueError('Unaccepted simultaneous reference.')
    values={}
    for name,parts in (('outer64x1',coarse),('outer64x2',fine),('unsplit',chart)):
        field='cold_count_history' if name=='unsplit' else 'state_history'
        a,b=(np.asarray(x[field]) for x in parts)
        d=(b-b[0])-(a-a[0])
        values[name]=d if name=='unsplit' else d[::4]
    linear=2*values['outer64x2']-values['outer64x1']
    rates={k:float(v[-1]/t) for k,v in values.items()}
    return dict(means=rates,exposure=at.tolist(),cumulative_effect={k:v.tolist() for k,v in values.items()},
                mean_change=float((values['outer64x2']-values['outer64x1'])[-1]/t),
                inner_linear_cancelled_mean=float(linear[-1]/t),
                inner_linear_cancelled_minus_unsplit=float(linear[-1]/t-rates['unsplit']),
                maximum_reference_history_gap_before_after=[
                    float(np.max(abs(values[k]-values['unsplit']))) for k in ('outer64x1','outer64x2')],
                conservation=ledgers,
                scope='Same outer width; cancellation assumes a leading inner term, not an order proof.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('coarse','fine','chart'):
        p.add_argument('--'+k,nargs=2,required=True,help='Native then aperture-local raw outputs.')
    a=p.parse_args()
    print(json.dumps(compare(*([json.loads(Path(r).read_text()) for r in paths]
                                for paths in vars(a).values())),indent=2,allow_nan=False))
