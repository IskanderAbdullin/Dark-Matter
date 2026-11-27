"""Measure what refining only the first original interleaving interval recovers."""
import argparse
import json
from pathlib import Path

import numpy as np

from early_pole_joined_temporal_refinement_reference import conservation


def compare(coarse,hybrid,fine):
    controls=('projection','duration','grid','aperture_grid','hard_edge_sets','strength',
              'same_soft_strength','compton_strength','conversion_strength','screening_update',
              'method','pair','response_smear')
    durations=np.arange(65)*fine[0]['duration']/64
    results={}
    for name,pair in (('coarse16x4',coarse),('hybrid19',hybrid),('fine64x1',fine)):
        ledgers=[]
        for x,ref in zip(pair,fine,strict=True):
            if any(x[k]!=ref[k] for k in controls) or x['screening_update']!='frozen':
                raise ValueError('Mismatched setup.')
            for k in ('hard_moments','cold_bins','instantaneous_source_by_q'):
                np.testing.assert_array_equal(x['snapshots'][0][k],ref['snapshots'][0][k])
            c=conservation(x)
            if name=='hybrid19':
                expected=[1]*4+[4]*15
                if x['macro_step_weights']!=expected or x['steps']!=19:
                    raise ValueError('Wrong hybrid schedule.')
                np.testing.assert_allclose(x['macro_exposure'],
                                           durations[np.r_[0,np.cumsum(expected)]],atol=0.,rtol=2e-15)
                np.testing.assert_array_equal(x['state_history'][:5],ref['state_history'][:5])
                allsteps=[s for k in x['elastic'] for s in k['elastic_substeps']]
                if len(allsteps)!=64 or x['elastic_substep_energy_budget']!=5e-15:
                    raise ValueError('Changed inner step count or tolerance.')
                peak=max(abs(s['rho_drift']) for s in allsteps)
                if peak>=x['elastic_substep_energy_budget']:
                    raise ValueError('A local elastic gate failed.')
                np.testing.assert_allclose([s['h'] for s in allsteps],ref['duration']/64,atol=0.,rtol=1e-15)
                c['peak_abs_macro_signed_elastic_defect']=c.pop('largest_step_energy_defect')
                c['peak_abs_individual_elastic_defect']=peak
            ledgers.append(c)
        a,b=pair
        results[name]=dict(mean_shift=(b['scalar_gain']-a['scalar_gain'])/a['duration'],
                           passes_first_interval_replay=name=='hybrid19',
                           outer_intervals=a['steps'],closure=ledgers)
    c,m,f=(results[k]['mean_shift'] for k in ('coarse16x4','hybrid19','fine64x1'))
    def signed(pair, idx):
        a,b=pair
        return (b['state_history'][idx]-b['state_history'][0])-(a['state_history'][idx]-a['state_history'][0])
    observed_startup_change=signed(fine,4)-signed(coarse,1)
    prediction=c+observed_startup_change/fine[0]['duration']
    return dict(cases=results,percent_gap_recovered=100*(m-c)/(f-c),
                hybrid_minus_fine_mean=m-f,
                pasted_startup_mean_prediction=prediction,
                actual_minus_pasted_startup_prediction=m-prediction,
                predicted_percent_gap_recovered_from_old_runs=100*(prediction-c)/(f-c),
                scope='Matches first interval sample exactly before continuing at coarser outer width; no general error bound.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('coarse','hybrid','fine'):
        p.add_argument('--'+k,nargs=2,required=True,help='Native then local.')
    a=p.parse_args()
    print(json.dumps(compare(*([json.loads(Path(f).read_text()) for f in paths]
                                 for paths in vars(a).values())),indent=2,allow_nan=False))
