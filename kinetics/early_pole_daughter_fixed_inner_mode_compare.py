"""Resolve held outer-width sensitivity by cold mode at common saved exposures."""
import argparse
import json
from pathlib import Path

import numpy as np

from early_pole_daughter_cold_mode_effect_reference import project
from early_pole_daughter_subcycle_width_compare import small_ledger


def compare(profiles):
    x={}
    controls=('projection','duration','grid','aperture_grid','hard_edge_sets',
              'strength','same_soft_strength','compton_strength','conversion_strength',
              'screening_update','method','pair','response_smear')
    ref=profiles['outer64'][0]
    for key,(a,b) in profiles.items():
        if any(z[k]!=ref[k] for z in (a,b) for k in controls):
            raise ValueError('Different fine physical controls.')
        n=a['steps']
        expected=dict(outer16=16,outer32=32,outer64=64)[key]
        if n!=expected or n!=b['steps'] or a['screening_update']!='frozen':
            raise ValueError('Unexpected outer grid or screening.')
        if key=='outer64':
            if a['elastic_energy_budget']!=5e-15 or b['elastic_energy_budget']!=5e-15:
                raise ValueError('Mismatched fine elastic budget.')
            for z in (a,b):
                if len(z['elastic'])!=64:
                    raise ValueError('Missing fine elastic receipts.')
                np.testing.assert_allclose([q['h'] for q in z['elastic']],ref['duration']/64,
                                           rtol=1e-15,atol=0.)
            x[key]=project(a,b)
        else:
            ledger=[small_ledger(z) for z in (a,b)]
            x[key]=project(dict(a,elastic_energy_budget=None),dict(b,elastic_energy_budget=None))
            x[key]['conservation']=ledger
    for index in (0,1):
        ref=profiles['outer64'][index]
        for z in [v[index] for v in profiles.values()]:
            for f in ('hard_moments','cold_bins','instantaneous_source_by_q'):
                np.testing.assert_array_equal(z['snapshots'][0][f],ref['snapshots'][0][f])
    duration=profiles['outer64'][0]['duration']
    common=sorted(set.intersection(*[set(d['exposure']) for d in x.values()]))
    if not common or common[-1]!=duration:
        raise ValueError('Missing identical terminal sample.')
    aligned={key:np.array([d['mean_mode_effects'][d['exposure'].index(t)]
                           for t in common]) for key,d in x.items()}
    a,b,c=(aligned[k] for k in ('outer16','outer32','outer64'))
    gap=a-c
    diff=b-c
    views={}
    for k,h in [('16_minus_64',gap),('32_minus_64',diff)]:
        magnitude=float(np.sum(abs(h[-1])))
        views[k]=dict(terminal_gap_by_mode=h[-1].tolist(),
                      terminal_l1_gap=magnitude,
                      middle_two_mode_fraction_of_l1_gap=float(sum(abs(h[-1,3:5]))/magnitude)
                      if magnitude else None,
                      maximum_sampled_mode_l1_gap=float(np.max(np.sum(abs(h),axis=1))))
    return dict(common_exposure=common,terminal_mode_effects={k:v[-1].tolist() for k,v in aligned.items()},
                terminal_negative_modes={k:np.flatnonzero(v[-1]<0).tolist() for k,v in aligned.items()},
                first_to_last_mode_indices='zero-based',cross_width=views,
                terminal_cancellation={k:d['cancellation_fraction'][-1] for k,d in x.items()},
                aligned_mean_effects={k:v.tolist() for k,v in aligned.items()},
                scope='Saved mode count differences only, same original cold mesh and fixed inner width.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('outer16','outer32','outer64'):
        p.add_argument('--'+k,nargs=2,required=True,help='Native then local prepared output.')
    q=p.parse_args()
    print(json.dumps(compare({k:[json.loads(Path(f).read_text()) for f in y]
                              for k,y in vars(q).items()}),indent=2,allow_nan=False))
