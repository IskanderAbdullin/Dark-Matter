"""Compare finite early-window rates with each elastic channel removed separately."""
import argparse
import json

import numpy as np

from early_pole_joined_temporal_refinement_reference import select,conservation


def compare(runs):
    base=runs['none']
    controls=('grid','projection','duration','steps','method','pair','compton_strength',
              'conversion_strength','hard_edge_sets','aperture_grid','response_smear')
    patterns=dict(none=(0.,0.,'frozen'),unlike=(1.,0.,'full'),same=(0.,1.,'frozen'),
                  both=(1.,1.,'full'))
    ledgers={}
    paths={}
    for key,r in runs.items():
        if any(base[k]!=r[k] for k in controls):
            raise ValueError('Unmatched controls in factorial comparison.')
        if (r['strength'],r['same_soft_strength'],r['screening_update'])!=patterns[key]:
            raise ValueError('Unexpected elastic toggle pattern.')
        for k in ('hard_moments','cold_bins','instantaneous_source_by_q'):
            np.testing.assert_array_equal(r['snapshots'][0][k],base['snapshots'][0][k])
        assert r['snapshots'][0]['exposure']==0.
        ledgers[key]=conservation(r)
        p=np.asarray(r['state_history'])
        assert p.shape==(r['steps']+1,) and np.isfinite(p).all() and p[-1]-p[0]==r['scalar_gain']
        paths[key]=p-p[0]
    nominal={k:float(p[-1]/base['duration']) for k,p in paths.items()}
    interaction=paths['both']-paths['same']-paths['unlike']+paths['none']
    return dict(projection=base['projection'],duration=base['duration'],steps=base['steps'],
                interval_mean_rates=nominal,
                signed_factorial_interaction_in_mean_rate=float(interaction[-1]/base['duration']),
                terminal_source_by_toggle={k:float(sum(r['snapshots'][-1]['instantaneous_source_by_q']))
                                           for k,r in runs.items()},
                histories={k:p.tolist() for k,p in paths.items()},ledgers=ledgers,
                scope='Matched finite early window; separating interactions conditional on this split and not a decomposition of the continuous kinetics.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--projection',choices=('steps','cell-conservative'),default='steps')
    for a in ('none','unlike','same','both'):
        p.add_argument('--'+a,required=True)
    args=vars(p.parse_args())
    projection=args.pop('projection')
    print(json.dumps(compare({k:select(v,projection) for k,v in args.items()}),indent=2,allow_nan=False))
