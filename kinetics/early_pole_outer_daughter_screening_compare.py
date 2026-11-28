"""Separate two elastic topologies from their active response refresh."""
import argparse
import json
from pathlib import Path

from early_pole_joined_temporal_refinement_reference import conservation


def compare(active,held):
    keys=('grid','projection','duration','steps','pair','method','strength',
          'same_soft_strength','conversion_strength','compton_strength',
          'response_smear','elastic_energy_budget','hard_edge_sets','aperture_grid')
    data={}
    for topology in ('unlike','both'):
        a,b=active[topology],held[topology]
        for x,y in zip(a,b,strict=True):
            if any(x[k]!=y[k] for k in keys):
                raise ValueError(f'Mismatched controls for {topology}.')
            if x['screening_update']!='full' or y['screening_update']!='frozen':
                raise ValueError(f'Mismatched active/held response for {topology}.')
            if any(x['snapshots'][0][k]!=y['snapshots'][0][k] for k in
                   ('hard_moments','cold_bins','instantaneous_source_by_q')):
                raise ValueError(f'Mismatched initial conditions for {topology}.')
        case={}
        for profile,(native,outer) in (('active',a),('held',b)):
            if any(native[k]!=outer[k] for k in keys):
                raise ValueError(f'Mismatched preparations for {topology}, {profile}.')
            n,o=(x['scalar_gain']/x['duration'] for x in (native,outer))
            case[profile]=dict(native_mean=n,outer_mean=o,outer_effect=o-n,
                               closure=[conservation(q) for q in (native,outer)])
        case['active_minus_held_outer_effect']=(case['active']['outer_effect']
                                                 -case['held']['outer_effect'])
        data[topology]=case
    return dict(cases=data,scope='Deterministic conditional finite source sensitivity to response refresh.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for profile in ('active','held'):
        p.add_argument('--'+profile,nargs=4,required=True,
                       help='Native unlike, outside unlike, native both, outside both, in order.')
    q=p.parse_args()
    def read(paths):
        v=[json.loads(Path(p).read_text()) for p in paths]
        return dict(unlike=v[:2],both=v[2:])
    print(json.dumps(compare(read(q.active),read(q.held)),indent=2,allow_nan=False))
