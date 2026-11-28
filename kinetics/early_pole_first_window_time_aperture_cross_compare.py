"""Separate macro time sensitivity from hard aperture depth in a crossed control."""
import argparse
import json
from pathlib import Path

import numpy as np
from early_pole_joined_temporal_refinement_reference import conservation


def compare(inputs):
    out={}
    for projection, raw in inputs.items():
        base=raw[(16,1)]
        contrasts={}
        for (n,d),doc in raw.items():
            overlap=('grid','projection','duration','aperture_grid','pair','method',
                     'strength','same_soft_strength','compton_strength','conversion_strength',
                     'screening_update','response_smear','elastic_energy_budget')
            if any(base[k]!=doc[k] for k in overlap):
                raise ValueError(f'Mismatched controls for {projection}, {(n,d)}.')
            if doc['steps']!=n or doc.get('aperture_subdivision',1)!=d:
                raise ValueError(f'Mismatched labels for {projection}, {(n,d)}.')
            for k in ('hard_moments','cold_bins','instantaneous_source_by_q'):
                np.testing.assert_array_equal(raw[(16,d)]['snapshots'][0][k],doc['snapshots'][0][k])
            if doc['hard_edge_sets']!=raw[(16,d)]['hard_edge_sets']:
                raise ValueError(f'Mismatched hard support across the time factor at {d}.')
            xa,xb=(x['snapshots'][0] for x in (base,doc))
            np.testing.assert_array_equal(xa['cold_bins'],xb['cold_bins'])
            np.testing.assert_allclose(xa['charged_and_photon_counts'],xb['charged_and_photon_counts'],
                                       atol=2e-14,rtol=0)
            contrasts[(n,d)]=doc['scalar_gain']/doc['duration']
        a,b,c,d=(contrasts[k] for k in ((16,1),(16,4),(32,1),(32,4)))
        out[projection]=dict(
            levels={f'{n}x{h}':dict(mean=value,closure=conservation(raw[(n,h)]))
                    for (n,h),value in contrasts.items()},
            aperture_effect=[b-a,d-c],time_effect=[c-a,d-b],
            interaction=(d-c)-(b-a))
    paired={k:100*(1-out['cell-conservative']['levels'][k]['mean']/v['mean'])
            for k,v in out['steps']['levels'].items()}
    return dict(cases=out,projection_contrast_percent=paired,
                scope='A time by hard-aperture finite sensitivity at one cold mesh and elastic quadrature.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('native','continuous'):
        p.add_argument('--'+k,nargs=4,required=True,help='16x1, 16x4, 32x1, 32x4, in order.')
    q=p.parse_args()
    def read(paths):
        return {k:json.loads(Path(f).read_text()) for k,f in
                zip(((16,1),(16,4),(32,1),(32,4)),paths,strict=True)}
    print(json.dumps(compare({'steps':read(q.native),'cell-conservative':read(q.continuous)}),
                     indent=2,allow_nan=False))
