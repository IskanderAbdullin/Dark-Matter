"""Compare the finite first-window result at hard aperture depths one, two and four."""
import argparse
import gzip
import json
from pathlib import Path

import numpy as np
from early_pole_joined_temporal_refinement_reference import conservation


def compare(baseline, refined, reference):
    controls = ('grid', 'projection', 'duration', 'steps', 'aperture_grid',
                'pair', 'method', 'screening_update', 'response_smear',
                'strength', 'same_soft_strength', 'compton_strength',
                'conversion_strength', 'elastic_energy_budget')
    a = baseline['snapshots'][0]
    before = baseline['hard_edge_sets']
    info = {}
    for depth, doc in zip((1,2,4), (baseline,*refined), strict=True):
        if any(doc[k] != baseline[k] for k in controls):
            raise ValueError(f'Mismatched controls at depth {depth}.')
        if doc.get('aperture_subdivision',1) != depth:
            raise ValueError(f'Mismatched aperture label at depth {depth}.')
        b = doc['snapshots'][0]
        np.testing.assert_array_equal(a['cold_bins'], b['cold_bins'])
        np.testing.assert_allclose(a['charged_and_photon_counts'],
                                   b['charged_and_photon_counts'], atol=2e-14, rtol=0)
        expected = reference['initial_only_depth'][str(depth)]
        np.testing.assert_allclose(expected['by_mode'], b['instantaneous_source_by_q'],
                                   atol=2e-14, rtol=1e-12)
        subdivided = doc['hard_edge_sets']
        for x,y in zip(before,subdivided,strict=True):
            if x[0]!=y[0] or x[-1]!=y[-1] or not set(x)<=set(y):
                raise ValueError(f'Hard endpoints or nested anchors changed at depth {depth}.')
        before = subdivided
        info[depth] = dict(mean=doc['scalar_gain']/doc['duration'],
                           relative_to_depth1_ppm=1e6*(doc['scalar_gain']/baseline['scalar_gain']-1),
                           initial=sum(b['instantaneous_source_by_q']),
                           terminal=sum(doc['snapshots'][-1]['instantaneous_source_by_q']),
                           hard_cells=[len(e)-1 for e in subdivided],
                           closure=conservation(doc))
    return info


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--native',nargs=3,required=True)
    p.add_argument('--continuous',nargs=3,required=True)
    p.add_argument('--static',required=True)
    a = p.parse_args()
    with gzip.open(a.static,'rt') as f:
        references=json.load(f)
    out={}
    for k,paths in (('steps',a.native),('cell-conservative',a.continuous)):
        docs=[json.loads(Path(p).read_text()) for p in paths]
        out[k]=compare(docs[0],docs[1:],references[k]['analysis'])
    contrast={d:100*(1-out['cell-conservative'][d]['mean']/out['steps'][d]['mean'])
              for d in (1,2,4)}
    print(json.dumps(dict(cases=out,contrast_percent=contrast,
                          scope='Finite hard aperture depths at held elastic quadrature and macro count.'),
                     indent=2,allow_nan=False))
