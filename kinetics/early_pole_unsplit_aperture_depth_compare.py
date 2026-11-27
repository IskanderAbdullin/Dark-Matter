"""Compare held-response simultaneous integrations at two hard aperture depths."""
import argparse
import gzip
import json
from pathlib import Path

import numpy as np


def compare(original, refined, depth4_bundle):
    receipt={}
    controls=('duration','projection','report_panels','rtol','base_atol','history_exposure')
    projections=('steps','cell-conservative')
    for p,a,b in zip(projections,original,refined,strict=True):
        z=depth4_bundle[p]
        if any(a[k]!=b[k] for k in controls) or a['projection']!=p:
            raise ValueError('Mismatched chart controls at fixed preparation.')
        if a.get('aperture_subdivision',1)!=1 or b.get('aperture_subdivision')!=4:
            raise ValueError('Expected exact aperture labels 1 and 4.')
        if z['projection']!=p or z['aperture_subdivision']!=4 or z['duration']!=b['duration']:
            raise ValueError('Wrong prepared support archive.')
        for x in (a,b):
            if (x.get('prepared_charged_swap') or x.get('elastic_channels',True) is not True
                    or x.get('selected_elastic_channel','both')!='both'):
                raise ValueError('Only the full held chart at the original projection is supported.')
            if not (abs(x['energy_defect'])<3e-12 and abs(x['entropy_defect'])<3e-11
                    and abs(x['hard_number_defect'])<3e-12):
                raise ValueError('A simultaneous output fails original gates.')
        h0,h1=(np.asarray(x['cold_count_history']) for x in (a,b))
        if h0.shape!=h1.shape or h0.shape!=(a['report_panels']+1,) or h0[0]!=h1[0]:
            raise ValueError('Mismatched sample grid or initial cold count.')
        np.testing.assert_allclose(h1[0],sum(z['snapshots'][0]['cold_bins']),rtol=0.,atol=1e-16)
        rates=np.array([(h[-1]-h[0])/b['duration'] for h in (h0,h1)])
        np.testing.assert_allclose(rates,[a['mean_cold_rate'],b['mean_cold_rate']],rtol=2e-15,atol=0.)
        receipt[p]=dict(means=rates.tolist(),mean_depth_shift=float(rates[1]-rates[0]),
                        relative_depth_change_ppm=float(1e6*(rates[1]/rates[0]-1)),
                        max_absolute_sampled_cold_count_change=float(np.max(abs(h1-h0))),
                        cumulative_depth_effect=(h1-h0).tolist(),
                        energy_defects=[a['energy_defect'],b['energy_defect']],
                        entropy_defects=[a['entropy_defect'],b['entropy_defect']],
                        calls=[a['calls'],b['calls']])
    n,c=(np.asarray(receipt[p]['means']) for p in projections)
    d=100*(1-c/n)
    return dict(cases=receipt,contrast_percent=d.tolist(),
                contrast_change_percentage_points=float(d[1]-d[0]),
                scope='Finite held-screening hard support control, not an active split comparison or a convergence bound.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--original',nargs=2,required=True,help='Original native, then continuous.')
    p.add_argument('--refined',nargs=2,required=True,help='Depth4 native, then continuous.')
    p.add_argument('--depth4-bundle',required=True)
    q=p.parse_args()
    def read(paths):
        return [json.loads(Path(f).read_text()) for f in paths]
    with gzip.open(q.depth4_bundle,'rt') as f:
        bundle=json.load(f)['branches']
    print(json.dumps(compare(read(q.original),read(q.refined),bundle),indent=2,allow_nan=False))
