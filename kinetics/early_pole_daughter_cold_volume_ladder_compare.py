"""Compare volume-centroid sampling with independent cold-volume quadrature."""
import argparse
import json
from pathlib import Path

import numpy as np


def compare(ladder,plain,split):
    out={}
    for depth,rows in ladder['hard_depths'].items():
        p=[next(x for x in v if str(x['hard_depth'])==depth) for v in plain]
        g=[next(x for x in v if str(x['hard_depth'])==depth) for v in split]
        r=g[-1]
        if any(v['initial_pole']!=r['initial_pole'] or v['grouped_seed_error']>5e-18 for v in p+g):
            raise ValueError('Comparison changed the scalar initial data.')
        if any(v.get('angular_edge_modes') for v in p):
            raise ValueError('The no-cut runs must not include extra geometric knots.')
        knots=np.asarray(r['angular_edge_modes'])
        for v in g:
            np.testing.assert_allclose(v['angular_edge_modes'],knots,atol=1e-13,rtol=0.)
        if len(knots)==0:
            raise ValueError('No genuine geometric knots discovered.')
        val=r['direct_scalar_volume_rate']
        if val>=0.:
            raise ValueError('Expected resolved frozen suppression.')
        split_group=[v['grouped_rate'] for v in g]
        expected=np.asarray(r['grouped_rate'])
        off=np.asarray([v['direct_pair'] for v in rows])-val
        out[depth]=dict(reference_paired_rate=val,extra_scalar_angular_knots=knots.tolist(),
                       cold_ladder_rates=rows,lumped_sample_errors_vs_reference=off.tolist(),
                       lumped_relative_errors_abs=(abs(off)/abs(val)).tolist(),
                       absolute_peak_quadrature_group_delta_with_geometric_knots=float(np.max(abs(np.asarray(split_group)-expected))),
                       absolute_peak_quadrature_total_delta_with_geometric_knots=max(
                           abs(x['direct_scalar_volume_rate']-val) for x in g),
                       plain_quadratures=p,geometric_quadratures=g)
    return dict(scope='Within-window frozen cold source and archived volume density, no revised cosmological evolution.',
                hard_depths=out)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('ladder')
    p.add_argument('--plain',nargs='+',required=True)
    p.add_argument('--split',nargs='+',required=True)
    q=p.parse_args()
    r=lambda f:json.loads(Path(f).read_text())
    print(json.dumps(compare(r(q.ladder),[r(v) for v in q.plain],[r(v) for v in q.split]),indent=2,allow_nan=False))
