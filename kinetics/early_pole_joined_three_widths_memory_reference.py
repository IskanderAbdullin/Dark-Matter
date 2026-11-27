"""Audit two hard inputs on a common early window with three exact nested widths."""
import argparse
import gzip
import json
from pathlib import Path

import numpy as np

from early_pole_joined_history_memory_reference import compare as memory
from early_pole_joined_initial_window_refinement_reference import scores
from early_pole_joined_temporal_refinement_reference import select


def run(coarse,middle,native,continuous):
    out={}
    paths={}
    for p,f in [('steps',native),('cell-conservative',continuous)]:
        states=[select(path,p) for path in (coarse,middle,f)]
        out[p]=scores(states)
        paths[p]=states[-1]
    if any(x['resolution_divisors']!=[1,2,4] for x in out.values()):
        raise ValueError('Require successive factors of two in the matched time widths.')
    a,b=paths.values()
    m=memory(a,b)
    gap=np.asarray(a['state_history'])-np.asarray(b['state_history'])
    assert gap[0]==0.
    width=a['duration']/a['steps']
    later=np.flatnonzero(np.diff(gap)<0.)+1
    at=int(np.argmax(gap))
    m.update(cumulative_count_gap_at_first_coarse_window=float(gap[4]),
             maximum_signed_count_gap=float(gap[at]),exposure_at_maximum=at*width,
             first_interval_with_continuous_rate_above_native=(
                 float(later[0]*width) if len(later) else None))
    return dict(scope='Conditional fixed-scale early path and prepared-input memory.',
                temporal=out,memory=m,branches=paths)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for a in ('coarse','middle','native','continuous','output'):
        p.add_argument(a)
    a=p.parse_args()
    d=run(a.coarse,a.middle,a.native,a.continuous)
    b=(json.dumps(d,sort_keys=True,separators=(',',':'),allow_nan=False)+'\n').encode()
    path=Path(a.output)
    path.write_bytes(gzip.compress(b,mtime=0))
    assert json.loads(gzip.decompress(path.read_bytes()))==d
