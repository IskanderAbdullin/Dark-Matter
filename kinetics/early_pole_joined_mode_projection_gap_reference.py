"""Separate coherent and cancelling source shifts in matched scalar snapshots."""
import argparse
import json
import numpy as np
from early_pole_joined_history_memory_reference import load

def compare(first,second):
    if first['projection'] == second['projection']:
        raise ValueError('Use two distinct initial projections.')
    for key in ('grid','steps','duration','aperture_grid','hard_edge_sets'):
        if first[key]!=second[key]:
            raise ValueError(f'Incompatible snapshot {key}.')
    rows=[]
    for a,b in zip(first['snapshots'],second['snapshots'],strict=True):
        if a['step']!=b['step'] or a['exposure']!=b['exposure']:
            raise ValueError('Exposure samples do not align.')
        x,y=(np.asarray(u['instantaneous_source_by_q']) for u in (a,b))
        if x.shape!=y.shape or x.ndim!=1:
            raise ValueError('Scalar modes do not align.')
        d=x-y
        net=float(d.sum())
        absolute=float(np.abs(d).sum())
        rows.append(dict(step=a['step'],exposure=a['exposure'],shift_by_q=d.tolist(),
                         negative=int(np.count_nonzero(d<0.)),positive=int(np.count_nonzero(d>0.)),
                         net=net,absolute=absolute,cancellation_fraction=1.-abs(net)/absolute if absolute else None))
    return dict(grid=first['grid'],projection_order=[first['projection'],second['projection']],rows=rows)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('first')
    p.add_argument('second')
    a=p.parse_args()
    print(json.dumps(compare(load(a.first),load(a.second)),indent=2,allow_nan=False))
