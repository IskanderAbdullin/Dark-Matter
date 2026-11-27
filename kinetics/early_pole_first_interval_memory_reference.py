"""Measure how rapidly the initial hard reconstruction differences decay in the cold source."""
import argparse
import json

import numpy as np

from early_pole_joined_temporal_refinement_reference import conservation,select


def compare(a,b):
    keys=('grid','aperture_grid','duration','steps','method','compton_strength',
          'conversion_strength','strength','same_soft_strength','pair','screening_update')
    if a['projection']!='steps' or b['projection']!='cell-conservative':
        raise ValueError('Expected ordered native and continuous hard inputs.')
    if any(a[k]!=b[k] for k in keys) or a.get('elastic_energy_budget')!=b.get('elastic_energy_budget'):
        raise ValueError('Physical controls and numerical budget must match.')
    ledgers=[conservation(r) for r in (a,b)]
    n=a['steps']
    h=a['duration']/n
    a0,b0=(np.asarray(r['state_history']) for r in (a,b))
    if not len(a0)==len(b0)==n+1 or a0[0]!=b0[0]:
        raise ValueError('Initial count or history length is not matched.')
    x,y=(np.diff(r)/h for r in (a0,b0))
    if not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)) or min(x)<=0.:
        raise ValueError('Use finite positive native averaged rates.')
    gap=(x-y)/x
    diff=a0-b0
    if diff[-1]==0.:
        raise ValueError('Cannot normalize the accumulated contrast by zero.')
    idx=sorted(set([1,n//16,n//8,n//4,n//2,3*n//4,n]))
    rows=[]
    for i in idx:
        if i<1:
            continue
        rows.append(dict(index=i,time=i*h,native_rate=float(x[i-1]),continuous_rate=float(y[i-1]),
                         relative_mean_rate_gap=float(gap[i-1]),relative_max_later_mean_rate_gap=float(np.max(abs(gap[i-1:]))),
                         cumulative_fraction_of_signed_gap=float(diff[i]/diff[-1])))
    def tail_above(eps):
        over=np.flatnonzero(abs(gap)>=eps)
        last=int(over[-1]+1) if len(over) else 0
        return None if last==n else last*h
    f={str(e):tail_above(e) for e in (.10,.05,.01,.005,.002)}
    q=np.asarray(a['snapshots'][-1]['instantaneous_source_by_q'])
    p=np.asarray(b['snapshots'][-1]['instantaneous_source_by_q'])
    return dict(steps=n,width=h,rows=rows,matched_conservation=ledgers,
                total_signed_count_gap=float(diff[-1]),all_remaining_rate_gaps_below_after=f,
                terminal_instantaneous_rate_difference=float((q.sum()-p.sum())/q.sum()),
                terminal_instantaneous_rate_gap_each_mode=(q-p).tolist(),
                scope='Memory within only the matched first coarse interval; averaged rates are not instantaneous slopes.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('native')
    p.add_argument('continuous')
    a=p.parse_args()
    print(json.dumps(compare(select(a.native,'steps'),select(a.continuous,'cell-conservative')),
                     indent=2,allow_nan=False))
