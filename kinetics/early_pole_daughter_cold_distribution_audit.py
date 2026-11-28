"""Distinguish short production increments from normalized populated cold spectra."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


def metrics(a,b,nodes):
    a,b,q=(np.asarray(x,dtype=float) for x in (a,b,nodes))
    if a.shape!=q.shape or b.shape!=q.shape or np.any(a<0.) or np.any(b<0.):
        raise ValueError('Shape distances require nonnegative amounts on one ordered support.')
    total=math.fsum(a)
    ab=math.fsum(b)
    if min(total,ab)<=0. or np.any(np.diff(q)<=0.):
        raise ValueError('Missing positive totals or ordered nodes.')
    d=b-a
    dt=math.fsum(d)
    change=(d-a*(dt/total))/ab
    mean=math.fsum(a*q)/total
    shift=math.fsum((q-mean)*change)
    np.testing.assert_allclose(math.fsum(change),0.,rtol=0.,atol=1e-16)
    return dict(native_count=total,local_count=ab,paired_count_excess=dt,
                relative_count_change=dt/total,
                normalized_total_variation=0.5*math.fsum(abs(change)),
                normalized_q_transport_distance=math.fsum(abs(np.cumsum(change)[:-1])*np.diff(q)),
                native_mean_q_over_T=mean,
                relative_normalized_mean_q_shift=shift/mean,
                count_delta_by_mode=d.tolist())


def compare(a,b):
    for k in ('duration','projection','report_panels','rtol','base_atol','history_exposure','cold_momentum_nodes_over_T'):
        if a[k]!=b[k]:
            raise ValueError('Mismatched chart and mode controls.')
    if a.get('prepared_charged_swap') or b.get('prepared_charged_swap')!='local-daughter':
        raise ValueError('Expected native and local paths.')
    for q in (a,b):
        if not (abs(q['energy_defect'])<3e-12 and abs(q['entropy_defect'])<3e-11
                and abs(q['hard_number_defect'])<3e-12):
            raise ValueError('Unexpected original ledger failure.')
    x,y=(np.asarray(z['cold_mode_count_history']) for z in (a,b))
    np.testing.assert_array_equal(x[0],y[0])
    q=a['cold_momentum_nodes_over_T']
    n=len(q)
    if n not in (8,32) or x.shape!=y.shape or x.shape!=(a['report_panels']+1,n):
        raise ValueError('Unexpected recorded mode grid.')
    v=metrics(x[-1]-x[0],y[-1]-y[0],q)
    u=metrics(x[-1],y[-1],q)
    d=np.asarray(v['count_delta_by_mode'])
    bsize=n//8
    center=math.fsum(d[4*bsize:6*bsize])
    return dict(increment_shape=v,terminal_shape=u,
                share_of_count_excess_in_original_bins_5_and_6=center/v['paired_count_excess'])


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('cold8','cold32'):
        p.add_argument('--'+k,nargs=2,required=True,help='Native then local mode history.')
    q=p.parse_args()
    print(json.dumps(dict(cases={k:compare(*[json.loads(Path(f).read_text()) for f in v])
                                  for k,v in vars(q).items()},
                          scope='Within each matching held cold support only, not relic or later transfer.'),
                     indent=2,allow_nan=False))
