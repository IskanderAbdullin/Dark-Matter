"""Compare signed paired window means at independently replayed hard depths."""
import argparse
import gzip
import json
from pathlib import Path

import numpy as np

from early_pole_daughter_two_window_compare import check


def read(k):
    if str(k).endswith('.gz'):
        with gzip.open(k,'rt') as o:
            return json.load(o)
    return json.loads(Path(k).read_text())


def compare(rows, depths=(1,2,4,8)):
    if len(rows)!=len(depths) or any(a['exposure']!=rows[0]['exposure'] for a in rows):
        raise ValueError('Matched report grid and each declared depth are required')
    t=np.array(rows[0]['exposure'])
    if t.shape!=(33,) or t[16]*2!=t[-1] or not np.all(np.diff(t)>0):
        raise ArithmeticError('Expected two contiguous original windows.')
    counts=np.array([a['paired_count'] for a in rows])
    modes=np.array([a['paired_mode_sum_by_interval'] for a in rows])
    means=np.array([[a['mean_shifts'][x] for x in ('first','second')] for a in rows])
    np.testing.assert_allclose(np.c_[counts[:,16], counts[:,-1]-counts[:,16]]/t[16], means,
                               rtol=4e-14,atol=2e-16)
    np.testing.assert_allclose(modes[:,:2].sum(axis=-1)/t[16],means,rtol=4e-14,atol=2e-16)
    if any(max(a['shared_prefix_max_abs_error'])>2e-18 for a in rows):
        raise ArithmeticError('A measured first-window continuation failed the prefix check.')
    delta=np.diff(means,axis=0)
    mode=np.diff(modes[:,:2],axis=0)
    panel=np.diff(np.diff(counts,axis=1),axis=0)
    return dict(depths=list(depths),means=means.tolist(),adjacent_differences=delta.tolist(),
                positive_changes_in_modes=np.count_nonzero(mode>0,axis=-1).tolist(),
                positive_changes_in_panels=np.count_nonzero(panel.reshape(-1,2,16)>0,axis=-1).tolist(),
                largest_absolute_mode_change=np.max(np.abs(mode),axis=-1).tolist(),
                crossings=[x['negative_to_nonnegative_brackets'] for x in rows])


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('prior_receipts',nargs=3)
    p.add_argument('--short-8',nargs=2,required=True,help='Native and local daughter')
    p.add_argument('--long-8',nargs=2,required=True,help='Native and local daughter')
    o=p.parse_args()
    measured=check([read(x) for x in o.short_8],[read(x) for x in o.long_8],8)
    print(json.dumps(compare([read(x) for x in o.prior_receipts]+[measured]),indent=2,allow_nan=False))
