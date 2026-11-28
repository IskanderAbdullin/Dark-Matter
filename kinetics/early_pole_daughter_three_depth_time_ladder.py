"""Joint depth/time contrasts after the three independently verified paired windows."""
import argparse
import json
from pathlib import Path

import numpy as np


def calculate(a):
    if len(a)!=3 or any(v['exposure']!=a[0]['exposure'] for v in a):
        raise ValueError("Common exposure grid and depths 1, 2 and 4 are required.")
    t=np.asarray(a[0]['exposure'])
    n=len(t)//2
    if t[-1]!=2*t[n]:
        raise ValueError("The declared windows are not equal.")
    rates=np.array([[v['mean_shifts'][k] for k in ('first','second')] for v in a])
    counts=np.array([x['paired_count'] for x in a])
    if not np.allclose(rates,np.c_[counts[:,n],counts[:,-1]-counts[:,n]]/t[n],rtol=2e-14,atol=2e-16):
        raise ValueError("Per-window counts did not reproduce the rate.")
    delta=np.diff(rates,axis=0)
    modes=np.asarray([x['paired_mode_sum_by_interval'] for x in a])
    np.testing.assert_allclose(modes.sum(axis=-1)[:,:2],rates*t[n],atol=1e-20,rtol=5e-14)
    pm=np.diff(np.diff(counts,axis=1),axis=0)>0
    px=np.diff(modes,axis=0)>0
    ratio=modes[2,1,:]/modes[1,1,:]
    return dict(scope="Three matched hard-depth/held-time controls; no asymptotic convergence assumed.",
                depths=[1,2,4],window_means=rates.tolist(),
                adjacent_depth_corrections=delta.tolist(),
                adjacent_second_minus_first_depth_correction=np.diff(delta,axis=1).ravel().tolist(),
                magnitude_ratio_of_adjacent_corrections_by_window=[
                    float(abs(delta[1,j]/delta[0,j])) if delta[0,j] else None for j in range(2)],
                panels_with_positive_adjacent_change=[[int(x[:n].sum()),int(x[n:].sum())] for x in pm],
                modes_with_positive_adjacent_change=px.sum(axis=-1).tolist(),
                late_depth4_to_depth2_mode_ratio_range=[float(ratio.min()),float(ratio.max())],
                adjacent_change_in_paired_counts=np.diff(counts,axis=0).tolist(),
                negative_to_nonnegative_brackets=[v['negative_to_nonnegative_brackets'] for v in a])


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('receipts',nargs=3)
    v=p.parse_args()
    print(json.dumps(calculate([json.loads(Path(p).read_text()) for p in v.receipts]),
                     indent=2,allow_nan=False))
