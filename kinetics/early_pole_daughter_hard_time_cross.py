"""Mixed hard-depth and hold-window effect in the matched 32-mode model."""
import argparse
import json
from pathlib import Path

import numpy as np


def cross(h1,h4):
    if h1['exposure']!=h4['exposure']:
        raise ValueError("Unmatched time sampling.")
    t=np.array(h1['exposure'])
    e=t[-1]/2
    i=len(t)//2
    if t[i]!=e:
        raise ValueError("The hold windows are not contiguous equal halves.")
    d1,d4=(np.asarray(r['paired_count']) for r in (h1,h4))
    hard=d4-d1
    paired=[]
    for r in (h1,h4):
        z=np.array(r['paired_mode_sum_by_interval'])
        if z.shape!=(3,32) or np.max(abs(z[0]+z[1]-z[2]))>1e-18:
            raise ValueError("Missing matched mode ledgers.")
        paired.append(z)
    g=paired[1]-paired[0]
    if np.max(abs(g.sum(axis=1)-np.array([hard[i],hard[-1]-hard[i],hard[-1]])))>3e-17:
        raise ValueError("Mode contrast did not reproduce integrated totals.")
    panel=np.diff(hard)/np.diff(t)
    if np.max(abs(np.r_[hard[0],hard[0]+np.cumsum(panel*np.diff(t))]-hard))>3e-20:
        raise ValueError("Per-panel hard contrasts do not telescope.")
    means=dict(first=float(hard[i]/e),second=float((hard[-1]-hard[i])/e),whole=float(hard[-1]/(2*e)))
    return dict(scope="Matched two-window 32-mode hard-depth contrast; not FRW.",
                hard_shift_in_paired_mean=means,mixed_correction=float(means['second']-means['first']),
                hard_shift_paired_counts_at_times=hard.tolist(),exposure=t.tolist(),
                hard_shift_panel_means=panel.tolist(),
                hard_shift_positive_panels_by_window=[int(np.count_nonzero(v>0)) for v in (panel[:i],panel[i:])],
                hard_shift_second_window_panel_magnitudes_decrease=bool(np.all(np.diff(abs(panel[i:]))<0)),
                hard_shift_count_modes_positive=[int(np.count_nonzero(z>0)) for z in g],
                hard_shift_count_modes_negative=[int(np.count_nonzero(z<0)) for z in g],
                hard1_receipt=h1,hard4_receipt=h4)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('hard1')
    p.add_argument('hard4')
    a=p.parse_args()
    print(json.dumps(cross(*(json.loads(Path(k).read_text()) for k in (a.hard1,a.hard4))),
                     indent=2,allow_nan=False))
