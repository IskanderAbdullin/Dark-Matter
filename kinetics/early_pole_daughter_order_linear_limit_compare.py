"""Relate fixed-orbit quadrature and full-slow held order contrasts."""
import argparse
import json
from pathlib import Path

import numpy as np


def compare(linear, full):
    spec = np.asarray(linear['cold_source_by_mode'])
    if spec.shape != (2,65,len(linear['fixed_seed'])) or not np.all(np.isfinite(spec)):
        raise ValueError('Incomplete fixed-orbit mode histories.')
    histories = {}
    for n,k in (('16',4),('32',2)):
        p = full['results'][n]['individual_outer_comparison']
        if p['exposure'][-1]/64 != linear['leaf_width']:
            raise ValueError('Different fixed elastic leaf widths.')
        r = linear['cases'][n]
        # For units of cold counts, quadrature points are multiplied by
        # a full macro exposure before cumulative summation.
        d = spec[1]-spec[0]
        end, mid = d[::k], d[k//2::k]
        d0 = ((end[:-1]+end[1:])/2).sum(axis=1)
        d1 = mid.sum(axis=1)
        phase = p['exposure'][-1]/int(n)
        a,b = (np.asarray(p['cases'][j]['paired_count_history'])
               for j in ('elastic_centered','slow_centered'))
        # Downsample to the same 17 boundaries as the full comparison.
        smooth = np.r_[0.,np.cumsum(phase*(d1-d0))][::int(n)//16]
        measured = b-a
        t = p['exposure'][-1]
        if smooth[-1]==0. or measured[-1]==0.:
            raise ValueError('Unresolved trajectory separation.')
        histories[n] = dict(increment_points=p['exposure'],
                            linear_alternate_minus_original_count=smooth.tolist(),
                            full_alternate_minus_original_count=measured.tolist(),
                            linear_mean_order_gap=r['alternate_midpoint_mean']-r['original_trapezoid_mean'],
                            full_mean_order_gap=float(measured[-1]/t),
                            full_over_linear_end_gap=float(measured[-1]/smooth[-1]),
                            max_abs_count_gap_difference=float(max(abs(smooth-measured))))
    return dict(scope='Tangent at zero slow coupling versus full slow field, not an evaluation of their difference at small positive coupling.',
                paired_histories=histories)


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('linear')
    p.add_argument('full')
    q=p.parse_args()
    r=lambda f:json.loads(Path(f).read_text())
    print(json.dumps(compare(r(q.linear),r(q.full)),indent=2,allow_nan=False))
