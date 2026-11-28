"""Compare original, refreshed-cut, and startup-cut scalar source projections."""
import argparse
import json
from pathlib import Path

import numpy as np


def compare(m8,m16,f8,f16):
    runs=(m8,m16,f8,f16)
    if [r['quadrature_order'] for r in runs]!=[8,16,8,16]:
        raise ValueError('Expected both cut controls at the same two checked Gauss orders.')
    t=np.asarray(m16['exposure'])
    if any(r['exposure']!=t.tolist() for r in runs) or f8['frozen_angular_knots']!=f16['frozen_angular_knots']:
        raise ValueError('Mismatched time samples or held startup cuts.')
    tags=('native','local-daughter')
    data={}
    for r in runs:
        for k in tags:
            if len(r['arms'][k])!=len(t):
                raise ValueError('Unexpected stored trajectory size.')
    for k in tags:
        native=m16['arms'][k]
        roots=np.asarray([x['angular_modes'] for x in native])
        if roots.shape!=(len(t),7):
            raise ValueError('Seven consistent sorted angular crossings are expected.')
        rows={}
        old=np.asarray([x['original_mode_rates'] for x in native])
        for j,r in zip(('move8','move16','fixed8','fixed16'),runs,strict=True):
            arms=r['arms'][k]
            np.testing.assert_array_equal([x['original_mode_rates'] for x in arms],old)
            np.testing.assert_array_equal([x['angular_modes'] for x in arms],roots)
            values=np.asarray([x['mode_grouped_rate'] for x in arms])
            if values.shape!=old.shape or not np.all(np.isfinite(values)):
                raise ValueError('Wrong scalar aperture quadrature output.')
            rows[j]=values
        rows['old']=old
        data[k]=dict(
            max_abs_sorted_root_motion=float(np.max(abs(roots-roots[0]))),
            max_abs_moving8_to_16_mode_rate=float(np.max(abs(rows['move8']-rows['move16']))),
            max_abs_fixed8_to_16_mode_rate=float(np.max(abs(rows['fixed8']-rows['fixed16']))),
            max_abs_fixed_to_moving_at16_mode_rate=float(np.max(abs(rows['fixed16']-rows['move16']))),
            rates=rows)
    n,l=(data[k]['rates'] for k in tags)
    pair={k:(l[k]-n[k]).sum(axis=-1) for k in n}
    brackets=lambda a:[t[j:j+2].tolist() for j in np.flatnonzero((a[:-1]<0.)&(a[1:]>=0.))]
    return dict(scope='Matched already-evolved snapshots; only evaluate different integrations of the instantaneous source. No time extrapolation.',
                exposure=t.tolist(),arms={k:{j:v for j,v in d.items() if j!='rates'} for k,d in data.items()},
                updated_reference_trajectory_rates_by_mode={k:d['rates']['move16'].tolist() for k,d in data.items()},
                paired_rate_series={k:v.tolist() for k,v in pair.items()},
                negative_to_nonnegative_brackets={k:brackets(v) for k,v in pair.items()},
                max_absolute_paired_moving_minus_original=float(np.max(abs(pair['move16']-pair['old']))),
                max_absolute_paired_fixed_minus_moving=float(np.max(abs(pair['fixed16']-pair['move16']))),
                original_sign_changes_vs_moving=int(np.count_nonzero(np.sign(pair['old'])!=np.sign(pair['move16']))),
                moving_sign_changes_vs_fixed=int(np.count_nonzero(np.sign(pair['move16'])!=np.sign(pair['fixed16']))))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('inputs',nargs=4,help='Moving 8 and 16, then frozen-startup 8 and 16.')
    a=p.parse_args()
    print(json.dumps(compare(*(json.loads(Path(f).read_text()) for f in a.inputs)),indent=2,allow_nan=False))
