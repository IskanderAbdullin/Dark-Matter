"""Cross-check independently found frozen aperture packets and integral values."""
import argparse
import json
from pathlib import Path

import numpy as np


def compare(default, low, high):
    out=[]
    for x,a,b in zip(default,low,high,strict=True):
        if (a['root_scan'],b['root_scan'])!=(129,257):
            raise ValueError('Compare the independent 129/257 root scans.')
        settings=(x['hard_depth'],x['cold_modes'])
        if any(settings!=(v['hard_depth'],v['cold_modes']) for v in (a,b)):
            raise ValueError('Different fixtures compared.')
        n=x['sample_counter']//12
        if x['sample_counter']%12 or any(v['independently_located_packet_count']!=n for v in (a,b)):
            raise ValueError('A route located a different admissible packet count.')
        reference=np.asarray(x['aperture_terms_by_mode_and_cell'])[:,:,0].sum(axis=1)
        y,z=(np.asarray(v['adaptive_direct_by_mode']) for v in (a,b))
        residual=float(np.max(abs(z-reference)))
        root=float(np.max(abs(z-y)))
        total=float(b['adaptive_direct_rate']-x['pointwise_archive_paired_total'])
        if max(residual,root,abs(total))>1e-11:
            raise ValueError('Independent routes do not agree within the stated check.')
        out.append(dict(hard_depth=settings[0],cold_modes=settings[1],packet_count=n,
                        adaptive_vs_default_total_difference=total,
                        adaptive_vs_default_max_mode_error=residual,
                        adaptive_129_vs_257_max_mode_error=root))
    return dict(scope='Independent verification of frozen affine dressed energy integration and packet location, not a continuum or physical error bound.',
                checked_fixtures=out,adaptive129=low,adaptive257=high)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('default','scan129','scan257'):
        p.add_argument(k)
    v=p.parse_args()
    q=lambda f:json.loads(Path(f).read_text())
    print(json.dumps(compare(*(q(f) for f in vars(v).values())),indent=2,allow_nan=False))
