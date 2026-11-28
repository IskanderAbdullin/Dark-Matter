"""Derive a sharp conditional L1 ranking from the two actual terminal receipts."""
import argparse
import hashlib
import json
from math import fsum
from pathlib import Path

import numpy as np


def loss(box):
    a,b=np.asarray(box).T
    if np.any(b<a) or not np.all(np.isfinite([a,b])):
        raise ArithmeticError('Invalid finite-node affine residual interval')
    # Every cell interval can be reached independently as a convex
    # interpolation of its two fixed-count/mean extreme occupation profiles.
    lower=np.where(a*b<=0.,0.,np.minimum(abs(a),abs(b)))
    return [fsum(lower),fsum(np.maximum(abs(a),abs(b)))]


def run(reconstruction,hull):
    rr=reconstruction.read_bytes()
    vv=hull.read_bytes()
    r,h=json.loads(rr),json.loads(vv)
    if ((r['x'],r['native_sha256'],r['fixed_optics_mass'],r['fixed_clock_H'])!=
            (h['x'],h['native_sha256'],h['common_m'],h['common_clock_H'])):
        raise ArithmeticError('Compare only the same verified actual endpoints and source clock')
    rows=[]
    for a,b in zip(r['paired_actual_results'],h['paired_chart_results'],strict=True):
        native=float(a['common_reference_scale'])
        d=np.asarray(a['all_cell_errors'][0])
        step=np.asarray(a['all_cell_errors'][3])
        S=np.asarray([x['allowed_net'] for x in b['local_rows']])
        v=np.asarray([x['actual_net'] for x in b['local_rows']])
        q=[x['q'][0] for x in b['local_rows']]+[b['local_rows'][-1]['q'][1]]
        if (a['input_sha256']!=b['input_sha256'] or q!=r['scalar_cells']
                or len(d)!=len(v) or len(step)!=len(v) or native<=0.):
            raise ArithmeticError('Changed source states or native warm cells')
        if abs(fsum(v)/b['net']-1)>3e-13 or abs(fsum(d)-(b['net']/native-1.))>3e-13:
            raise ArithmeticError('Independent aggregate source ledgers do not close')
        band=(S-v[:,None])/native
        out={}
        for key,e,scale in (('volume96',d,1.),('volume64',d-step,1+fsum(step))):
            if not scale>0.:
                raise ArithmeticError('Nonpositive alternate reference')
            box=(band+e[:,None])/scale
            window=loss(box)
            measured=fsum(abs(e))/scale
            if not window[0]-2e-13<=measured<=window[1]+2e-13:
                raise ArithmeticError('Measured residual not inside independently certified loss window')
            out[key]=dict(allowed_relative_local_net_L1=window,actual_relative_local_net_L1=measured)
        rows.append(dict(chart=a['chart'],source_sha256=a['input_sha256'],denominator_Volume96=native,checks=out))
    G,M=rows
    ans={}
    for key in G['checks']:
        a=G['checks'][key]['allowed_relative_local_net_L1'][0]
        b=M['checks'][key]['allowed_relative_local_net_L1'][1]
        gap=a-b
        if gap<=2e-12:
            raise ArithmeticError('No ranking independent of allowed node moments')
        ans[key]=dict(minimum_gap_G7_minus_mixed255=gap,
                      equal_independent_normalized_L1_change_per_reference_strictly_below_to_preserve_gap=gap/2)
    return dict(x=.05,source_sha256=[hashlib.sha256(rr).hexdigest(),hashlib.sha256(vv).hexdigest()],
                chart_results=rows,ranking=ans,
                scope='Sharp L1 intervals on fixed evolved nodes under each stored count/mean and a specified Volume64 or Volume96 reference. The ranking need not hold for arbitrary continuum reconstructions or time histories.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('reconstruction',type=Path)
    p.add_argument('hull',type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.reconstruction,a.hull),indent=2,allow_nan=False))
