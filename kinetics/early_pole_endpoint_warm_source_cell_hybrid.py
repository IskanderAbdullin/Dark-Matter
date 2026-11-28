"""Train a frozen warm-cell mask at one snapshot and check another."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np


def payload(d):
    cases={r['rule']:np.array(r['warm_cell_forward_inverse_net']) for r in d['cases']}
    edges=d['warm_native_edges']
    if len(edges)<17 or not np.all(np.diff(edges)>0) or not d['mesh_reference']['global_L1_over_total']:
        raise ArithmeticError('Invalid paired reference support or envelope')
    mats=[cases[k] for k in ('W1','Gauss2','Volume64')]
    if any(a.shape!=(len(edges)-1,3) or not np.all(np.isfinite(a)) for a in mats):
        raise ArithmeticError('Invalid warm-cell data')
    return [a[:,2] for a in mats]


def run(training,holdout):
    if (training['x'],holdout['x'])!=(.01,.05) or training['warm_native_edges']!=holdout['warm_native_edges']:
        raise ArithmeticError('Need ordered early/held-out endpoints on the exact same native edges')
    edges=training['warm_native_edges']
    p,g,r=payload(training)
    err=abs(p-r)
    rank=np.argsort(-err,kind='stable')
    n=len(rank)
    out=[]
    for count in (0,1,3,5,8,12,n):
        select=rank[:count]
        row=dict(refined_native_cells=count,selected_q=[edges[i:i+2] for i in select])
        for d,name in ((training,'early'),(holdout,'held_out')):
            p,g,r=payload(d)
            val=p.copy()
            val[select]=g[select]
            e=val-r
            scale=fsum(r)
            if abs(scale/d['ref_net_pressure']-1)>2e-14 or scale<=0:
                raise ArithmeticError('Mismatched scalar normalization')
            row[name]=dict(relative_signed_error=fsum(e)/scale,relative_L1_residual=fsum(abs(e))/scale)
        out.append(row)
    return dict(train_x=training['x'],test_x=holdout['x'],total_warm_cells=n,rows=out,
                scope='Mix already-resolved frozen affine sums by cell, not an evolved adaptive chart or an error bound on individual cells.')


if __name__=='__main__':
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument('early',type=Path)
    a.add_argument('holdout',type=Path)
    z=a.parse_args()
    print(json.dumps(run(json.loads(z.early.read_text()),json.loads(z.holdout.read_text())),indent=2,allow_nan=False))
