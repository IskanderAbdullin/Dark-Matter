"""Attribute three-clock G9 count and pressure differences to original cells."""
import argparse
import hashlib
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_cell_endpoint_startup_pilot import drive_checkpoint_digest
from early_pole_endpoint_gauss_warm_comparison import chart_for, gauss
from early_pole_endpoint_scalar_subshell_chain_receipt import read
from early_pole_endpoint_scalar_subshell_pilot import ScalarSubshellChart,SPECS,load


def run(fine, coarse, ultrafine, third):
    (F,hf),(A,ha),(B,hb)=map(read,(fine,coarse,ultrafine))
    raw=third.read_bytes()
    T=json.loads(raw)
    if (T['steps']!=[1024,2048,4096] or (T['x_start'],T['x_end'])!=(.02,.04) or hf!=T['input_hashes'][1]
            or B['payload_sha256']!=T['checkpoint_payload_sha256'] or F['warm_gauss_order']!=9
            or [v['completed'] for v in (A,B)]!=[1024,4096]
            or any(v['payload_sha256']!=drive_checkpoint_digest(v) for v in (A,B))):
        raise ArithmeticError('Exact published third-clock evidence and original two checkpoints required')
    base=ScalarSubshellChart(load('fine')[0],SPECS['fine'],2,4,1e-5)
    _,d,repeats=chart_for(F)
    cards=[dict(F,steps=1024,nodal_occupation=A['occupation']),F,dict(F,steps=4096,nodal_occupation=B['occupation'])]
    Z=[gauss(q,base,T['mass'],.04,expected_steps=k) for q,k in zip(cards,T['steps'],strict=True)]
    if not np.array_equal(np.array([a for a,_,_ in Z]),T['terminal_moments']):
        raise ArithmeticError('Each full first two moments and pressure must reproduce exactly')
    cells=np.array([p for _,_,p in Z])
    result=[]
    for i in (0,1):
        den=Z[i+1][0][:2]
        q=(cells[i]-cells[i+1])/den
        sums=np.asarray([fsum(p) for p in q.T])
        absolute=np.asarray([fsum(abs(p)) for p in q.T])
        cold=repeats==1
        if np.max(abs(sums-T['previous_over_next_minus_one'][i][:2]))>2e-14:
            raise ArithmeticError('Cell differences do not sum to measured clock shift')
        if min(absolute)<=0.:
            raise ArithmeticError('Expected a nonzero resolved contrast in each observable')
        groups={}
        for label,keys in (('cold',cold),('refined_warm',~cold)):
            groups[label]=dict(signed=[fsum(p) for p in q[keys].T],
                               share_of_absolute=[fsum(abs(p))/a for p,a in zip(q[keys].T,absolute,strict=True)])
        ids=np.argsort(-abs(q),axis=0)[:4]
        cells_top=[[
            dict(q_bounds=d.ds.ed[k:k+2].tolist(),signed=q[k,j].item(),relative_absolute_share=(abs(q[k,j])/absolute[j]).item())
            for k in ids[:,j]] for j in (0,1)]
        result.append(dict(clocks=T['steps'][i:i+2],signed_relative_changes=sums.tolist(),absolute_l1_changes=absolute.tolist(),
                           groups=groups,largest_cells_by_field=cells_top))
    return dict(x=.04,fields=['number','pressure'],rows=result,
                input_sha256=[hf,ha,hb,hashlib.sha256(raw).hexdigest()],
                scope='Within one archived G9 chart and common mass, with exact fixed start. Not an error bound, a different collision trajectory, or a continuum measurement.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('fine','coarse','ultrafine','third'):
        p.add_argument(k,type=Path)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
