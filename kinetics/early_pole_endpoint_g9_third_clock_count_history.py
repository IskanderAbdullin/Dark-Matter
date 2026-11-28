"""Replay count sensitivity at the eight common clocks from exact completed G9 cards."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from early_pole_endpoint_scalar_subshell_chain_receipt import read


def run(report,comparison,third):
    p,ph=read(report)
    u,v=[a.read_bytes() for a in (comparison,third)]
    c,t=map(json.loads,(u,v))
    ch=hashlib.sha256(u).hexdigest()
    if (t['steps']!=[1024,2048,4096] or p['coarse_steps']!=1024 or t['input_hashes'][-1]!=ch
            or c['rows'][0]['chart']!='G9' or c['rows'][0]['sha256'][1]!=ph):
        raise ArithmeticError('Exact original report required at both accepted clocks')
    A,B,C=[np.array([[q['x'],q['scalar_yield']] for q in x]) for x in (*p['checkpoint_trajectories'],t['checks'])]
    if any(q.shape!=(8,2) for q in [A,B,C]) or not np.array_equal(A[:,0],B[:,0]) or not np.array_equal(A[:,0],C[:,0]):
        raise ArithmeticError('Clock samples differ across integrations')
    rows=[]
    for x,q,r,s in zip(A[:,0],A[:,1],B[:,1],C[:,1],strict=True):
        if min(q,r,s)<=0.:
            raise ArithmeticError('Expected positive archived total count')
        a,b=np.log(q/r),np.log(r/s)
        rows.append(dict(x=x,scalar_yield=[q,r,s],coarse_over_fine_minus_one=[q/r-1.,r/s-1.],
                         adjacent_log_refinement_order=np.log2(abs(a/b)) if a*b>0. else None))
    if max(abs(np.array(rows[-1]['coarse_over_fine_minus_one'])-np.array(t['previous_over_next_minus_one'])[:,0]))>2e-14:
        raise ArithmeticError('Final count-ratio replay mismatch')
    return dict(steps=t['steps'],rows=rows,input_sha256=[ph,ch,hashlib.sha256(v).hexdigest()],
                scope='Eight matched count samples on a single fixed-start G9 segment, not a global leading-error or order proof.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('report','comparison','third'):
        p.add_argument(k,type=Path)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
