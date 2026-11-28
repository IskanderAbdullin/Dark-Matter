"""Independently tie both same-start D8 point responses to the varied-k baseline."""
import argparse
import json
from pathlib import Path

import numpy as np

from early_pole_endpoint_scalar_subshell_chain_receipt import read


def run(cards,clock,archive,states):
    out=[]
    for g,a,b,s in zip(('fine','third'),cards,clock,states,strict=True):
        if (a['grid'],b['grid'],a['x'],b['start'],a['k'],b['k'],a['clock_steps'],b['steps'])!=(
                g,g,.05,.04,.0015,.0015,[1024,2048],[1024,2048]):
            raise ArithmeticError('Not the two same-restart D8 chart records')
        if a['sha256_uncompressed_states']!=list(s) or not 0<=b['maximum_archived_trace_absolute_difference']<=5e-13:
            raise ArithmeticError('Accepted point hash or recovery tolerance changed')
        n=b['scalar_count_coarse_fine_relative']
        if not -1e-3<n<0:
            raise ArithmeticError('Unresolved last-segment number gate')
        group=[v for v in archive['grids'] if v['grid']==g]
        if len(group)!=1 or [v['q'] for v in a['rows']]!=([3.5e-6] if g=='fine' else [3.5e-6,1.5e-5]):
            raise ArithmeticError('Expected one fine and two third matched q values')
        for r in a['rows']:
            q=r['q']
            ref=[z for z in group[0]['rows'] if z['q']==q and abs(z['k']-.0015)<1e-15]
            if len(ref)!=1:
                raise ArithmeticError('Missing independent fixed k reference')
            low,high=[np.array(r[k]) for k in ('coarse','fine')]
            if min(low)<=0 or min(high)<=0 or not 0<=r['max_joint_mesh_defect']<=2e-6:
                raise ArithmeticError('Need both actual positive coefficients on resolved meshes')
            check=np.max(abs(high/np.array(ref[0]['coefficients'])-1))
            d=low/high-1
            if check>3e-10 or not np.array_equal(d,r['coarse_fine_relative_difference']):
                raise ArithmeticError('Fixed k point did not reproduce independent weights or paired ratios')
            out.append(dict(grid=g,q=q,relative_coarse_fine_coefficients=d.tolist(),
                         relative_coarse_fine_scalar_count=n,
                         forward_relative_shift_over_count_relative_shift=float(abs(d[0]/n)),
                         maximum_fine_point_crosscheck_relative_difference=check))
    return dict(x=.05,k=.0015,rows=out,scope='Late same-restart diagnostic at three benchmark points, not a whole-history or band error estimate.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--points',nargs=2,type=Path,required=True)
    p.add_argument('--clock',nargs=2,type=Path,required=True)
    p.add_argument('--archive',type=Path,required=True)
    p.add_argument('--states',nargs=4,type=Path,required=True)
    a=p.parse_args()
    dig=[read(f)[1] for f in a.states]
    doc=lambda f:json.loads(f.read_text())
    print(json.dumps(run(list(map(doc,a.points)),list(map(doc,a.clock)),doc(a.archive),[dig[:2],dig[2:]]),indent=2,allow_nan=False))
