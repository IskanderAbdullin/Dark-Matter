"""Merge only complete compatible panel logs for the normalized signed band."""
import argparse
import gzip
import json
from math import fsum
from pathlib import Path

import numpy as np


def reduce(path,g,mesh):
    rows={}
    op=gzip.open if path.suffix=='.gz' else open
    with op(path,'rt',encoding='utf-8') as stream:
        for line in stream:
            row=json.loads(line)
            if row['i'] in rows:
                raise ArithmeticError('Repeated/refitted panel')
            if row['grid']!=g or row['mesh']!=mesh or row['common_thermal_defect']>3e-13:
                raise ArithmeticError('Panel scope or thermal gate failed')
            rows[row['i']]=row
    size=rows[0]['total_panels']
    if set(rows)!=set(range(size)):
        raise ValueError(f'{g} mesh {mesh}: {len(rows)} of {size} completed')
    rows=[rows[i] for i in range(size)]
    if rows[0]['bounds'][0]!=3e-6 or rows[-1]['bounds'][-1]!=1e-5 or any(
            a['bounds'][1]!=b['bounds'][0] for a,b in zip(rows[:-1],rows[1:],strict=True)):
        raise ArithmeticError('Finite frequency coverage is not contiguous')
    if any(r['orders']!=rows[0]['orders'] or r['total_panels']!=size for r in rows):
        raise ArithmeticError('Panel metadata changed')
    if any(not np.all(np.isfinite(r['parts'])) for r in rows):
        raise ArithmeticError('Non-finite signed contribution')
    raw=np.array([r['parts'] for r in rows])
    both=np.array([[fsum(raw[:,k,j]) for j in range(3)] for k in range(2)])
    other=both[0]+.85*(both[1]-both[0])
    return both,other,len(rows),rows


def run(root,packed=False):
    out=[]
    for g in ('fine','third'):
        a,b=[reduce(root/(f'witness-{g}-band-mesh{mesh}.jsonl'+('.gz' if packed else '')),g,mesh) for mesh in (0,1)]
        if [r['bounds'] for r in a[3]]!=[r['bounds'] for r in b[3]]:
            raise ArithmeticError('Coarse and refined panels diverge')
        dif=float(np.max(abs(b[0]-a[0])))
        if dif>4e-8 or np.any(np.sign(b[0][:,2])!=np.sign(a[0][:,2])):
            raise ArithmeticError(('Unresolved frequency-band quadrature',g,dif))
        native,endpoint=b[0][:,2]
        if endpoint==native:
            raise ArithmeticError('The witness interpolation is degenerate')
        cross=-native/(endpoint-native)
        before=-a[0][0,2]/(a[0][1,2]-a[0][0,2])
        full_panel=np.asarray([z['parts'] for z in b[3]])[:,:,2]
        following=full_panel[:,0]+.85*(full_panel[:,1]-full_panel[:,0])
        if b[1][2]<=0. or (native<0)!=(g=='fine'):
            raise ArithmeticError(('Finite-band sign pattern not reproduced',g,native,b[1][2]))
        out.append(dict(grid=g,panels=b[2],native=b[0][0].tolist(),witness_085=b[1].tolist(),
                        endpoint=b[0][1].tolist(),mesh_absolute_defect=dif,
                        zero_mixture_in_this_band=float(cross) if 0<cross<1 else None,
                        absolute_root_change_on_refinement=float(abs(cross-before)) if 0<cross<1 else None,
                        native_positive_and_negative_panel_counts=[int(sum(full_panel[:,0]>0)),int(sum(full_panel[:,0]<0))],
                        witness_positive_and_negative_panel_counts=[int(sum(following>0)),int(sum(following<0))],
                        scope='Signed vacuum remainder on a finite frequency interval, no gauge, no physical rate.'))
    return out


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('root',nargs='?',type=Path,default=Path(__file__).parent/'results')
    p.add_argument('--packed',action='store_true')
    o=p.parse_args()
    print(json.dumps(run(o.root,o.packed),indent=2,allow_nan=False))
