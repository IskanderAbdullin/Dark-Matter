"""Enclose monotone momentum probes using only fixed native cold-panel masses."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_endpoint_cold_non_nested_cdf_contrast import audit
from early_pole_endpoint_cold_normalized_profile_envelope import fit_state,geometry


def bounds(edges,area,m):
    area=np.asarray(area,dtype=float)
    edges=np.asarray(edges,dtype=float)
    if area.ndim!=1 or not np.all(np.isfinite(area)) or not fsum(area)>0. or not np.isfinite(m):
        raise ArithmeticError('Require positive finite cold response and a finite fixed mass')
    u=area/fsum(area)
    if (len(u)+1!=len(edges) or np.any(np.diff(edges)<=0.) or min(u)<0. or not m>0.
        or abs(fsum(u)-1.)>3e-14 or min(edges)<0.):
        raise ArithmeticError('Invalid native response histogram')
    if not np.all(np.isfinite(edges)) or not np.all(np.isfinite(u)):
        raise ArithmeticError('Nonfinite native support')
    bands=[np.array([fsum(z*u) for z in (v[:-1],v[1:])])
           for v in (edges,edges**2,edges/np.hypot(edges,m))]
    bands[1]=np.sqrt(bands[1])
    return np.asarray(bands)


def pair(left,right,m):
    shape=[]
    for (e,A),(f,B) in zip(left,right,strict=True):
        r={}
        for name,i,j in (('forward_number',0,0),('inverse_number',0,1),('forward_excess',2,0)):
            a,b=[bounds(q,v[:,i,j],m) for q,v in ((e,A),(f,B))]
            r[name]=dict(left=a.tolist(),right=b.tolist(),sharp_right_minus_left=(b-a[:,::-1]).tolist())
        shape.append(r)
    a,b=shape
    err=max(float(np.max(abs(np.asarray(a[k]['sharp_right_minus_left'])-z['sharp_right_minus_left']))) for k,z in b.items())
    if err>2e-12:
        raise ArithmeticError('Paired support moments did not resolve')
    return b,err


def run(root,archive):
    tag='2026-11-28'
    fixture=json.loads((root/f'early-pole-cold-fine-third-d2-native-cdf-{tag}.json').read_text())
    m=fixture['gates'][0]['full_pole_and_snapshot']['scalar_mass']
    a,b=[audit([root/f'{g}-d2-endpoint-continuum-cold-{i}.jsonl.gz' for i in (1,2)],
               root/f'early-pole-cell-endpoint-{g}-d2-x05-state-{tag}.json.gz',archive,g) for g in ('fine','third')]
    if [a[1],b[1]]!=fixture['states_sha256_uncompressed'] or [a[2],b[2]]!=fixture['gates']:
        raise ArithmeticError('The archived matched native reference changed')
    z,e=pair(a[0],b[0],m)
    states=[([a[1],b[1]],z,e,('fine_D2','third_D2'))]
    G=[]
    for n in (2,4,8,16):
        f=[root/f'fine-charge2-scalar{n}-endpoint-continuum-cold-{i}.jsonl.gz' for i in (1,2)]
        s=root/f'early-pole-fine-charge2-scalar{n}-x05-state-{tag}.json.gz'
        records=fit_state(s,archive,n,f)
        G.append(([geometry(x) for x in f],records[-1],n))
    for (a,u,j),(b,v,k) in zip(G[:-1],G[1:],strict=True):
        z,e=pair(a,b,m)
        states.append(([u,v],z,e,(f'S{j}',f'S{k}')))
    return dict(x=.05,reference_mass=m,observables=['mean_q','rms_q','mean_speed_at_fixed_mass'],
                pairs=[dict(levels=labels,states_sha256_uncompressed=sha,rows=z,maximum_paired_mesh_absolute_defect=e)
                       for sha,z,e,labels in states],
                scope='Sharp independent placement intervals at fixed response masses; speed at one common numerical mass, not particle streaming or a physical trajectory.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('root','archive'):
        p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args()
    print(json.dumps(run(a.root,a.archive),indent=2,allow_nan=False))
