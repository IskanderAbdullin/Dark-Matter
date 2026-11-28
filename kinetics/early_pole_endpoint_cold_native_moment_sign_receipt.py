"""Check direct sign-targeted cells on two meshes without inferring subpanel shapes."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from early_pole_endpoint_cold_native_placement_refinement import CHANNELS, metrics, project
from early_pole_endpoint_cold_normalized_profile_envelope import geometry
from early_pole_endpoint_continuum_cold_band_pilot import refine_panels
from early_pole_endpoint_continuum_cold_band_receipt import run as receipt


def run(root, archive, previous, design):
    raw=previous.read_bytes()
    x=json.loads(raw)
    w=design.read_bytes()
    d=json.loads(w)
    if d['accepted_source_sha256']!=hashlib.sha256(raw).hexdigest():
        raise ArithmeticError('The sign plan does not bind the accepted predecessor')
    plan=d['plans'][1]
    cards=x['refined_pair_receipts']
    m=cards[0]['full_pole_and_snapshot']['scalar_mass']
    accept,old,new=[],[],[]
    for g,r,p in zip(('fine','third'),cards,plan['split_plans_by_grid'],strict=True):
        base=[root/f'{g}-d2-placement-eight-endpoint-continuum-cold-{j}.jsonl.gz' for j in (1,2)]
        end=[root/f'{g}-d2-placement-sign-endpoint-continuum-cold-{j}.jsonl.gz' for j in (1,2)]
        if receipt(base,archive)!=r:
            raise ArithmeticError('Source pair differs from the accepted earlier band')
        z=receipt(end,archive)
        if z['panels']-r['panels']!=len(p) or z['full_pole_and_snapshot']!=r['full_pole_and_snapshot']:
            raise ArithmeticError('Conditional sign experiment changed the native pole or counts')
        before=list(map(geometry,base))
        after=list(map(geometry,end))
        for (e,_),(f,_) in zip(before,after,strict=True):
            if not np.array_equal(refine_panels(e,p)[0],f):
                raise ArithmeticError('The joins do not implement the exact locked plan')
        old.append(before)
        new.append(after)
        accept.append(z)
    out=[]
    for ((e,A),(f,B)),((g,G),(h,H)) in zip(zip(*old,strict=True),zip(*new,strict=True),strict=True):
        rows={}
        for label,i,j in CHANNELS:
            a,b=G[:,i,j],H[:,i,j]
            u,tu,du=project(e,A[:,i,j],g,a)
            v,tv,dv=project(f,B[:,i,j],h,b)
            if max(du,dv)>3e-8:
                raise ArithmeticError('New child totals fail their measured parent')
            val=metrics(g,h,a,b,m)
            cond=metrics(g,h,u,v,m)
            base=metrics(e,f,A[:,i,j],B[:,i,j],m)
            ww=np.array(val['mean_rms_speed_shift'])
            pp=np.array(cond['mean_rms_speed_shift'])
            bb=np.array(base['mean_rms_speed_shift'])
            eps=2e-12
            if np.max(abs(ww-pp))>eps or np.any(pp[:,1]>bb[:,1]+eps) or np.any(pp[:,0]<bb[:,0]-eps):
                raise ArithmeticError('Empirical or parent-conservative bounds did not remain nested')
            sign=['third_below_fine' if max(a[1],b[1]) < -eps else
                  'third_above_fine' if min(a[0],b[0]) > eps else
                  'unresolved' for a,b in zip(ww,pp,strict=True)]
            rows[label]=dict(actual=val,projection=cond,prior=base,signed_order=sign,
                             maximum_relative_parent_drift=[du,dv],panel_projection_TV=[tu,tv])
        out.append(rows)
    for a,b in zip(out[0].values(),out[1].values(),strict=True):
        if a['signed_order']!=b['signed_order'] or np.max(abs(np.array(a['actual']['mean_rms_speed_shift'])-
                np.array(b['actual']['mean_rms_speed_shift'])))>2e-12:
            raise ArithmeticError('Both meshes must independently resolve a reported sign')
    return dict(previous_sha256=hashlib.sha256(raw).hexdigest(),design_sha256=hashlib.sha256(w).hexdigest(),
                fields=['mean_q','RMS_q','mean_speed_at_fixed_mass'],pairs=accept,metrics=out[1],
                scope='Fixed native response on one retained band; mass distributions inside the final subpanels are otherwise unrestricted.')


if __name__=='__main__':
    a=argparse.ArgumentParser(description=__doc__)
    for k in ('root','archive','previous','design'):
        a.add_argument('--'+k,type=Path,required=True)
    p=a.parse_args()
    print(json.dumps(run(p.root,p.archive,p.previous,p.design),indent=2,allow_nan=False))
