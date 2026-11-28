"""Verify an actual matched D2 adjacent window before recomputing shape order."""
import argparse
import hashlib
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_endpoint_continuum_cold_band_pilot import run as topology, FIRST_SIDE_BAND
from early_pole_endpoint_continuum_cold_band_receipt import run as base_receipt, read as raw
from early_pole_endpoint_cold_native_placement_refinement import CHANNELS, metrics
from early_pole_endpoint_cold_normalized_profile_envelope import geometry
from early_pole_endpoint_cold_native_mass_speed import evaluate as native_speed
from early_pole_endpoint_scalar_subshell_chain_receipt import read


def side(path,state,mesh,native,expected):
    rows=raw(path)
    t=topology(state,mesh,None,topology_only=True,mirror_azimuth=True,span=FIRST_SIDE_BAND)
    e=t.pop('edges')
    if (not rows or len(rows)!=expected or len(rows)!=t['total']
            or any(t[k]!=native[k] for k in ('grid','time_x','scalar_mass','omega','residue'))):
        raise ArithmeticError('Do not join unmatched endpoints or partial sidebands')
    for i,row in enumerate(rows):
        if (row['panel']!=i or row['bounds']!=e[i:i+2]
                or any(row[k]!=v for k,v in t.items())):
            raise ArithmeticError('Raw sideband does not reproduce its exact topology and hash')
    if (e[0],e[-1])!=FIRST_SIDE_BAND or np.min(np.diff(e))<=0.:
        raise ArithmeticError('The sideband is not the planned adjacent interval')
    A=np.asarray([r['coeff'] for r in rows])
    if A.shape!=(len(rows),3,3) or not np.all(np.isfinite(A)) or np.min(A[:,:,:2])<=0.:
        raise ArithmeticError('Expected positive resolved gross coefficients in each panel')
    m=t['scalar_mass']
    denom=np.maximum(1e-65,np.max(A[:,:,:2],axis=-1))
    if (np.max(abs(A[:,:,0]-A[:,:,1]-A[:,:,2])/denom)>4e-13
            or np.max(abs(m*A[:,0]+A[:,2]-A[:,1])/denom[:,1,None])>4e-13):
        raise ArithmeticError('Sideband ledger failed quantum or centered closure')
    return np.array(e),A


def run(root,archive,predecessor,accepted):
    x,y=map(Path.read_bytes,(predecessor,accepted))
    original,cert=json.loads(x),json.loads(y)
    if hashlib.sha256(x).hexdigest()!=cert['previous_sha256'] or original['original_grid_labels']!=['fine','third']:
        raise ArithmeticError('Adjacent test not bound to the actual refined native baseline')
    labels=('fine','third')
    tag='2026-11-28'
    start=[]
    for g,old,count,h in zip(labels,cert['pairs'],(10,14),original['design']['states_sha256_uncompressed'],strict=True):
        c,ch=read(root/f'early-pole-cell-endpoint-{g}-d2-x05-state-{tag}.json.gz')
        if ch!=h:
            raise ArithmeticError('Actual archive behind the refined sign changed')
        lower=[root/f'{g}-d2-placement-sign-endpoint-continuum-cold-{i}.jsonl.gz' for i in (1,2)]
        if base_receipt(lower,archive)!=old:
            raise ArithmeticError('Retained baseline did not reproduce its accepted pair')
        b=[geometry(r) for r in lower]
        u=[side(root/f'{g}-d2-endpoint-continuum-adjacent-cold-{j}.jsonl.gz',c,j,
                old['full_pole_and_snapshot'],count) for j in (1,2)]
        if not np.array_equal(u[0][0],u[1][0]):
            raise ArithmeticError('Two grids diverged in the adjacent topology')
        a,a1=[np.sum(z[1],axis=0) for z in u]
        if np.min(a1)<=0. or np.max(abs(a/a1-1.))>6e-5:
            raise ArithmeticError('Whole sideband gross and signed ledger have not resolved')
        start.append((b,u,old))
    m=start[0][2]['full_pole_and_snapshot']['scalar_mass']
    n=start[1][2]['full_pole_and_snapshot']['scalar_mass']
    out=[]
    speed=[]
    for j in (0,1):
        a,b=[z[0][j] for z in start]
        c,d=[z[1][j] for z in start]
        combined=[]
        for (g,G),(h,H) in ((a,c),(b,d)):
            if g[-1]!=h[0]:
                raise ArithmeticError('Adjacent sideband touches the retained boundary exactly')
            combined.append((np.r_[g,h[1:]],np.concatenate((G,H))))
        (e,A),(f,B)=combined
        speed.append(native_speed(e,A,f,B,m,n))
        row={}
        for k,i,s in CHANNELS:
            z=metrics(e,f,A[:,i,s],B[:,i,s],m)
            v=np.asarray(z['mean_rms_speed_shift'])
            order=['third_lower' if p[1]<-2e-12 else 'third_higher' if p[0]>2e-12
                   else 'unresolved' for p in v]
            z['signed_order']=order
            z['sideband_fraction_by_chart']=[fsum(l[1][:,i,s])/fsum(g[:,i,s])
                                            for l,g in ((c,A),(d,B))]
            row[k]=z
        out.append(row)
    for k in out[0]:
        a,b=out[0][k],out[1][k]
        if a['signed_order']!=b['signed_order'] or np.max(abs(np.asarray(a['mean_rms_speed_shift'])-b['mean_rms_speed_shift']))>2e-10:
            raise ArithmeticError('Compare resolved signs or non-resolution on both meshes')
    for label in speed[0]:
        for key,a in speed[0][label].items():
            b=speed[1][label][key]
            if (a!=b if key=='signed_order' else np.max(abs(np.asarray(a)-b))>2e-10):
                raise ArithmeticError('Native-mass adjacent comparison is not resolved across meshes')
    return dict(native_mass_mean_speed=speed[1],accepted_prior_sha256=hashlib.sha256(y).hexdigest(),retained_plus_sideband_q=[1e-6,FIRST_SIDE_BAND[1]],
                new_sideband_q=list(FIRST_SIDE_BAND),actual_native_momentum_results=out[1],
                scope='Two directly integrated adjacent D2 intervals on exactly the previous held states, without intra-panel interpolation or extrapolating past the new upper boundary.')


if __name__=='__main__':
    a=argparse.ArgumentParser(description=__doc__)
    for k in ('root','archive','predecessor','accepted'):
        a.add_argument('--'+k,type=Path,required=True)
    p=a.parse_args()
    print(json.dumps(run(p.root,p.archive,p.predecessor,p.accepted),indent=2,allow_nan=False))
