"""Select prospective cold midpoint splits targeted at an unresolved moment sign."""
import argparse
import hashlib
import json
from math import fsum
from pathlib import Path

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp

from early_pole_endpoint_cold_native_placement_refinement import CHANNELS
from early_pole_endpoint_continuum_cold_band_receipt import run as receipt
from early_pole_endpoint_cold_normalized_profile_envelope import geometry


def targets(a,b,m,n):
    """Use equal split mass only as a design proxy, never a measured bound."""
    tags=[]
    gaps=[]
    gains=[]
    for label,i,j in CHANNELS:
        old=[]
        want=[]
        for side,(e,M) in enumerate((a,b)):
            f=e[:-1]+np.diff(e)/2
            u=M[:,i,j]
            prob=u/fsum(u)
            def probes(t):
                return np.array([t,t*t,t/np.hypot(t,m)])
            z=probes(e)
            half=probes(f)
            old.append([fsum(p*prob) for p in (z[:,1:] if side else z[:,:-1])])
            want.append(.5*prob*(z[:,1:]-half if side else half-z[:,:-1]))
        gaps.extend(np.subtract(*old[::-1]))
        gains.extend(np.concatenate(want,axis=1))
        tags.extend([label+'/'+k for k in ('mean_q','mean_q_squared','mean_speed')])
    G,V=np.asarray(gains),np.asarray(gaps)
    precision=np.maximum(np.abs(V),np.array([1e-10,1e-15,1e-8]*len(CHANNELS)))
    k=G.shape[1]
    score=G/precision[:,None]
    target=V/precision
    z=milp(np.r_[np.zeros(k),1.],integrality=np.r_[np.ones(k),0.],
           bounds=Bounds(np.r_[np.zeros(k),-np.inf],np.r_[np.ones(k),np.inf]),
           constraints=LinearConstraint(np.vstack((np.r_[np.ones(k),0.],
               np.c_[score,np.ones(len(V))])),np.r_[n,target],np.r_[n,np.repeat(np.inf,len(V))]),
           options={'mip_rel_gap':0.})
    if not z.success:
        raise ArithmeticError('No finite midpoint target plan')
    mask=np.rint(z.x[:k])
    if abs(np.sum(mask)-n)>0 or np.max(abs(mask-z.x[:k]))>1e-7:
        raise ArithmeticError('Non-integral or wrong-size midpoint plan')
    return tags,V,V-G@mask,mask


def run(root, accepted, archive):
    raw=accepted.read_bytes()
    d=json.loads(raw)
    if d['original_grid_labels']!=['fine','third'] or d['x']!=.05:
        raise ArithmeticError('Require accepted two-mesh eight-split endpoint receipt')
    cards=d['refined_pair_receipts']
    paired_by_grid=[[root/f'{g}-d2-placement-eight-endpoint-continuum-cold-{i}.jsonl.gz'
                     for i in (1,2)] for g in ('fine','third')]
    if [receipt(paths,archive) for paths in paired_by_grid]!=cards:
        raise ArithmeticError('The source no longer reproduces its accepted two-mesh receipt')
    m=cards[0]['full_pole_and_snapshot']['scalar_mass']
    output=[]
    for mesh in (1,2):
        a,b=[geometry(paths[mesh-1]) for paths in paired_by_grid]
        if [len(a[0])-1,len(b[0])-1]!=[q['panels'] for q in cards]:
            raise ArithmeticError('Accepted placement topology does not match')
        rows=[]
        for n in (4,8):
            labels,old,new,mask=targets(a,b,m,n)
            k=len(a[0])-1
            split=[[[float(edges[i]),float(edges[i]+(edges[i+1]-edges[i])/2),float(edges[i+1])]
                    for i in np.flatnonzero(v)] for edges,v in ((a[0],mask[:k]),(b[0],mask[k:]))]
            rows.append(dict(selected_original_panels=n,split_plans_by_grid=split,
                             unresolved_upper_before=old.tolist(),equal_child_mass_proxy_upper_after=new.tolist(),
                             fields=labels))
        output.append(rows)
    if [[r['split_plans_by_grid'] for r in q] for q in output][0]!=[r['split_plans_by_grid'] for r in output[1]]:
        raise ArithmeticError('Direct pairs do not give a stable sign-target design')
    return dict(accepted_source_sha256=hashlib.sha256(raw).hexdigest(),plans=output[1],scope='Targets from a balanced-child proxy only; signs require new two-mesh integrations.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True)
    p.add_argument('--accepted',type=Path,required=True)
    p.add_argument('--archive',type=Path)
    a=p.parse_args()
    arc=a.archive
    if arc is None:
        arc=a.root/'finite-mass-qed-fixed-pole-cold-scalar-band-2026-11-28.json'
        if not arc.exists():
            arc=Path('../../dm-0046-scalar-q-band-stage/kinetics/results')/arc.name
    print(json.dumps(run(a.root,a.accepted,arc),indent=2,allow_nan=False))
