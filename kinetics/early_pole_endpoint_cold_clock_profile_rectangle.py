"""Bound the normalized profile change over the two-by-two endpoint rectangle."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_endpoint_cold_matched_clock_profile import common_partition,fitted
from early_pole_endpoint_cold_normalized_profile_envelope import geometry


def row(q,books,i,j):
    arr=[]
    for e,z in books:
        k=np.searchsorted(q,(e[:-1]+e[1:])/2)-1
        if np.any(k<0) or np.any(k>=len(q)-1) or np.any(e[:-1]<q[k]) or np.any(e[1:]>q[k+1]):
            raise ArithmeticError('Original q panel crosses a shared bin')
        p=np.bincount(k,weights=z[:,i,j],minlength=len(q)-1)
        arr.append(p/fsum(p))
    a,b,c,d=arr # S2 coarse, S4 coarse, S2 fine, S4 fine
    def TV(a,b):
        return fsum(abs(a-b))/2
    delta=np.asarray([b-a,d-c])
    k=np.argmax(abs(delta),axis=1)
    within=(TV(a,c),TV(b,d))
    all_cross=np.asarray([[TV(u,v) for v in (b,d)] for u in (a,c)])
    if np.ptp(all_cross)>sum(within)+1e-15:
        raise ArithmeticError('TV triangle bound failed')
    return dict(S2_and_S4_single_chart_clock_TV=list(within),all_cross_TV=all_cross.tolist(),
                maximum_clock_mixing_uncertainty=sum(within),
                signed_S4_minus_S2_by_bin=delta.tolist(),largest_matched_bin_index=k.tolist(),
                ratio_of_clock_uncertainty_to_lowest_cross_TV=sum(within)/np.min(all_cross))


def run(pairs,states,archive,prior):
    records=[fitted(p,s,archive,n,t) for p,s,n,t in zip(pairs,states,(2,4,2,4),(1024,1024,2048,2048),strict=True)]
    if prior['ordered_state_sha256_uncompressed']!=[r[-1] for r in records]:
        raise ArithmeticError('Mismatched original normalization fixture')
    names=('forward_number','inverse_number','forward_excess')
    ret=[]
    for k in (0,1):
        books=[geometry(p[k]) for p in pairs]
        q=common_partition(books)
        if q.tolist()!=prior['common_q_edges']:
            raise ArithmeticError('The exact display bins changed')
        ret.append(dict(zip(names,[row(q,books,i,j) for i,j in ((0,0),(0,1),(2,0))],strict=True)))
    low,hi=ret
    err=max(np.max(abs(np.asarray(low[k]['all_cross_TV'])-np.asarray(hi[k]['all_cross_TV']))) for k in hi)
    if err>2e-9:
        raise ArithmeticError('Clock surface unresolved across accepted meshes')
    for k,u in hi.items():
        a=u['all_cross_TV']
        for z,ref in ((a[0][0],prior['coarse']),(a[1][1],prior['fine'])):
            if abs(z-ref[k]['original_S2_S4_TV'])>2e-12:
                raise ArithmeticError('The native matched TV diagonals do not close')
    return dict(x=.05,rows=hi,common_edges=prior['common_q_edges'],absolute_paired_mesh_difference=float(err),
                scope='Normalized raw profiles and metric triangle bounds on fixed common q bins. Within-final-clock stability does not address earlier history.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('s2_coarse','s4_coarse','s2_fine','s4_fine'):
        p.add_argument('--'+k.replace('_','-'),nargs=2,type=Path,required=True)
    p.add_argument('--states',nargs=4,type=Path,required=True)
    p.add_argument('--archive',type=Path,required=True)
    p.add_argument('--profile',type=Path,required=True)
    a=p.parse_args()
    print(json.dumps(run((a.s2_coarse,a.s4_coarse,a.s2_fine,a.s4_fine),a.states,a.archive,json.loads(a.profile.read_text())),indent=2,allow_nan=False))
