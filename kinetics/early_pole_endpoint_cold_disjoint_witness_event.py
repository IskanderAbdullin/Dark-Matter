"""Construct event witnesses from the complete non-nested cold panel capacity cover."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_endpoint_cold_non_nested_cdf_contrast import audit


def segments(q,keep):
    ix=np.flatnonzero(np.diff(np.r_[False,keep,False].astype(int)))
    return np.asarray(q)[ix].reshape(-1,2).tolist()


def witness(e,f,a,b,prices):
    q=np.union1d(e,f)
    i,j=[np.searchsorted(edges,q[:-1],side='right')-1 for edges in (e,f)]
    price=[np.asarray(w,dtype=float) for w in prices]
    if any(z.shape!=(len(w),) for z,w in zip(price,(a,b),strict=True)) or any(np.max(abs(p-np.round(p)))>2e-9 for p in price):
        raise ArithmeticError('Expected an independently auditable near-integer capacity cover')
    pc,qc=[p>.5 for p in price]
    if not np.all(pc[i]|qc[j]):
        raise ArithmeticError('The integral capacity cover leaves a positive-length intersection uncovered')
    P,Q=[np.asarray(v)/fsum(v) for v in (a,b)]
    # E is union of uncovered left panels; every overlapping right panel
    # must be covered. Bound ingress into that physical set by whole masses.
    ri=np.zeros(len(Q),dtype=bool)
    li=np.zeros(len(P),dtype=bool)
    ri[j[~pc[i]]]=True
    li[i[~qc[j]]]=True
    if np.any(ri&~qc) or np.any(li&~pc):
        raise ArithmeticError('Uncovered exact witness interval crosses a non-cover panel')
    p0,p1=fsum(P[~pc]),fsum(Q[ri])
    q0,q1=fsum(Q[~qc]),fsum(P[li])
    v=1.-fsum(P[pc])-fsum(Q[qc])
    d=[p0-p1,q0-q1]
    if min(d)+2e-14<v:
        raise ArithmeticError('Direct event mass lower bound is worse than its cover')
    return dict(fine_event_segments=segments(e,~pc),third_event_segments=segments(f,~qc),
                fine_event_minimum_fine_minus_third=d[0],third_event_minimum_third_minus_fine=d[1],
                fine_and_third_exact_own_event_masses=[p0,q0],
                uncertain_opposite_event_mass_upper=[p1,q1],minimum_gap_from_cover=v)


def run(fine,third,states,archive,fixture):
    if fixture['grids']!=['fine','third'] or fixture['steps']!=2048 or fixture['charge_subcells']!=2:
        raise ArithmeticError('Expected matched native D2 result')
    a,b=[audit(p,s,archive,g) for p,s,g in zip((fine,third),states,('fine','third'),strict=True)]
    if [a[1],b[1]]!=fixture['states_sha256_uncompressed'] or [a[2],b[2]]!=fixture['gates']:
        raise ArithmeticError('Witness events do not use the accepted exact source profiles')
    out={}
    defect=0.
    for key,i,j in (('forward_number',0,0),('inverse_number',0,1),('forward_excess',2,0)):
        z=fixture['observables'][key]['sharp_compatible_TV']
        u,v=[witness(e,f,A[:,i,j],B[:,i,j],z['native_panel_capacity_prices'])
             for ((e,A),(f,B)) in zip(a[0],b[0],strict=True)]
        low,high=z['minimum_TV_bounds']
        for k in ('fine_event_minimum_fine_minus_third','third_event_minimum_third_minus_fine','minimum_gap_from_cover'):
            defect=max(defect,abs(u[k]-v[k]))
            if not low-2e-12<=v[k]<=high+2e-12:
                raise ArithmeticError('A real momentum event disagrees with the independently bracketed optimum')
        out[key]=v
    if defect>2e-9:
        raise ArithmeticError('Paired meshes disagree on disjoint event witnesses')
    return dict(grids=['fine','third'],charge_subcells=2,x=.05,
                states_sha256_uncompressed=fixture['states_sha256_uncompressed'],observables=out,
                maximum_paired_event_bound_defect=defect,
                scope='Selected unions of full native bins; bound holds for every absolutely continuous pair consistent with the panel masses.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('fine','third','states'):
        p.add_argument('--'+k,type=Path,nargs=2,required=True)
    for k in ('archive','cdf'):
        p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args()
    print(json.dumps(run(a.fine,a.third,a.states,a.archive,json.loads(a.cdf.read_text())),indent=2,allow_nan=False))
