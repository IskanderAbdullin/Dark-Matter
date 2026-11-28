"""Measure one same-restart cold chart at two final time steps."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_endpoint_cold_normalized_profile_envelope import geometry
from early_pole_endpoint_cold_scalar_band_ladder import verify_pair
from early_pole_endpoint_scalar_subshell_chain_receipt import read


def profile(pair, i, j):
    a,b=pair
    e=np.intersect1d(a[0],b[0])
    if [e[0],e[-1]]!=[1e-6,1.5e-5] or len(e)<6:
        raise ArithmeticError('Need exact retained common support')
    select=[]
    for q,y in pair:
        jn=np.searchsorted(e,(q[1:]+q[:-1])/2)-1
        if np.any(q[:-1]<e[jn]) or np.any(q[1:]>e[jn+1]):
            raise ArithmeticError('Profile panel crosses the shared bin')
        data=np.bincount(jn,weights=y[:,i,j],minlength=len(e)-1)
        select.append(data/fsum(data))
    r=fsum(abs(select[0]-select[1]))/2
    return r,e.tolist()


def run(coarse,fine,states,recovery,archive):
    start,low,high,original_report=states
    sha=[read(p)[1] for p in (start,high,original_report)]
    if recovery['sha256_source']!=sha or recovery['steps']!=[1024,2048] or recovery['x']!=.05 or recovery['max_trace_defect']>5e-13:
        raise ArithmeticError('The original restart trace does not match')
    r,s=[verify_pair(lines,state,archive,2,steps=n)
         for lines,state,n in zip((coarse,fine),(low,high),(1024,2048),strict=True)]
    u,v=(q['diagnostic'] for q in (r,s))
    A,B=(np.asarray(q['raw_number_scalar_energy_and_excess']) for q in (u,v))
    def get(record):
        p=record['full_pole_and_snapshot']
        return [p[k] for k in ('scalar_mass','omega','residue')]
    if [get(u),get(v)]!=recovery['optical_parameters_coarse_fine']:
        raise ArithmeticError('The integrated cold bands did not use the recovered pole and chart')
    name=['forward_number','inverse_number','forward_excess']
    tables=[]
    for k in (0,1):
        G=list(map(geometry,(coarse[k],fine[k])))
        card=[profile(G,i,j) for i,j in ((0,0),(0,1),(2,0))]
        if any(q[1]!=card[0][1] for q in card):
            raise ArithmeticError('Different channel q support')
        tables.append((card[0][1],dict(zip(name,(q[0] for q in card),strict=True))))
    (e,X),(f,Y)=tables
    d=max(abs(X[k]-Y[k]) for k in name)
    if e!=f or d>2e-9:
        raise ArithmeticError('Normalized clock change is unresolved')
    return dict(x=.05,chart='fine-D2/S2',coarse_fine_steps=[1024,2048],
                original_restart_and_fine_sha256=sha,verification=[r,s],
                same_restart_yield_relative_change=recovery['scalar_yield_coarse_fine'][0]/recovery['scalar_yield_coarse_fine'][1]-1,
                coarse_over_fine_minus_one=(A/B-1).tolist(),
                exactly_common_q_edges=f,normalized_profile_TV=Y,paired_mesh_maximum_TV_difference=d,
                scope='Final .04-to-.05 clock test on one S2 chart with fixed restart and paired q meshes; no earlier error bound.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--coarse',nargs=2,type=Path,required=True)
    p.add_argument('--fine',nargs=2,type=Path,required=True)
    p.add_argument('--states',nargs=4,type=Path,required=True,metavar=('START','COARSE','FINE','FINE_REPORT'))
    p.add_argument('--recovery',type=Path,required=True)
    p.add_argument('--archive',type=Path,required=True)
    a=p.parse_args()
    print(json.dumps(run(a.coarse,a.fine,a.states,json.loads(a.recovery.read_text()),a.archive),indent=2,allow_nan=False))
