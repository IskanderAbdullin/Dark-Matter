"""Separate the explicit scalar factor from the rate kernel on the whole cold band."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_endpoint_cold_matched_clock_profile import fitted
from early_pole_endpoint_cold_normalized_profile_envelope import geometry

FIELDS=((0,0),(0,1),(2,0))


def crossed(books,cards):
    if not np.array_equal(cards[0][0].ds.ed,cards[1][0].ds.ed):
        raise ArithmeticError('Each affected scalar cell must exist at both clocks')
    ed=cards[0][0].ds.ed
    C=[]
    for i,(q,K) in enumerate(books):
        k=np.searchsorted(ed,(q[:-1]+q[1:])/2)-1
        if np.any(k<0) or np.any(k>=len(ed)-1) or np.any(q[:-1]<ed[k]) or np.any(q[1:]>ed[k+1]):
            raise ArithmeticError('A raw panel crosses the shared native scalar cell')
        f=[c[1][k] for c in cards]
        row=[]
        for j in (0,1):
            x=[]
            for moment,inv in FIELDS:
                if inv and min(f[i])<=0:
                    raise ArithmeticError('Cannot infer inverse kernel from zero underlying occupation')
                kernel=K[:,moment,inv]/(f[i] if inv else 1+f[i])
                p=kernel*(f[j] if inv else 1+f[j])
                if np.any(p<0) or (i==j and np.max(abs(p-K[:,moment,inv])/np.maximum(K[:,moment,inv],1e-90))>4e-14):
                    raise ArithmeticError('The measured native response does not reconstruct')
                x.append(fsum(p))
            row.append(x)
        C.append(row)
    Z=np.asarray(C)
    if np.any(Z<=0):
        raise ArithmeticError('Log fractions need positive gross response')
    L=np.log(Z)
    total=L[1,1]-L[0,0]
    kernel=((L[1,0]-L[0,0])+(L[1,1]-L[0,1]))/2
    pop=((L[0,1]-L[0,0])+(L[1,1]-L[1,0]))/2
    if np.max(abs(total-kernel-pop))>2e-14 or np.any(total==0):
        raise ArithmeticError('Paired replacement paths do not close')
    return Z,dict(increasing_fine_minus_coarse_total_log=total.tolist(),
                  all_bare_kernel_log=kernel.tolist(),explicit_population_log=pop.tolist(),
                  explicit_population_share_of_total_log=(pop/total).tolist(),
                  mixed_kernel_population_log_interaction=(L[1,1]-L[1,0]-L[0,1]+L[0,0]).tolist())


def run(paths,states,archive,fixture):
    if fixture['charts']!=[2,2,4,4] or fixture['step_counts']!=[1024,2048,1024,2048]:
        raise ArithmeticError('Wrong matched restart receipt')
    c=[fitted(p,s,archive,n,t) for p,s,n,t in zip(paths,states,(2,2,4,4),(1024,2048,1024,2048),strict=True)]
    rec=fixture['full_state_panel_receipts']
    if [a[-1] for a in c]!=[a['state_sha256_uncompressed'] for a in rec]:
        raise ArithmeticError('Original matched endpoints changed')
    out=[]
    for i,n in enumerate((2,4)):
        k=2*i
        row=[]
        for t in (0,1):
            row.append(crossed([geometry(p[t]) for p in paths[k:k+2]],c[k:k+2]))
        (U,a),(V,b)=row
        err=float(np.max(abs(U/V-1)))
        if err>6e-5:
            raise ArithmeticError('Alternative gross integrated cells do not resolve')
        for j,M in enumerate((V[0,0],V[1,1])):
            tally=rec[k+j]['diagnostic']['raw_number_scalar_energy_and_excess']
            v=np.array([tally[x][y] for x,y in FIELDS])
            if np.max(abs(M/v-1))>5e-14:
                raise ArithmeticError('Original response and paired band integral differ')
        out.append(dict(scalar_subcells=n,paired_mesh_maximum_relative_defect=err,matrix_kernel_by_population=V.tolist(),**b))
    return dict(x=.05,k=.0015,clock_order=['coarse','fine'],field_order=['forward_number','inverse_number','forward_excess'],
                results=out,scope='Exact explicit scalar factor replacement on independently held kernels. Kernel bundles optical and charged state changes; no dynamical causality or extrapolation.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('s2_coarse','s2_fine','s4_coarse','s4_fine'):
        p.add_argument('--'+n.replace('_','-'),nargs=2,type=Path,required=True)
    p.add_argument('--states',nargs=4,type=Path,required=True)
    p.add_argument('--archive',type=Path,required=True)
    p.add_argument('--matched',type=Path,required=True)
    a=p.parse_args()
    print(json.dumps(run((a.s2_coarse,a.s2_fine,a.s4_coarse,a.s4_fine),a.states,a.archive,json.loads(a.matched.read_text())),indent=2,allow_nan=False))
