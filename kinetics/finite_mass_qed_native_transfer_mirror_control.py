"""Check reflected azimuth quadrature on both signs of the retained transfer."""
import argparse
import json
from time import process_time

import numpy as np

from finite_mass_qed_timelike_pole_angular_reference import angular_knots
from finite_mass_qed_timelike_pole_batched_absorption_reference import load,pole,radial


def probe(grid,timing=False):
    d,ch,_=load(grid)
    m=d['checkpoints'][-1]['scalar_pole_over_T']
    k=.0015
    w,R=pole(ch,'T',k)
    q0=np.sqrt((w-m)*(w+m))
    L,H=.0015,.004
    if not L < q0 < H:
        raise ArithmeticError('There is no interior transfer zero')
    rows=[]
    times=np.zeros(2)
    for q in (L+(q0-L)/4,(L+q0)/2,q0-1e-7,q0+1e-7,(q0+H)/2,H-(H-q0)/4):
        es=np.hypot(q,m)
        z,_=angular_knots(ch,k,w,es,q,allow_sum_boundaries=True)
        mids=(z[:-1]+z[1:])/2.
        close=[max(-.999999,float(mids[0])),.0,min(.999999,float(mids[-1]))]
        for j,c in enumerate(close):
            row=dict(q=q, transfer_energy=es-w,cosine=c)
            for n,a in ((5,12),(8,18)):
                v=[]
                for flag in ((j%2==1),(j%2==0)):
                    start=process_time()
                    out=radial(ch,k,w,R,es,c,'T',n,a,qscalar=q,
                               allow_future_spacelike=True,mirror_azimuth=flag)
                    times[int(flag)]+=process_time()-start
                    v.append(out)
                u,t=v if j%2==0 else v[::-1]
                glob=max(abs(u[0][0])+abs(u[0][1]),1e-85)
                err=float(np.max(abs(u[0]-t[0]))/glob)
                if u[1]!=t[1] or err>2e-12 or not np.all(np.isfinite([u[0],t[0]])):
                    raise ArithmeticError('Reflected azimuth does not reproduce the full stencil')
                row[f'energy_{n}_azimuth_{a}']=dict(support=bool(u[1]),error=err,
                                                    full_number_coefficients=u[0].tolist())
            rows.append(row)
    for sign in (-1,1):
        if not any(sign*r['transfer_energy']>0. and r['energy_8_azimuth_18']['support'] for r in rows):
            raise ArithmeticError('No nonempty channel tested for one transfer sign')
    out=dict(grid=grid,frozen_band_q=[L,H],mass=m,pole_omega=w,zero_transfer_q=q0,
             probe_count=len(rows),cases=rows,
             maximum_error=max(v['error'] for r in rows for name,v in r.items() if name.startswith('energy_')),
             scope='Local radial stencils across both transfer signs only. The reflected angular stencil is the same quadrature, not a new independent mesh or complete rate receipt.')
    if timing:
        out['optional_process_seconds_full_and_mirror']=times.tolist()
        out['optional_cpu_ratio_full_over_mirror']=times[0]/times[1]
    return out


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--timing',action='store_true')
    a=p.parse_args()
    print(json.dumps([probe(g,a.timing) for g in ('fine','third')],indent=2,allow_nan=False))
