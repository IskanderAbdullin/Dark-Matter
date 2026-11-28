"""Isolate scalar mass, photon pole, and charged histograms at cold q."""
import json
from itertools import product

import numpy as np

from finite_mass_qed_timelike_pole_angular_reference import angle
from finite_mass_qed_timelike_pole_batched_absorption_reference import (
    load,pole,spacelike_parent_minimum)


GRIDS=('fine','third')
K=.0015


def load_all():
    out={}
    for g in GRIDS:
        d,ch,_=load(g)
        m=d['checkpoints'][-1]['scalar_pole_over_T']
        w,R=pole(ch,'T',K)
        out[g]=(d,ch,m,w,R)
    if [x.m for x in out['fine'][1]]!=[x.m for x in out['third'][1]]:
        raise ArithmeticError('Hybrid audit requires fixed charged masses')
    return out


def scalar_occup(d,q):
    j=np.searchsorted(d['scalar_momentum_edges'],q,side='right')-1
    return d['endpoint_scalar_occupations'][j]


def one(fixture,chgrid,mgrid,phgrid,q):
    d,ch,_,_,_=fixture[chgrid]
    _,_,m,_,_=fixture[mgrid]
    _,_,_,w,R=fixture[phgrid]
    Es=np.hypot(m,q)
    t=Es-w
    M,mb=[x.m for x in ch]
    ts=np.array([K+q,np.hypot(K,q),K-q])
    low=angle(ch,K,w,R,Es,'T',3,(3,8),q)[0]
    high=angle(ch,K,w,R,Es,'T',6,(6,12),q)[0]
    defect=float(max(abs(high-low)/np.maximum(abs(high),1e-50)))
    if defect>8e-5 or min(high)<=0.:
        raise ArithmeticError(('Unresolved point response',q,chgrid,mgrid,phgrid,defect,high))
    emin=[spacelike_parent_minimum(M,mb,t,v,(t-v)*(t+v)) for v in ts]
    return dict(charged_histogram=chgrid,scalar_mass=mgrid,photon_pole=phgrid,
                q=q,scalar_energy=Es,photon_energy=w,transfer_energy=t,
                radial_thresholds_cosine_minus_one_zero_plus_one=emin,
                raw_A_B_D=high.tolist(),mesh_defect=defect,
                scalar_histogram_occupations={k:scalar_occup(fixture[k][0],q) for k in GRIDS})


def Shapley3(data,col):
    from itertools import permutations
    key=lambda x:tuple(x[k] for k in ('charged_histogram','scalar_mass','photon_pole'))
    values={key(x):np.log(x['raw_A_B_D'][col]) for x in data}
    amount=np.zeros(3)
    for way in permutations(range(3)):
        state=['third']*3
        v=values[tuple(state)]
        for i in way:
            state[i]='fine'
            u=values[tuple(state)]
            amount[i]+=(u-v)/6
            v=u
    total=values[('fine',)*3]-values[('third',)*3]
    if abs(sum(amount)-total)>3e-13:
        raise ArithmeticError('The ordering-averaged log changes did not close')
    return amount.tolist()


def run():
    f=load_all()
    out=[]
    for q in (2.5e-6,3.5e-6,4.5e-6,6.0e-6,8.2e-6):
        rows=[one(f,a,b,c,q) for a,b,c in product(GRIDS,repeat=3)]
        decomposition=[Shapley3(rows,j) for j in (0,1,2)]
        out.append(dict(q=q,log_effect_order=['charged_histogram','scalar_mass','photon_pole'],
                        forward_inverse_and_affine_log_effect=decomposition,
                        rows=rows))
    return out


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
