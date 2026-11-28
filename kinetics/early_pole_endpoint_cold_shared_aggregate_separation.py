"""Fixed projection checks for joint aggregate rates with common S2 moments."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_endpoint_cold_joint_absolute_profile import make
from early_pole_endpoint_cold_normalized_profile_envelope import fit_state, geometry

W={
    'forward_number': [1.,0.,0.],
    'inverse_number': [0.,1.,0.],
    'forward_excess': [0.,0.,1.],
    'coupled_three': [.69030,.71998,-.07156],
    'coupled_numbers': [.90785,-.41929,0.],
}


def support(card, c, mode):
    start,_,*families=card
    p=float(fsum(start.sum(axis=1)*c))
    n=len(c)
    values=[x.reshape(len(x),n,-1).sum(axis=2)@c for x in families[mode]]
    return p+fsum(min(v) for v in values),p+fsum(max(v) for v in values)


def compare(left,right,mode):
    scale=left[1].sum(axis=1)
    ans={}
    if min(scale)<=0:
        raise ArithmeticError('Need positive reference totals')
    for name,v in W.items():
        w=np.asarray(v)
        ranges=[support(c,w/scale,mode) for c in (left,right)]
        lo,hi=ranges
        gap=max(lo[0]-hi[1],hi[0]-lo[1],0.)
        ans[name]=dict(weights=v, S2_S4_ranges=ranges, disjoint_gap=gap,
                       guaranteed_linf_fraction_at_least=gap/fsum(abs(w)))
    return ans


def run(left,right,states,archive):
    charts=[fit_state(p,archive,n,files) for p,n,files in zip(states,(2,4),(left,right),strict=True)]
    d,f,m,_=charts[0]
    lo,hi=1e-6,1.5e-5
    parent=np.r_[lo,d.ds.ed[(d.ds.ed>lo)&(d.ds.ed<hi)],hi]
    shared=parent,f[np.searchsorted(d.ds.ed,(parent[:-1]+parent[1:])/2)-1],m
    out=[]
    normalizers=[]
    for t in (0,1):
        Z=list(map(geometry,(left[t],right[t])))
        bins=np.intersect1d(Z[0][0],Z[1][0])
        if len(bins)<5 or [bins[0],bins[-1]]!=[lo,hi]:
            raise ArithmeticError('Incomplete common support')
        a,b=[make(*z,s,bins,shared) for z,s in zip(Z,charts,strict=True)]
        out.append({k:compare(a,b,j) for j,k in enumerate(('count_only','count_and_mean'))})
        normalizers.append(a[1].sum(axis=1))
    low,high=out
    defect=max(abs(x-y) for mode,a in low.items() for name,b in a.items()
               for p,q in zip(b['S2_S4_ranges'],high[mode][name]['S2_S4_ranges'],strict=True)
               for x,y in zip(p,q,strict=True))
    if defect>2e-9 or max(abs(normalizers[0]/normalizers[1]-1))>2e-9:
        raise ArithmeticError('Paired meshes disagree about aggregate response ranges')
    return dict(x=.05,k=.0015, fixed_S2_mass=m, original_S2_normalizers=normalizers[1].tolist(),
                maximum_mesh_absolute_endpoint_defect=defect, **high,
                scope='Fixed optimizer-free separating projections of finite positive panel families; not a continuum or physical bound.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('left','right','states'):
        p.add_argument('--'+k,nargs=2,type=Path,required=True)
    p.add_argument('--archive',type=Path,required=True)
    a=p.parse_args()
    print(json.dumps(run(a.left,a.right,a.states,a.archive),indent=2,allow_nan=False))
