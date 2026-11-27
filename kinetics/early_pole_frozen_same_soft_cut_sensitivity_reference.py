"""Prescribed screened Q-cut sensitivity of the initial scalar source tangent."""
import argparse
import json
import numpy as np

from early_pole_live_aperture_witness_reprojection_reference import initial_projection
from early_pole_native_aperture_witness_reference import PROFILES
from early_pole_same_flavor_screened_pole_reference import SelfSoft


def run(grid,subdivision,level):
    profile=next(x for x in json.loads(PROFILES.read_text()) if x['grid']==grid)
    C,_,Y0,Y1,*_=initial_projection(profile,subdivision,16)
    scalar=slice(C.j0,C.j0+C.n)
    out=[]
    for Q in (.04,.08,.16):
        qual=[(1,1,2),(1,2,3)][level]
        K=SelfSoft(C.M,grid,qual,Q)
        rs=[]
        for mix in (0.,.97):
            Y=(1-mix)*Y0+mix*Y1
            G,F,E,*_=C.M.weights(Y,True)
            D,S,W=K.collision(G,F,E)
            z=np.r_[Y.ravel(),C.ns,0.,0.]
            tangent=np.r_[D.ravel(),np.zeros(C.n+2)]
            ref=cold_source(C,z,scalar)
            h=(2e-14,1e-13,4e-13)
            ds=[(cold_source(C,z+step*tangent,scalar)-cold_source(C,z-step*tangent,scalar))/(2*step)
                for step in h]
            assert np.all(np.isfinite(ds))
            assert np.max(abs(np.array(ds)-ds[0]))<.04*max(1.,abs(ds[0])),ds
            parts=[]
            for group in C.M.r[1:]:
                part=np.zeros_like(D)
                part[group]=D[group]
                vv=np.r_[part.ravel(),np.zeros(C.n+2)]
                parts.append((cold_source(C,z+h[1]*vv,scalar)-cold_source(C,z-h[1]*vv,scalar))/(2*h[1]))
            assert abs(sum(parts)-ds[1])<.03*max(1.,abs(ds[1]))
            rs.append(dict(mix=mix,scalar_source=ref,finite_difference_steps=h,
                           scalar_source_time_tangents=ds,entropy=S,elastic_energy_residual=W,
                           parent_and_daughter_tangents=parts))
        out.append(dict(Q_cut=Q,packet_counts=[int(x[-1].size) for x in K.graph],rows=rs))
        del K
    return dict(grid=grid,hard_subdivision=subdivision,scalar_bins=16,Q_quadrature=qual,cases=out)


def cold_source(C,z,scalar):
    return float(C.rhs(z,False)[scalar].sum())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('grid',choices=('fine','third'))
    p.add_argument('--subdivision',type=int,default=3)
    p.add_argument('--level',type=int,choices=(0,1),default=0)
    a=p.parse_args()
    print(json.dumps(run(a.grid,a.subdivision,a.level),indent=2,allow_nan=False))
