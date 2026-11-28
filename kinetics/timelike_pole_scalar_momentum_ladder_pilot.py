"""Response of fixed photon pole A/B to resolving more than one cold scalar shell."""
import argparse
import json
import numpy as np

from finite_mass_qed_frozen_transverse_plasmon_gate_reference import Q
from finite_mass_qed_timelike_pole_angular_reference import angle
from finite_mass_qed_timelike_pole_batched_absorption_reference import load,pole,radial


def run(g):
    d,ch,_=load(g)
    m=d['checkpoints'][-1]['scalar_pole_over_T']
    q=np.array(d['scalar_momentum_edges'])
    n=d['endpoint_scalar_occupations']
    points=(0.,Q/2,Q,Q*2,4e-5,1e-4,3e-4,8e-4,.0012)
    for p,k in (('T',.0015),('T',.00256),('L',.0012),('L',.00256)):
        w,R=pole(ch,p,k)
        for z in points:
            E=np.hypot(m,z)
            if not E<w:raise ArithmeticError('Not on the past-directed energy arm')
            low,high=((6,(6,12)),(10,(14,24))) if z>=8e-4 else ((3,(3,8)),(6,(6,12)))
            A,_,_,_=angle(ch,k,w,R,E,p,*low,z)
            B,S,N,D=angle(ch,k,w,R,E,p,*high,z)
            e=max(abs(A-B))/max(1e-25,max(abs(B)))
            if e>3e-6:raise ArithmeticError((g,p,k,z,e,A,B))
            n0=None if z==0. else n[np.searchsorted(q,z)-1]
            if z==0.:
                variance=max(float(max(abs(radial(ch,k,w,R,E,c,p,*high[1],qscalar=0.)[0]-B))/
                                   max(1e-25,max(abs(B)))) for c in (-1.,-.25,.4,1.))
                if variance>5e-13:raise ArithmeticError('A scalar at rest must remove angular dependence')
            else:
                variance=None
            nB=B[1]/(B[2]+B[0]/n0) if n0 and B[1]>0 and B[2]+B[0]/n0>0 else None
            print(json.dumps(dict(grid=g,pol=p,k=k,q=z,E_scalar=E,omega=w,scalar_occup=n0,
                                  raw_A_B_D=B.tolist(),isotropic_conditional_population=nB,
                                  accepted_panels=N,exact_angle_defect=D,mesh_defect=e,meshes=[low,high],
                                  at_rest_angle_defect=variance),allow_nan=False),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('grid',choices=('fine','third'))
    run(p.parse_args().grid)
