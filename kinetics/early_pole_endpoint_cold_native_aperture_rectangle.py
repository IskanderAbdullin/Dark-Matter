"""Bound parent energy thresholds on the finite native fine/third cold configurations."""
import argparse
import json
from pathlib import Path

import numpy as np

from early_pole_endpoint_stencil_remap_pilot import load
from finite_mass_qed_timelike_pole_batched_absorption_reference import spacelike_parent_minimum


def threshold(M,m,t,v):
    s=(t-v)*(t+v)
    A=M*M+s-m*m
    if not (t<0. and v>abs(t) and A>0.):
        raise ArithmeticError('Require the past-directed spacelike opening with positive A')
    E=spacelike_parent_minimum(M,m,t,v,s)
    u=np.sqrt((E-M)*(E+M))
    k=v*u+E*t-A/2.
    if not E>M or abs(k)>2e-14*max(abs(v*u),abs(E*t),A):
        raise ArithmeticError('The parent threshold does not satisfy the independent collinear equation')
    d=v*E/u+t
    if d<=0.:
        raise ArithmeticError('The admissible parent threshold is not monotone')
    return E,-(E-t)/d,-(u+v)/d


def run(fixture):
    poles=[v['full_pole_and_snapshot'] for v in fixture['gates']]
    if [u['grid'] for u in poles]!=['fine','third'] or fixture['charge_subcells']!=2 or fixture['steps']!=2048:
        raise ArithmeticError('Use the accepted matched native D2 CDF fixture')
    ratio=[load(p['grid'])[0]['mass_ratio'] for p in poles]
    if ratio[0]!=ratio[1]:
        raise ArithmeticError('Matched reaction masses are needed for ordering')
    M,K=.05,.0015
    m=M*ratio[0]
    out=[]
    for q in (2.5e-6,3.5e-6,4.7e-6,8.2e-6):
        values=np.zeros((2,2,3,3))
        transferred=np.zeros((2,2))
        for i in (0,1):
            for j in (0,1):
                t=np.hypot(poles[i]['scalar_mass'],q)-poles[j]['omega']
                transferred[i,j]=t
                for z,c in enumerate((-1.,0.,1.)):
                    v=np.sqrt(K*K+q*q-2*K*q*c)
                    values[i,j,z]=threshold(M,m,t,v)
        E=values[:,:,:,0]
        if (not np.all(np.diff(E,axis=2)>0.) or np.any(values[:,:,:,1:]>=0.)
                or not transferred[1,1]<transferred[0,0]):
            raise ArithmeticError('The aperture ordering and derivative gates did not match')
        # At fixed q and c the same reaction masses are used. By monotonicity
        # in t the smaller native third transfer increases E for all c.
        dm=(E[1,0]-E[0,0]+E[1,1]-E[0,1])/2.
        dw=(E[0,1]-E[0,0]+E[1,1]-E[1,0])/2.
        if np.max(abs(dm+dw-E[1,1]+E[0,0]))>2e-14:
            raise ArithmeticError('Independent geometry replacement paths do not reproduce their total opening')
        out.append(dict(q=q,transfer_energy_rectangle=transferred.tolist(),
                        parent_energy_and_derivatives=values.tolist(),
                        symmetric_scalar_gap_shift=dm.tolist(),symmetric_photon_energy_shift=dw.tolist()))
    return dict(grids=['fine','third'],x=.05,k=K,reaction_mass_ratio=ratio[0],cosine_samples=[-1.,0.,1.],
                axes=['scalar_gap','photon_energy','cosine'],
                components=['E_min','d_E_min_d_transfer_energy','d_E_min_d_transfer_magnitude'],
                states_sha256_uncompressed=fixture['states_sha256_uncompressed'],points=out,
                scope='Exact collinear threshold and monotone angular endpoints, not an integral or an assumed Boltzmann rate.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('fixture',type=Path)
    a=p.parse_args()
    print(json.dumps(run(json.loads(a.fixture.read_text())),indent=2,allow_nan=False))
