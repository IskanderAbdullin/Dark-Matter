"""Scalar-occupation weighted, signed bare soft aperture at fixed frequencies."""
import argparse
import gzip
import json
from pathlib import Path
import numpy as np
from scipy.special import expit
from early_pole_frozen_yukawa_photoproduction_source_reference import parent_input
from finite_mass_qed_parent_integrated_signed_soft_reference import Source
from early_pole_decay_mode_clock_reference import Q
from early_pole_native_aperture_witness_reference import occup


def build(g,d,r):
    a=dict(d)
    p,changed=parent_input(g,d,1.)
    f0=occup((p[:-1]+p[1:])/2,np.asarray(d['endpoint_parent_momentum_edges']),
             np.asarray(d['endpoint_parent_occupations']))
    f=f0+r*(changed-f0)
    a.update(endpoint_parent_momentum_edges=p,endpoint_parent_occupations=f)
    return Source(a)


def net(k,T,ns,s,bounds,base):
    A,B,X=s.signed(k,*bounds)[0]
    a,b,x=s.signed(-k,*bounds)[0]
    inv=1/ns
    v=A-b+inv*A-base
    stim=X+x+inv*(A+a)-2*base
    return np.array([v,stim, v+stim/np.expm1(k/T)],dtype=np.longdouble)


def thermal_gate(es):
    err=0.
    for T in (.4,.9,1.3):
        for ep in (.2,1.4,3.):
            for k in (1e-9,3e-6,1e-5):
                nS,nK=1/np.expm1(es/T),1/np.expm1(k/T)
                P=expit((.13-ep)/T)
                for sign in (-1,1):
                    b=expit((.13-(ep-es-sign*k))/T)
                    A,B=P*(1-b),(1-P)*b
                    gain=A*(1+nS)*(nK+(sign==1))
                    loss=B*nS*(nK+(sign==-1))
                    err=max(err,abs(gain/loss-1.))
    if err>3e-13:raise ArithmeticError(('Common thermal exact directional null',err))
    return err


def run(g,ks,which):
    root=Path(__file__).parent/'results'
    with gzip.open(root/f'early-pole-spectrum-{g}-m06-2026-11-26.json.gz') as f:d=json.load(f)
    n=d['endpoint_scalar_occupations'][np.searchsorted(d['scalar_momentum_edges'],Q)-1]
    if n<=0:raise ArithmeticError('Outside occupied archive')
    mix=(0.,.85,1.)
    src=[build(g,d,r) for r in mix]
    thermal=float(thermal_gate(float(src[0].Es)))
    order=(24,8,16),(36,12,24),(64,18,36)
    for mesh in which:
        bounds=order[mesh]
        zero=[]
        for s in src:
            a,b,x=s.signed(0.,*bounds)[0]
            zero.append(x+a/n)
            if abs(x-s.B*s.stim)>6e-12:
                raise ArithmeticError('Signed leading Born mismatch')
        for k in ks:
            out=np.array([net(np.longdouble(str(k)),.9,n,s,bounds,a) for s,a in zip(src,zero,strict=True)])
            linear=out[0]+.85*(out[2]-out[0])
            defect=float(np.max(abs(out[1]-linear)))
            if defect>1e-11:raise ArithmeticError((g,mesh,k,defect))
            print(json.dumps(dict(grid=g,mesh=mesh,order=bounds,omega=k,mix=mix,scalar_occupation=n,
                                  leading_residues=np.asarray(zero,dtype=float).tolist(),affine_defect=defect,thermal_pointwise_defect=thermal,
                                  signed_vacuum_stimulus_total=out.astype(float).tolist())),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('grid',choices=('fine','third'))
    p.add_argument('ks',type=float,nargs='+')
    p.add_argument('--meshes',type=int,nargs='+',choices=(0,1,2),default=[0,1])
    o=p.parse_args()
    run(o.grid,o.ks,o.meshes)
