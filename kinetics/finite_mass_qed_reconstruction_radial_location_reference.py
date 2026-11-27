"""Locate where conservative smoothing changes an isolated bare vertex."""
import gzip,json
from math import fsum
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre
from early_pole_decay_mode_clock_reference import Q,intervals
from finite_mass_qed_charged_closed_radial_pv_reference import angular_kernel
from finite_mass_qed_charged_overlap_asymptotic_reference import solve_axes
from finite_mass_qed_conservative_virtual_reconstruction_reference import Smoother

D=np.longdouble
LADDER=(3e-4,.02,.25)


def audit(d,kind,E,lo,hi,witness):
    M=d['x_end'];m=d['mass_ratio']*M;s=d['checkpoints'][-1]['scalar_pole_over_T']
    es=np.hypot(D(Q),D(s))
    k,ma,t0=('parent',M,D(E)) if kind=='parent' else ('fermion',m,D(E)-es)
    pp=np.asarray(d[f'endpoint_{k}_momentum_edges'],dtype=D)
    old=np.asarray(d[f'endpoint_{k}_occupations'],dtype=D)
    rec=Smoother(pp,old,ma,.05)
    edges=np.hypot(pp,D(ma))
    start,end=float(edges[0]-t0),float(edges[-1]-t0)
    center=np.searchsorted(edges,t0)-1
    ref=angular_kernel(M,m,s,E,kind)
    odd=solve_axes(M,m,s,E,kind)[1]
    shift0=(rec.value(t0)-old[center])*np.sqrt((t0-ma)*(t0+ma))/2*odd
    step_breaks=[float(x-t0) for x in edges]
    rec_breaks=[float(x-t0) for x in rec.boundaries]
    ap=[t+q for t in (lo-E,hi-E) for q in (0.,-1e-7,1e-7,-1e-5,1e-5)]
    narrow=[t for r in np.geomspace(1e-12,max(abs(start),end),48) for t in (-r,r)]
    coords=sorted({start,end,0.,*[x for x in (*rec_breaks,*step_breaks,*ap,*narrow,
                                              *LADDER,*(-np.asarray(LADDER))) if start<x<end]})

    def mesh(n):
        x,w=roots_legendre(n)
        y=np.pi*(x+1)/2
        k=w*np.sin(y)*np.pi/4
        z=np.sin(y/2)**2
        out=[[],[],[],[]]
        angular=[[],[]]
        for a,b in zip(coords[:-1],coords[1:],strict=True):
            cell=np.searchsorted(LADDER,abs((a+b)/2))
            x=D(a)+(b-a)*z.astype(D)
            v=[]
            for delta in x:
                e=t0+delta
                p=np.sqrt(max(0.,(e-ma)*(e+ma)))
                oi=np.searchsorted(edges,e)-1
                value=rec.value(e)-old[oi]
                wr=ref(delta)*p*value/2
                if cell==0:
                    wr-=shift0/delta
                v.append(float(wr))
            part=(b-a)*(k@v)
            out[cell].append(part)
            support=lo-E<(a+b)/2<hi-E
            angular[int(support)].append(part)
        return np.array([fsum(x) for x in out+angular])

    a,b=mesh(38),mesh(60)
    defect=float(max(abs(a-b)))
    expected=witness['conservative'][1]['bare_virtual']-witness['bare_steps_and_linear']['step']['value']
    closure=max(abs(fsum(b[:4])-expected),abs(fsum(b[4:])-expected))
    if max(defect,closure)>3e-5:
        raise ArithmeticError(f'Radial contribution zoning did not recover the measured shift: {defect}, {closure}.')
    return dict(cut=kind,total_conservative_minus_steps=fsum(b[:4]),
                signed_radial_band_contributions=b[:4].tolist(),
                polar_aperture_exterior_and_interior=b[4:].tolist(),resolution_defect=defect,
                independent_full_difference_defect=closure)


def run():
    root=Path(__file__).parent/'results'
    with open(root/'finite-mass-qed-conservative-virtual-reconstruction-2026-11-27.json') as f:
        wit=json.load(f)
    out=[]
    for g in ('fine','third'):
        with gzip.open(root/f'early-pole-spectrum-{g}-m06-2026-11-26.json.gz') as f:d=json.load(f)
        M=d['x_end'];m=d['mass_ratio']*M;s=d['checkpoints'][-1]['scalar_pole_over_T']
        es,lo,hi,*_=intervals(M,s,m,Q)
        w=next(z for z in wit['grids'] if z['grid']==g)
        out.append(dict(grid=g,cases=[audit(d,k,(lo+hi)/2,lo,hi,next(z for z in w['cases'] if z['cut']==k and not z['tuned_step']))
                                     for k in ('parent','daughter')]))
    return dict(scope='Local vs exterior decomposition of isolated positive bare virtual differences; no hot rate.',
                radial_band_edges_in_abs_delta_energy=[0,*LADDER,'archived_end'],rows=out)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
