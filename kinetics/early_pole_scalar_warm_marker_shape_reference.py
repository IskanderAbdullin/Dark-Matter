"""Audit warm shell shape, energy and inverse-energy moments before tail evolution."""
import argparse
import gzip
import json

import numpy as np
from scipy.optimize import brentq
from scipy.special import roots_legendre

from early_pole_archived_shared_cold_window_reference import Cold
from early_pole_joined_scalar_exterior_source_reference import band_kernel
from early_pole_native_aperture_witness_reference import GRID


def run(path,projection,order,compensate=False):
    with gzip.open(path,'rt') as fh:r=json.load(fh)['branches'][projection]
    grid=r['grid']
    with gzip.open(GRID/f'early-pole-spectrum-{grid}-m06-2026-11-26.json.gz','rt') as fh:a=json.load(fh)
    shot=r['snapshots'][-1]
    X=np.array(shot['hard_moments'])
    nc=np.array(shot['cold_bins'])
    C=Cold(grid,1,0,len(nc),True,exposure=r['duration'],reference_aperture_grid=r['aperture_grid'],
           archive_projection=projection)
    for x,y in zip(C.M.edge_sets,r['hard_edge_sets'],strict=True):
        np.testing.assert_array_equal(x,y)
    E,F,*_=C.M.weights(X,True)
    m0=C.scalar(nc)[-1]
    pe=np.array(a['scalar_momentum_edges'])
    selected=(C.spq>=.1)&(C.spq<=5.)&(C.spect>0.)
    q0=C.spq[selected]
    Ns=C.spect[selected]
    keep=np.r_[C.q,C.spq[~selected]]
    nk=np.r_[nc,C.spect[~selected]]
    lo,hi=pe[:-1][selected]**3,pe[1:][selected]**3
    nodes,weights=roots_legendre(order)
    z=np.cbrt((lo+hi)[:,None]/2+np.diff(np.c_[lo,hi],axis=1)*nodes/2)
    v=(hi-lo)[:,None]/(6*np.pi**2)*weights/2
    Eb=np.hypot(z,m0)
    E0=np.hypot(q0,m0)
    # Stable reduced kinetic moment about each exact archived marker.
    u=(z*z-q0[:,None]**2)/(Eb+E0[:,None])
    scale=np.ptp(Eb,axis=1)
    beta=np.empty(len(q0))
    tilted=np.zeros_like(z)
    for j,x in enumerate(u/scale[:,None]):
        w=weights/2
        def prop(b):
            a=b*x
            p=w*np.exp(a-a.max())
            return p/p.sum()
        beta[j]=brentq(lambda b:prop(b)@x,-200.,200.)
        tilted[j]=prop(beta[j])
    uniform=np.broadcast_to(weights[None,:]/2,z.shape)
    old_e=float(Ns@E0)
    old_j=float(Ns@np.reciprocal(E0))
    scope=[('volume_uniform',uniform),('energy_matched',tilted)]
    if compensate:
        scope.append(('global_nej_matched',tilted))
    base=C.scalar(nc)
    base_cold=sum(w*band_kernel(C,E,F,q,m0,f,10)[2] for q,f,w in zip(C.q,base[3],C.sw,strict=True))
    original_vol=v.sum(axis=1)
    orig_flux=sum(w*band_kernel(C,E,F,q,m0,N/w,10)[2] for q,w,N in zip(q0,original_vol,Ns,strict=True))
    out=[]
    for name,p in scope:
        weights_new=(Ns[:,None]*p).ravel()
        assert np.min(weights_new)>0.
        nodes_new=z.ravel()
        J0=weights_new@np.reciprocal(np.hypot(nodes_new,m0))
        dE=float(np.sum(Ns[:,None]*p*u))
        if name in ('energy_matched','global_nej_matched'):
            variance=np.sum(Ns[:,None]*p*u*u/(Eb*E0[:,None]**2))
            assert abs(dE)<2e-21 and J0>old_j
            assert abs(J0-old_j-variance)<1e-20
        seq=np.r_[keep,nodes_new]
        adjusted=nk.copy()
        ext={}
        if name=='global_nej_matched':
            allowed=(keep>1e-5)&(nk>0.)
            en=np.hypot(keep[allowed],m0)
            Q=np.array([np.ones(len(en)),en,1/en])
            A=Q*nk[allowed]
            norm=np.linalg.norm(A,axis=1)
            target=np.array([0.,-dE,old_j-J0])
            x=np.linalg.lstsq(A/norm[:,None],target/norm,rcond=1e-14)[0]
            assert max(abs(x))<.8
            delta=x*nk[allowed]
            adjusted[allowed]+=delta
            err=Q@delta-target
            assert np.max(abs(err))<3e-22 and np.all(adjusted>=0.)
            ext=dict(outer_adjustment_bins=int(len(x)),maximum_relative_occupation_adjustment=float(max(abs(x))),
                     total_abs_count_moved=float(np.sum(abs(delta))),
                     adjusted_other_marker_moments=(Q@delta).tolist(),moment_defect=err.tolist())
        total=np.r_[adjusted,weights_new]
        def gap(m):
            return m*m-.5*C.lam*(total@np.reciprocal(np.hypot(seq,m)))
        m=brentq(gap,1e-8,.03,xtol=1e-17)
        field=total@np.hypot(seq,m)-m**4/(2*C.lam)
        d_cold=sum(w*band_kernel(C,E,F,q,m,f,10)[2] for q,f,w in zip(C.q,base[3],C.sw,strict=True))-base_cold
        warm_flux=np.zeros(2)
        for q,pop,w in zip(nodes_new,weights_new,v.ravel(),strict=True):
            a=band_kernel(C,E,F,q,m,pop/w,10)
            warm_flux+=w*a[2:]
        out.append(dict(shape=name,partition_count_error=float(weights_new.sum()-Ns.sum()),
                        warm_energy_difference_at_fixed_pole=dE,
                        warm_inverse_energy_difference_at_fixed_pole=float(J0-old_j),
                        scalar_pole=m,pole_displacement=m-m0,
                        scalar_field_energy_difference=float(field-base[1]),
                        cold_rate_difference=float(d_cold),
                        warm_rate_and_energy=warm_flux.tolist(),**ext))
    return dict(grid=grid,projection=projection,order=order,marker_nodes=q0.tolist(),
                marker_warm_energy=old_e,marker_inverse_energy_moment=old_j,
                marker_pole=m0,marker_cold_rate=base_cold,marker_warm_rate=float(orig_flux),
                cases=out,scope='Warm shell reconstruction controls; fixed hard snapshots, no trajectory.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('path')
    p.add_argument('--projection',choices=('steps','cell-conservative'),default='steps')
    p.add_argument('--order',type=int,default=8)
    p.add_argument('--compensate',action='store_true')
    a=p.parse_args()
    if a.order<3:
        p.error('Use order at least 3.')
    print(json.dumps(run(**vars(a)),indent=2,allow_nan=False))
