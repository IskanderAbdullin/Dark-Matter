"""A positive angular-aligned discrete scalar state, not a time integration."""
import argparse
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy.special import roots_legendre

from early_pole_archived_shared_cold_window_reference import Cold
from early_pole_daughter_local_fraction_curvature_attribution import initial
from early_pole_frozen_dynamic_vlasov_reference import load


def regrid(C,raw,volume_edges,order):
    volume=np.asarray(volume_edges)
    if (volume[0]!=1e-6**3 or not np.isclose(volume[-1],1e-5**3,rtol=2e-15,atol=0.)
            or np.any(np.diff(volume)<=0.) or order not in (1,2,4)):
        raise ValueError('Expected increasing retained scalar volume edges and positive Gauss rule.')
    e=np.asarray(raw['scalar_momentum_edges'])**3
    count=np.maximum(raw['endpoint_scalar_signed_shell_number_over_T3'],0.)
    for v in e[(e>volume[0])&(e<volume[-1])]:
        if not np.any(volume==v):
            raise ValueError('A panel crosses a native scalar density discontinuity.')
    z,w=roots_legendre(order)
    a,b=volume[:-1],volume[1:]
    width=(b-a)/2
    midpoint=(a+b)/2
    idx=np.searchsorted(e,midpoint)-1
    if np.any((idx<0)|(idx>=len(count))):
        raise ValueError('Regridded nodes outside the retained cold archive.')
    weights=width[:,None]*w
    V=midpoint[:,None]+width[:,None]*z
    C.q=np.cbrt(V.ravel())
    C.sw=weights.ravel()/(6*np.pi**2)
    C.ns=(weights*(count[idx]/np.diff(e)[idx])[:,None]).ravel()
    if not np.all(np.isfinite(C.ns)&(C.ns>0.)&(C.sw>0.)):
        raise ValueError('A cold Gauss node is not admissible.')
    C.n=len(C.q)
    pole=np.hypot(C.q,C.m0)
    C.lam=2*C.m0**2/(C.ns@np.reciprocal(pole)+C.spect@np.reciprocal(np.hypot(C.spq,C.m0)))
    C.lambdas=[SimpleNamespace(n=C.M.decay.n,cells=C.M.decay.cells,q=q,nod=C.M.z,w=C.M.w) for q in C.q]
    return C


def check(archive,edges,order):
    C,seed,(a,_,b),_,dep=initial(archive)
    raw=load('results/early-pole-spectrum-fine-m06-2026-11-26.json.gz')
    mass0=C.scalar(seed)[-1]
    oldspec=C.spect.copy()
    cp=regrid(C,raw,edges,order)
    np.testing.assert_array_equal(oldspec,cp.spect)
    old_edge=np.geomspace(1e-6,1e-5,9)**3
    groups=np.searchsorted(old_edge,cp.q**3)-1
    if np.any((groups<0)|(groups>=8)):
        raise ValueError('Unexpected original shell group.')
    regroup=np.bincount(groups,weights=cp.ns,minlength=8)
    np.testing.assert_allclose(regroup,seed,rtol=0.,atol=8e-18)
    mass=cp.scalar(cp.ns)[-1]
    if abs(mass-mass0)>5e-18:
        raise ValueError('Discrete scalar gap not preserved.')
    rates=[Cold.rhs(cp,np.r_[y.ravel(),cp.ns,0.,0.],False)[cp.j0:-2] for y in (a,b)]
    return dict(scope='Cold positive integration nodes and initialized source only; not a solved candidate trajectory.',
                depth=dep,order=order,max_seed_group_difference=float(max(abs(regroup-seed))),
                gap_shift=float(mass-mass0),cold_momentum_nodes=cp.q.tolist(),
                scalar_volume_weights=cp.sw.tolist(),cold_nodes_count=cp.n,seed=cp.ns.tolist(),
                original_group=groups.tolist(),native_raw_rate=rates[0].tolist(),local_raw_rate=rates[1].tolist())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('native_archive')
    p.add_argument('refined_archive')
    p.add_argument('common_mesh')
    p.add_argument('--order',type=int,default=2)
    a=p.parse_args()
    edges=json.loads(Path(a.common_mesh).read_text())['shared_scalar_volume_edges']
    print(json.dumps([check(v,edges,a.order) for v in (a.native_archive,a.refined_archive)],
                     indent=2,allow_nan=False))
