"""Integrate the frozen-parent tagged gain as a window diagnostic, without scalar feedback."""
import argparse
import gzip
import json

import numpy as np
from scipy.optimize import brentq
from scipy.special import roots_legendre

from early_pole_archive_scalar_pole_factorization_reference import residue
from early_pole_native_aperture_witness_reference import GRID, occup


def run(grid, radial_order, q_order, exposure):
    with gzip.open(GRID/f'early-pole-spectrum-{grid}-m06-2026-11-26.json.gz','rt') as f:
        a=json.load(f)
    q0,q1=1e-6,1e-5
    M,m=.05,.003
    ms=a['checkpoints'][-1]['scalar_pole_over_T']
    qedge=np.asarray(a['scalar_momentum_edges'])
    f=np.maximum(0.,np.asarray(a['endpoint_scalar_occupations']))
    cuts=[q0,q1,*qedge[(q0<qedge)&(qedge<q1)]]
    for name,mass,center in [('fermion',m,(M*M-m*m-ms*ms)/2),
                             ('parent',M,(M*M-m*m+ms*ms)/2)]:
        for p in a[f'endpoint_{name}_momentum_edges']:
            e=np.hypot(p,mass)
            for s in (-1.,1.):
                def edge(q):
                    return e*np.hypot(q,ms)+s*p*q-center
                if edge(q0)*edge(q1)<0.:
                    cuts.append(brentq(edge,q0,q1,xtol=5e-18))
    cuts=np.unique(cuts)
    x,w=roots_legendre(q_order)
    h=np.diff(cuts)[:,None]/2
    q=((cuts[1:]+cuts[:-1])[:,None]/2+h*x).ravel()
    w=(h*w).ravel()*q*q*occup(q,qedge,f)/(2*np.pi**2)
    v0=0.
    for a0,a1,y in zip(qedge[:-1],qedge[1:],f,strict=True):
        lo,hi=max(a0,q0),min(a1,q1)
        if hi>lo:
            v0+=(hi**3-lo**3)*y/(6*np.pi**2)
    np.testing.assert_allclose(w.sum(),v0,rtol=5e-14)
    rows=[]
    for r,weight in zip(q,w,strict=True):
        x=residue(a,float(r),radial_order,np.array([0.,exposure]))
        initial=x['exact_archived_initial_gross_and_net_rate'][-1]
        gain=x['integrated_gross_and_net_gain'][-1][-1]
        final=x['instant_gross_and_net_rate'][-1][-1]
        rows.append(dict(q=r,scalar_count_weight=weight,initial_rate=initial,final_rate=final,
                         frozen_parent_gain=exposure*initial,transient_gain=gain))
        assert initial>0 and final>0 and gain>0.
    f0=sum(r['scalar_count_weight']*r['frozen_parent_gain'] for r in rows)
    f1=sum(r['scalar_count_weight']*r['transient_gain'] for r in rows)
    return dict(grid=grid,exposure=exposure,q_interval=[q0,q1],q_order=q_order,radial_order=radial_order,
                initial_count=v0,frozen_parent_gain=f0,transient_parent_gain=f1,relative_change=f1/f0-1,
                q_panel_edges=cuts.tolist(),quadrature_rows=rows)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('grid',choices=('fine','third'))
    p.add_argument('--radial-order',type=int,default=24)
    p.add_argument('--q-order',type=int,default=6)
    p.add_argument('--exposure',type=float,default=7.2e-5)
    a=p.parse_args()
    if not 0.<a.exposure<=7.2e-5 or min(a.radial_order,a.q_order)<2:
        p.error('Use positive finite exposure no larger than 7.2e-5 and order >=2.')
    print(json.dumps(run(**vars(a)),indent=2,allow_nan=False))
