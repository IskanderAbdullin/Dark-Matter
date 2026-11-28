"""Attribute elastic enhancement of a coupled slow mode to its conditional heat variance."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from boost_axis_elastic_heat_reference import integrate as scatter,leptons
from fd_portal_differential_reciprocity_reference import PositivePortal,integrate
from fd_portal_finite_capacity_modes_reference import dark
from fd_portal_reservoir_exponential_filter_reference import symmetric_root
from sm_dark_entropy_exchange_reference import StandardModel


def sample(eos,p,T,n0,cut):
    u,A=dark(T,0.)
    _,_,_,cv=eos.state(T)
    H=np.diag([1.,1./T])@np.linalg.inv(A)+np.diag([0.,1./(cv*T*T)])
    P,_=symmetric_root(H)
    _,J,_,_=integrate(p,T,0.,0.)
    S=np.diag([2.,T])
    K=S@J@S
    correl=K[1,0]/np.sqrt(K[0,0]*K[1,1])
    schur=K[1,1]-K[1,0]**2/K[0,0]
    by=np.asarray([T*scatter(T,T,0.,m,n=n0,cut=cut,linear=True) for m in leptons])
    b=np.linalg.eigvalsh(P@K@P)
    c=np.linalg.eigvalsh(P@(K+np.diag([0.,sum(by)]))@P)
    if schur<=0. or min(by)<=0. or min(b)<=0.:
        raise ArithmeticError('Conditional and flavor conductances must be positive')
    det_ratio=1.+sum(by)/schur
    closure=abs((c[0]/b[0])*(c[1]/b[1])/det_ratio-1.)
    check=abs(schur/K[1,1]-(1.-correl*correl))
    if max(closure,check)>2e-12:
        raise ArithmeticError('Eigenvalue and Schur determinant decompositions did not close')
    return dict(slow_fast_eigenvalue_ratios=(c/b).tolist(),timelike_correlation=correl,
                conditional_over_direct_timelike_heat=schur/K[1,1],
                elastic_individual_over_conditional_heat=(by/schur).tolist(),
                total_elastic_over_direct_heat=by.sum()/K[1,1],determinant_ratio=det_ratio,
                direct_spectrum_product_closure=closure)


def run(path,arch):
    content=arch.read_bytes()
    d=json.loads(content)
    h=hashlib.sha256(path.read_bytes()).hexdigest()
    if d['table_input_sha256']!=h:
        raise ArithmeticError('Archive and diagnostic EOS inputs must match')
    eos=StandardModel(path)
    p,q=[PositivePortal(o) for o in ((13,36,26),(17,44,32))]
    rows=[]
    for old in d['rows']:
        T=old['T_SM']
        a,b=[sample(eos,z,T,k,c) for z,k,c in ((p,18,29.),(q,25,36.))]
        delta=max(np.max(abs(np.asarray(a[k])/b[k]-1.)) for k in (
            'slow_fast_eigenvalue_ratios','conditional_over_direct_timelike_heat','elastic_individual_over_conditional_heat'))
        if (delta>2e-5 or np.max(abs(np.asarray(b['slow_fast_eigenvalue_ratios'])/old['combined_over_timelike']-1.))>2e-12):
            raise ArithmeticError('Two quadratures and the archived ratio do not reproduce the attribution')
        rows.append(dict(T_SM=T,maximum_relative_mesh_difference=delta,**b))
    return dict(rows=rows,archive_sha256=hashlib.sha256(content).hexdigest(),table_input_sha256=h,
                elastic_flavors=['electron','muon','tau'],
                scope='Exact local linear-response determinant decomposition within the stated vacuum portal. Not a nonlinear or freeze-out attribution.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('eos',type=Path)
    p.add_argument('modes',type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.modes),indent=2,allow_nan=False))
