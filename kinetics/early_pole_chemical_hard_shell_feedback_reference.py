"""Finite-volume hard-shell feedback using analytic outgoing-angle intervals."""
import argparse
import gzip
import json
from pathlib import Path

import numpy as np

from early_pole_live_parent_energy_reference import state_moments
from early_pole_scalar_angular_chemical_coupled_reference import AngularChemical
from early_pole_scalar_number_change_angular_reference import data
from early_pole_scalar_number_change_threshold_reference import LAMBDA, threshold


def one(s,extra,power):
    shape={(67,112,122):(65,90,72),(97,156,176):(95,125,108)}
    d=AngularChemical(shape[tuple(s['cells'])],s['mass_ratio'],0.,power)
    xm=s['x_end']
    fixed=np.sqrt(3)*s['checkpoints'][-1]['scalar_pole_over_T']
    s,meta=d.setup_hot_grid(s,extra,(.8*fixed,1.2*fixed))
    f=np.asarray(s['endpoint_scalar_occupations'])
    E,D,P,m=state_moments(d,xm,np.r_[f,s['endpoint_fermion_occupations'],s['endpoint_parent_occupations']])
    E=E[:d.ns]
    edges=d.d.ed
    vol=d.d.w
    H=np.sqrt(np.pi**2*106.75/90)*(86./xm)/2.435e18
    hs=H*np.sqrt((d.radiation+D+P)/d.radiation)
    q=d.d.p
    cold=np.flatnonzero((edges[:-1]>=1e-6)&(edges[1:]<=1e-5))
    I,M,lab,angle=data(tuple(cold),tuple(q),power)
    energy=E[I].sum(axis=1)
    S=energy[:,None]**2-lab**2
    beta=np.sqrt(1-4*m*m/S)
    angular_span=.5*lab*beta
    eb=np.hypot(m,edges)
    k=np.flatnonzero((eb[1:]>=np.min(energy[:,None]/2-angular_span))&
                  (eb[:-1]<=np.max(energy[:,None]/2+angular_span)))
    offset=k[0]
    h=len(k)
    cut=eb[offset:k[-1]+2]
    J=np.zeros((h,h))
    net=0.
    norm=threshold()['event_coefficient_times_lambda4_over_m8']*LAMBDA**4/m**8/hs
    assert min(f[I].ravel())>0.
    cf=np.expm1(np.log1p(1/f[I]).sum(axis=1))
    W=norm*M*np.prod(vol[I]*f[I]*m/E[I],axis=1)/len(angle)
    eps=0.
    for i,(a,b) in enumerate(zip(energy/2,angular_span,strict=True)):
        p=np.sort(np.r_[cut,2*a-cut])
        low=a-b
        high=a+b
        L=np.maximum(p[:-1],low[:,None])
        U=np.minimum(p[1:],high[:,None])
        dr=np.maximum(U-L,0.)
        if not dr.any(): continue
        sub=dr>0
        t=(L[sub]+U[sub])*.5
        j=np.searchsorted(cut,t)-1
        k1=np.searchsorted(cut,2*a-t)-1
        assert min(np.min(j),np.min(k1))>=0 and max(np.max(j),np.max(k1))<h
        factor=W[i]*(beta[i]/np.sqrt(.75)/(2*b))[np.nonzero(sub)[0]]*dr[sub]
        u,v=f[offset+j],f[offset+k1]
        net+=np.sum(factor*(1+u+v-cf[i]*u*v))
        deriv1=factor*(1-cf[i]*v)/vol[offset+j]
        deriv2=factor*(1-cf[i]*u)/vol[offset+k1]
        np.add.at(J,(j,j),deriv1);np.add.at(J,(k1,j),deriv1)
        np.add.at(J,(j,k1),deriv2);np.add.at(J,(k1,k1),deriv2)
        eps=max(eps,abs((t+(2*a-t))-2*a).max())
    eigen=np.linalg.eigvals(J)
    assert max(abs(eigen.imag))<1e-7*max(eigen.real)
    return dict(x=xm,angular_samples=2**power,hard_receiver_bins=h,
                exact_continuum_pair_energy_identity_error=eps,
                frozen_initial_formal_loss_per_log_a=float(2*net/(vol@f)),
                frozen_formal_hard_count_growth_eigenvalue_per_log_a=float(max(eigen.real)),**meta)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('spec')
    a=p.parse_args()
    path=Path(a.spec)
    raw=gzip.decompress(path.read_bytes()) if path.suffix=='.gz' else path.read_bytes()
    s=json.loads(raw)
    results=[]
    for hot,power in [(0,5),(36,5),(100,5),(300,5),(300,7),(700,7)]:
        results.append(one(s,hot,power))
    print(json.dumps(results,indent=2,allow_nan=False))
