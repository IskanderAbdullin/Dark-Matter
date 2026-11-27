"""A dilute scalar 2->2 collision clock on the retained Hartree snapshots.

A dimensional diagnostic, not a transport rate in a highly occupied medium.
"""
import argparse
import gzip
import json
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre


def sigma_v(q,m,lam,order):
    E=np.hypot(q,m)
    p=q[:,None]
    k=q[None,:]
    ep=E[:,None]
    ek=E[None,:]
    zw,ww=roots_legendre(order)
    bulk=p*p*ek/(ep+m)+m*k*k/(ek+m)
    out=np.zeros_like(bulk)
    for z,w in zip(zw,ww,strict=True):
        d=np.maximum(bulk-p*k*z,0.)
        s=4*m*m+2*d
        mol=np.sqrt(d*(2*m*m+d))/(ep*ek)
        out+=w*.5*lam**2*mol/(32*np.pi*s)
    return out


def analyze(d,M):
    row=d['checkpoints'][-1]
    x=row['x']
    q=np.asarray(d['scalar_momentum_nodes'])
    n=np.maximum(np.asarray(d['endpoint_scalar_signed_shell_number_over_T3']),0.)
    f=np.asarray(d['endpoint_scalar_occupations'])
    m=row['scalar_pole_over_T']
    lam=7.65175e-8
    T=M/x
    h_over_T=np.sqrt(np.pi**2*106.75/90)*T/2.435e18
    W=sigma_v(q,m,lam,24)
    V=sigma_v(q,m,lam,48)
    cc=(W@n)/h_over_T
    fine=(V@n)/h_over_T
    count=n.sum()
    assert np.max(abs(V-V.T))<1e-13*V.max()
    keep=n>1e-7*count
    assert max(abs(cc[keep]/fine[keep]-1))<1e-4
    ir=np.asarray(d['scalar_momentum_edges'])[1:]<=1e-5
    hot=np.asarray(d['scalar_momentum_edges'])[:-1]>=.1
    total=fine@n
    weight=fine*n/total
    return dict(x=x,parent_mass_GeV=M,temperature_GeV=T,lambda_scalar=lam,
                scalar_hartree_mass_over_T=m,scalar_count_over_T3=count,
                number_weighted_dilute_Born_collision_clock_over_Href=float(total/count),
                sub1e5_number_weighted_dilute_clock_over_Href=float(fine[ir]@n[ir]/n[ir].sum()),
                fraction_of_dilute_collisions_on_scalar_endpoints_q_ge_point1=float(weight[hot].sum()),
                lambda_times_peak_occupation=float(lam*f.max()),
                lambda_times_number_weighted_occupation=float(lam*(f@n)/count),
                number_over_in_medium_mass_cubed=float(count/m**3),
                maximum_quadrature_drift_active_cells=float(max(abs(cc[keep]/fine[keep]-1))))


def analytic_nr_check():
    q=np.array([1e-5,3e-5,8e-5])
    m=.1
    p=q[:,None]/m
    k=q[None,:]/m
    mean=((p+k)**3-abs(p-k)**3)/(6*p*k)
    lam=1e-7
    target=lam**2*mean/(128*np.pi*m*m)
    ans=sigma_v(q,m,lam,96)
    error=float(np.max(abs(ans/target-1)))
    assert error<5e-6
    return error


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('spectra',nargs='+')
    p.add_argument('--parent-mass-GeV',type=float,default=86.)
    a=p.parse_args()
    assert a.parent_mass_GeV>0
    results=dict(nonrelativistic_angle_average_relative_test=analytic_nr_check())
    for name in a.spectra:
        f=Path(name)
        raw=gzip.decompress(f.read_bytes()) if f.suffix=='.gz' else f.read_bytes()
        results[f.name]=analyze(json.loads(raw),a.parent_mass_GeV)
    print(json.dumps(results,indent=2,allow_nan=False))
