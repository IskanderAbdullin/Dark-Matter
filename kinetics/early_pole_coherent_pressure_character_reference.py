"""Cycle-averaged pressure and adiabatic action for a uniform scalar oscillator."""
import argparse
import gzip
import json
from pathlib import Path

import numpy as np
from scipy.special import ellipj, ellipk, roots_legendre

from early_pole_coherent_quartic_floquet_reference import background


def average(A,m,lam,order):
    V2=.5*m*m*A*A
    V4=lam*A**4/24
    rho=V2+V4
    O=np.sqrt(m*m+lam*A*A/6)
    k=lam*A*A/(12*O*O)
    z,w=roots_legendre(order)
    # A quarter period suffices for these even observables.
    u=(1+z)*ellipk(k)/2
    sn,cn,dn,_=ellipj(u,k)
    v=-A*O*sn*dn
    K=.25*w@(v*v)
    V=.5*w@(V2*cn**2+V4*cn**4)
    vir=.5*w@(2*V2*cn**2+4*V4*cn**4)
    tr=.5*w@(m*m*(A*cn)**2)
    eps=max(abs((K+V)/rho-1),abs(2*K/vir-1),abs((rho-3*(K-V)-tr)/rho))
    assert eps<3e-12
    # Frozen oscillator action independent of the phase quadrature.
    theta=(1+z)*np.pi/4
    e=np.cos(theta)**2*(V2+V4*(1+np.sin(theta)**2))
    action=np.pi*A*(w@(np.cos(theta)*np.sqrt(2*e)))
    return dict(w=float((K-V)/rho),rho=rho,action=action,Omega=O,trace_over_rho=float(tr/rho),
                numerical_residual=float(eps))


def one(s):
    x=s['x_end']
    rho=s['checkpoints'][-1]['integrated_scalar_energy_over_T4']
    m=1e-9*x
    lam=7.65175e-8
    A,O,_=background(rho,m,lam)
    seq=np.sqrt(12/lam)*m
    hi=average(A,m,lam,96)
    hi2=average(A,m,lam,192)
    eq=average(seq,m,lam,192)
    scale=float((hi2['action']/eq['action'])**(1/3))
    assert abs(hi2['action']/hi['action']-1)<1e-11
    H=np.sqrt(np.pi**2*106.75/90)*(86./x)/2.435e18
    return dict(x=x,mass_ratio=s['mass_ratio'],initial_coherent_amplitude_over_T=A,
                initial_amplitude_over_quadratic_quartic_equality= A/seq,
                coherent_pressure_over_energy_initial=hi2['w'],
                coherent_trace_over_energy_initial=hi2['trace_over_rho'],
                coherent_pressure_over_energy_at_equality=eq['w'],
                inferred_scale_factor_ratio_in_isolated_adiabatic_model=scale,
                fast_oscillator_reference_H_over_Omega_initial=H/O,
                shape_series=[dict(amplitude_over_equality=r,**average(r*seq,m,lam,96))
                              for r in (.03,.1,.3,1.,3.,10.,30.)])


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('spectra',nargs='+')
    a=p.parse_args()
    d={}
    for name in a.spectra:
        f=Path(name)
        raw=gzip.decompress(f.read_bytes()) if f.suffix=='.gz' else f.read_bytes()
        d[f.name]=one(json.loads(raw))
    print(json.dumps(d,indent=2,allow_nan=False))
