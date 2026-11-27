"""Spatial correlation and frozen-pole dephasing scales, not phase evolution."""
import argparse
import gzip
import json
from pathlib import Path

import numpy as np

from early_pole_coherent_quartic_floquet_reference import one as coherent_mode


def one(spec):
    s=spec['checkpoints'][-1]
    q=np.asarray(spec['scalar_momentum_nodes'])
    n=np.maximum(np.asarray(spec['endpoint_scalar_signed_shell_number_over_T3']),0.)
    m=s['scalar_pole_over_T']
    en=np.hypot(q,m)
    K=q*q/(en+m)
    H=n/en
    norm=H.sum()
    h=H/norm
    cold=np.asarray(spec['scalar_momentum_edges'])[1:]<=1e-5
    C=h[cold].sum()
    h0=h[cold]/C
    cs=K[cold]
    mean=h0@cs
    width=np.sqrt(h0@((cs-mean)**2))
    fc=coherent_mode(spec,86.)
    qres=fc['peak_q_over_T']
    mu=fc['maximum_floquet_exponent_over_Omega']*fc['oscillator_Omega_over_T']
    lam=2*np.pi/qres
    distances=[]
    times=[]
    for z in (1.,10.,100.):
        r=z*lam
        t=z/mu
        Csp=float(h@np.sinc(q*r/np.pi))
        # Exact sampled frozen pole sum and an independent cold-moment bound
        # that is insensitive to very rapidly oscillating hot terms.
        frozen=abs(np.sum(h*np.exp(1j*(K-mean)*t)))
        bound=max(0.,C*(1-.5*(width*t)**2)-(1-C))
        assert frozen+1e-13>=bound
        # sinc(x)>=1-x^2/6, and all omitted hot terms are >=-1.
        spbound=max(-1.,C*(1-r*r*(h0@(q[cold]**2))/6)-(1-C))
        assert Csp+1e-13>=spbound
        distances.append(dict(distance_in_conditional_resonant_wavelengths=z,
                              normalized_equal_time_2point=Csp,lower_bound=float(spbound)))
        times.append(dict(time_in_conditional_floquet_inverse_rates=z,
                          frozen_pole_temporal_envelope=float(frozen),lower_bound=float(bound)))
    return dict(x=spec['x_end'],mass_ratio=spec['mass_ratio'],m_over_T=m,
                cold_equal_time_covariance_fraction=float(C),cold_kinetic_spread_over_T=float(width),
                conditional_peak_wavenumber_over_T=qres,conditional_peak_amplitude_rate_over_T=mu,
                spatial=distances,temporal=times)


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
