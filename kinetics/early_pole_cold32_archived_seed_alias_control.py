"""Distinguish the actual joined seed from frozen diagnostic point sampling."""
import gzip,json
from math import fsum
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre

D=np.longdouble


def run():
    R=Path(__file__).parent/'results'
    def rd(name):
        with gzip.open(R/name) as f:return json.load(f)
    s=rd('early-pole-spectrum-fine-m06-2026-11-26.json.gz')
    h=rd('early-pole-unsplit-hard8-cold32-native-2026-11-27.json.gz')
    src=np.asarray(s['scalar_momentum_edges'],dtype=D)
    wt=np.asarray(s['endpoint_scalar_signed_shell_number_over_T3'],dtype=D)
    fp=np.asarray(s['endpoint_scalar_occupations'],dtype=D)
    m=D(s['checkpoints'][-1]['scalar_pole_over_T'])
    e=np.geomspace(1e-6,1e-5,33)
    q=np.cbrt((e[1:]**3+e[:-1]**3)/2)
    actual=np.array(h['cold_mode_count_history'][0])
    np.testing.assert_array_equal(q,h['cold_momentum_nodes_over_T'])
    N=[]; point=[]
    X,W=roots_legendre(21)
    exact_energy_weight=[]
    for l,r,qc in zip(e[:-1],e[1:],q,strict=True):
        values=[]
        exact=[]
        for a,b,n in zip(src[:-1],src[1:],wt,strict=True):
            lo,hi=max(a,D(l)),min(b,D(r))
            if not hi>lo:continue
            frac=(hi-lo)*(hi*hi+hi*lo+lo*lo)/((b-a)*(b*b+a*b+a*a))
            count=frac*max(0.,n)
            values.append(float(count))
            u=lo+(hi-lo)*(1+X)/2
            av=(W@(u*u/np.hypot(u,m)))/(W@(u*u))
            exact.append(float(count*av))
        N.append(fsum(values))
        exact_energy_weight.append(fsum(exact))
        k=np.searchsorted(src,qc)-1
        point.append(float(fp[k]*(r**3-l**3)/(6*np.pi**2)))
    N=np.array(N)
    if max(abs(N-actual))>3e-17:
        raise ArithmeticError('Independent shell overlap failed the actual retained H8 initialization.')
    tab=fsum(N/np.hypot(q,float(m)))
    continuum=fsum(exact_energy_weight)
    return dict(scope='Audit of actual H8 native cold32 initial selection only, not time-dependent node or continuum error.',
                exact_shell_overlap_total=fsum(N),retained_actual_scalar_number=fsum(actual),
                retained_modewise_overlap_max_error=float(max(abs(N-actual))),
                hypothetical_point_sample_number_error=fsum(point)/fsum(N)-1,
                initial_mode_center_inverse_energy_weight=tab,independent_continuum_inverse_energy_weight=continuum,
                actual_initial_inverse_energy_weight_relative_error=tab/continuum-1)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
