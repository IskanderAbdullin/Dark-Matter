"""Check exact derivatives of the past-spacelike parent threshold on fine late cards."""
import argparse
import json
from pathlib import Path

import numpy as np

from early_pole_endpoint_cold_threshold_enclosure import energy
from early_pole_endpoint_stencil_remap_pilot import load


def row(a,b,q,c):
    M=a['x']
    mf=M*load('fine')[0]['mass_ratio']
    k=a['k']
    m,w=a['scalar_mass'],a['frequencies'][1]
    v=np.hypot(q*np.sqrt((1-c)*(1+c)),q*c-k)
    def root(m,w):
        return energy(M,mf,w-np.hypot(m,q),v)
    E=root(m,w)
    u=w-np.hypot(m,q)
    p=np.sqrt((E-M)*(E+M))
    bottom=v*E/p-u
    dw=(E+u)/bottom
    dm=-m/np.hypot(m,q)*dw
    deta=1e-10
    nDw=(root(m,w+deta)-root(m,w-deta))/(2*deta)
    nDm=(root(m+deta,w)-root(m-deta,w))/(2*deta)
    check=float(max(abs(nDw/dw-1),abs(nDm/dm-1)))
    if check>2e-7 or min(dw,-dm)<=0:
        raise ArithmeticError('Implicit derivative disagrees with independently perturbed root')
    newm,neww=b['scalar_mass'],b['frequencies'][1]
    masses=(newm-m,neww-w)
    parts=dm*masses[0],dw*masses[1]
    exact=root(newm,neww)-E
    return dict(q=q,cosine=c,threshold_before=E,threshold_after=E+exact,derivative_w=float(dw),
                derivative_m=float(dm),derivative_relative_difference=check,
                linear_threshold_shift_from_mass_and_frequency=[float(k) for k in parts],
                exact_threshold_shift=exact,relative_linearization_difference=sum(parts)/exact-1)


def run(a,b):
    u,v=a[0],b[0]
    if (u['grid'],u['subcells'],v['grid'],v['subcells'],u['x'],v['x'],u['k'],v['k'])!=(
            'fine',4,'fine',8,.05,.05,.0015,.0015):
        raise ArithmeticError('Expected the fine four and eight late optical cards')
    results=[row(u,v,q,c) for q in (3.5e-6,1.5e-5) for c in (0.,.9)]
    return dict(points=results,scope='Local kinematic derivative at fitted poles; not a bound or an extension of the occupied tail.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('d4',type=Path)
    p.add_argument('d8',type=Path)
    a=p.parse_args()
    print(json.dumps(run(json.loads(a.d4.read_text()),json.loads(a.d8.read_text())),indent=2,allow_nan=False))
