"""Positive finite-recoil transverse dissipation kernel for any isotropic Fermi histogram."""

import gzip
import json
from pathlib import Path

import numpy as np

from finite_mass_qed_spacelike_exact_statistical_aperture_reference import cell
from finite_mass_qed_transverse_kinetic_noise_reference import archive
from finite_mass_qed_frozen_transverse_plasmon_gate_reference import Species


def kernel(m, k, w, p, f):
    """Return D/w, resolving the threshold layer without differences of Fermi functions."""
    if not (m >= 0 and k > 0 and 0 < w < k):
        raise ValueError("Require a nonnegative mass and positive spacelike frequency")
    p, f = np.asarray(p), np.asarray(f)
    if (p.ndim != 1 or f.ndim != 1 or len(p) != len(f)+1 or not len(f)
            or not np.all(np.isfinite(p)) or not np.all(np.isfinite(f))
            or p[0] < 0 or np.any(np.diff(p) <= 0) or np.any((f < 0) | (f > 1))):
        raise ValueError("Invalid histogram")
    a = (k-w)*(k+w)/(k*k)
    s = a*k*k
    # k*sqrt(1+4m^2/s), rationalize its difference from w when needed.
    b = np.hypot(k, 2*m/np.sqrt(a))
    e0 = (s+(2*m/np.sqrt(a))**2)/(2*(b+w))
    E = np.hypot(p,m)
    # Layer coordinates avoid forming e0+w when w is below one energy ulp.
    x = np.clip((E[:-1]-e0)/w, 0, 1)
    y = np.clip((E[1:]-e0)/w, 0, 1)
    dlow = (y-x) * (2*s + 4*a*w*(b*(x+y)/2 + w*(x*x+x*y+y*y)/3))
    # The layer's upper edge can round onto e0 for tiny w; its contribution to
    # this integral is then below one energy ulp. The layer primitive above stays resolved.
    upper = np.minimum(np.maximum(E[:-1],e0+w),E[1:])
    high = 4*a*(E[1:]-upper)*(E[1:]+upper)
    return float(f @ (dlow + high))


def run():
    data=[]
    for grid in ('fine','third'):
        r=Path(__file__).parent/'results'
        with gzip.open(r/f'early-pole-spectrum-{grid}-m06-2026-11-26.json.gz') as fd:
            d=json.load(fd)
        parent=d['x_end']
        species=[Species(m,d[f'endpoint_{label}_momentum_edges'],d[f'endpoint_{label}_occupations'])
                 for m,label in ((parent,'parent'),(parent*d['mass_ratio'],'fermion'))]
        comparison=[]
        limiting=[]
        ultrasoft_ratio=[]
        for u in (.0047,.25,.5,.9,.99,.9996):
            for k in (.0002,.0008,.00256,.01,.1,.5):
                w=k*u
                pos=sum(kernel(c.m,k,w,c.e,c.f) for c in species)
                old=sum(cell(c.m,k,w,c.e,c.f)[2] for c in species)
                if pos <= 0 or old <= 0:
                    raise ArithmeticError(('Expected occupied scattering aperture', grid, k, u))
                comparison.append(abs(pos*w/old-1))
            pos=sum(kernel(c.m,1e-12,u*1e-12,c.e,c.f) for c in species)
            soft=sum(archive(c,u,36)[1] for c in species)
            limiting.append(abs(pos/(4*soft)-1))
            k=1e-22
            w=k*u
            d_by_w=sum(kernel(c.m,k,w,c.e,c.f) for c in species)
            G,L,_=np.sum([cell(c.m,k,w,c.e,c.f) for c in species],axis=0)
            S,D=np.sum([archive(c,u,36)[:2] for c in species],axis=0)
            ultrasoft_ratio.append(abs((G+L)/(2*d_by_w)/(S/D)-1))
        data.append(dict(grid=grid,archived_exact_max_rel=max(comparison),
                         soft_limit_max_rel=max(limiting),ultrasoft_ratio_max_rel=max(ultrasoft_ratio)))
    if any(x['archived_exact_max_rel']>2e-9 or x['soft_limit_max_rel']>2e-11
           or x['ultrasoft_ratio_max_rel']>2e-11 for x in data):
        raise ArithmeticError(data)

    # Occupied high-energy islands are maximally population-inverted on their
    # lower flank. Scan delta-like one-cell distributions and random masks.
    rng=np.random.default_rng(2400)
    p=np.linspace(0,6,121)
    inverted=[]; cancellation=[]
    for m in (0.,.05,.3):
        for k in (.001,.1,.7):
            for u in (.02,.3,.8,.98):
                w=k*u
                for i in (0,1,3,11,22,49,87,119):
                    f=np.zeros(len(p)-1)
                    f[i]=1
                    pos=kernel(m,k,w,p,f)
                    raw=cell(m,k,w,p,f)[2]/w
                    if min(pos,raw) < -1e-12:
                        raise ArithmeticError((m,k,u,i,pos,raw))
                    inverted.append(pos)
                    if pos > 1e-7:
                        cancellation.append(abs(raw/pos-1.))
                f=rng.choice([0.,.1,.7,.97,1.],size=len(p)-1)
                for eps in (1e-8,1e-14,1e-22):
                    h=kernel(m,eps,u*eps,p,f)
                    if h<=0.:
                        raise ArithmeticError(('Ultrasoft positive kernel lost',m,eps,u))
    if min(inverted)!=0 or max(cancellation)>2e-8:
        raise ArithmeticError(('Synthetic aperture gate',min(inverted),max(cancellation)))
    return dict(scope="Free isotropic spacelike transverse cut only",
                archive=data, one_cell_samples=len(inverted),
                one_cell_nonzero_samples=sum(x>0 for x in inverted),
                one_cell_max_rel_direct=max(cancellation))


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
