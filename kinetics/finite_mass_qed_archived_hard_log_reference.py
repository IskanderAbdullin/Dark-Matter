"""Frozen real one-loop hard log shift on two archived isotropic shapes."""
import gzip
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad

from finite_mass_qed_shell_log_reference import gauss, integrated_kernel


def frozen_species(edges,f,m,cut,pole_mom,n):
    # Truncating the integration interval does not rebin or extrapolate f.
    first=np.searchsorted(edges,cut,side='right')-1
    edges=edges[first:].copy()
    edges[0]=cut
    f=f[first:]
    t,w=gauss(n)
    widths=np.diff(edges)
    q=edges[:-1,None]+widths[:,None]*t
    E=np.hypot(q,m)
    ww=widths[:,None]*w
    j=(ww*q*q/E).sum(axis=1)
    dj=-(ww*q*q/E*(m/E)**2).sum(axis=1)
    val=j@f
    K,D=integrated_kernel(edges,m,n=n)
    L=f@K@f
    DL=f@D@f
    max_mesh_scale=float(max(abs(K).ravel()))
    # Independent pointwise sum using weighted integrals only if p is a
    # segment endpoint; within a segment split exactly at p.
    changes=[]
    for p in pole_mom:
        point=0.
        for f0,l,r in zip(f,edges[:-1],edges[1:],strict=True):
            amp=lambda k: k/np.hypot(k,m)
            if l<=p<=r:
                plus=quad(lambda k: amp(k)*np.log(p+k),l,r,epsabs=2e-12)[0]
                minus=0.
                for sign,size in ((-1,p-l),(1,r-p)):
                    if size>0:
                        minus+=quad(lambda s: amp(p+sign*s),0,size,weight='alg-loga',
                                    wvar=(0.,0.),epsabs=2e-12)[0]
                k_int=plus-minus
            else:
                k_int=quad(lambda k: amp(k)*np.log1p(2*min(k,p)/abs(k-p)),l,r,epsabs=2e-12)[0]
            point+=f0*k_int
        changes.append(-m*m*point/(2*p))
    return dict(moment=float(val),moment_derivative=float(dj@f),log_integral=float(L),log_derivative=float(DL),
                pole_log_numerator=changes,max_log_matrix_coefficient=max_mesh_scale)


def run():
    rows=[]
    momenta=[.25,.5,1.,2.,5.]
    for grid in ('fine','third'):
        path=Path(__file__).parent/f'results/early-pole-spectrum-{grid}-m06-2026-11-26.json.gz'
        with gzip.open(path) as handle:
            d=json.load(handle)
        mp=d['x_end']; mb=mp*d['mass_ratio']
        for cut in (.1,.25,.5):
            item=[]
            for name,m in (('parent',mp),('fermion',mb)):
                edges=np.asarray(d[f'endpoint_{name}_momentum_edges'],float)
                f=np.asarray(d[f'endpoint_{name}_occupations'],float)
                a=frozen_species(edges,f,m,cut,momenta,20)
                b=frozen_species(edges,f,m,cut,momenta,36)
                drift=max(abs(a[k]/b[k]-1) for k in ('moment','moment_derivative','log_integral','log_derivative'))
                if drift>5e-9:
                    raise ValueError(f'Momentum quadrature mismatch {drift}.')
                a.update(quadrature_relative_defect=drift,species=name,mass=m)
                item.append(a)
            # Photon bath hard moments are an explicit comparator input, not a
            # field inferred from the archive. Cold/thermal photon response is
            # common and cancels from the parent-minus-daughter squared shift.
            jp=np.pi**2/6-quad(lambda q:q/np.expm1(q),0.,cut,epsabs=1e-13)[0]
            leading=(item[0]['moment']-item[1]['moment'])/np.pi**2
            delta=np.subtract(item[0]['pole_log_numerator'],item[1]['pole_log_numerator'])/np.pi**2
            L=-sum(v['log_integral'] for v in item)/4
            LJ=-sum(v['log_derivative'] for v in item)/4
            trace0=sum((jp+v['moment'])*v['moment_derivative'] for v in item)
            if leading<=0 or not all(delta<0):
                raise ValueError('Unexpected sign of the hard same-flavor difference or its log correction.')
            rows.append(dict(grid=grid,internal_momentum_min=cut,external_momenta=momenta,
                             species=item,leading_difference_over_g2=float(leading),
                             leading_plus_log_difference_over_g2=(leading+delta).tolist(),
                             log_change_relative_to_leading=(delta/leading).tolist(),
                             log_energy_over_g2=float(L/np.pi**4),
                             trace_factorizable_over_g2=float(trace0/np.pi**4),
                             trace_log_over_g2=float(LJ/np.pi**4)))
    return dict(scope='Frozen charge-even bare-on-shell hard shifts on a represented archive; not rates or evolved apertures.',
                reference_epoch_x=.05,rows=rows)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
