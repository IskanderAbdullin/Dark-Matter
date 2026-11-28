"""Frozen hard-exterior cross term of the symmetric massive log energy."""
import gzip
import json
from pathlib import Path

import numpy as np

from finite_mass_qed_shell_log_reference import gauss, integrated_kernel


def decompose(K,D,f,active,target,J):
    x=f*active
    y=f*~active
    vv=np.array([x,y])
    parts=vv@K@vv.T
    red=vv@D@vv.T
    total=f@K@f
    red_total=f@D@f
    if abs(parts.sum()/total-1)>2e-14 or abs(red.sum()/red_total-1)>2e-14:
        raise ValueError('The active/exterior log terms did not close.')
    # The diagonal and cross term are terms of the same quadratic
    # energy, but the cross term enters the active first variation.
    grad=np.stack([(K@x)[target]/J[target],(K@y)[target]/J[target]])/2
    whole=np.sum(grad,axis=0)
    if np.min(whole)<=0:
        raise ValueError('The comparison aperture lost its positive kernel weight.')
    ext_share=grad[1]/whole
    dx=.1*x*target
    derivative=2*(dx@K@y)
    h=.25
    # Independent finite difference with the exterior held fixed.
    plus=(x+h*dx+y)@K@(x+h*dx+y)-(x+h*dx)@K@(x+h*dx)
    minus=(x-h*dx+y)@K@(x-h*dx+y)-(x-h*dx)@K@(x-h*dx)
    fderr=abs((plus-minus)/(2*h)-derivative)/total
    if fderr>2e-14:
        raise ValueError('External term does not reproduce the active energy derivative.')
    return dict(probe_external_first_variation_scaled_defect=float(fderr),
                energy_component_fractions=[float(parts[0,0]/total),float(2*parts[0,1]/total),
                                            float(parts[1,1]/total)],
                trace_component_fractions=[float(red[0,0]/red_total),float(2*red[0,1]/red_total),
                                           float(red[1,1]/red_total)],
                fixed_comoving_probe_external_fraction=[float(ext_share.min()),float(ext_share.max())],
                probe_external_numerator_minmax=[float(grad[1].min()),float(grad[1].max())])


def run():
    result=[]
    for name in ('fine','third'):
        p=Path(__file__).parent/f'results/early-pole-spectrum-{name}-m06-2026-11-26.json.gz'
        with gzip.open(p) as f: d=json.load(f)
        for label in ('parent','fermion'):
            borders=np.asarray(d[f'endpoint_{label}_momentum_edges'])
            occ=np.asarray(d[f'endpoint_{label}_occupations'])
            m=d['x_end']*(d['mass_ratio'] if label=='fermion' else 1.)
            for cut in (.1,.25,.45):
                ends=np.unique(np.r_[borders[borders>cut],cut,.45,2.2,5.,10.,1.40,1.46])
                lo,hi=ends[:-1],ends[1:]
                mid=(lo+hi)/2
                f=occ[np.searchsorted(borders,mid)-1]
                probe=(lo>=1.40)&(hi<=1.46)
                z,w=gauss(24)
                q=lo[:,None]+(hi-lo)[:,None]*z
                J=(hi-lo)*((q*q/np.hypot(q,m))@w)
                K,D=integrated_kernel(ends,m,n=20)
                K2,D2=integrated_kernel(ends,m,n=32)
                qerr=max(np.max(abs(K/K2-1)),np.max(abs(D/D2-1)))
                if qerr>5e-8:
                    raise ValueError('Log matrix diagonal or shared corner quadrature was unresolved.')
                for high in (2.2,5.,10.):
                    keep=(lo>=.45)&(hi<=high)
                    row=decompose(K,D,f,keep,probe,J)
                    row.update(grid=name,species=label,hard_floor=cut,active_low=.45,active_high=high,
                               quadrature_relative_defect=float(qerr))
                    result.append(row)
    return dict(scope='Decomposition of the frozen same-species hard log energy, trace and probe variation, no external kinetic exchange.',
                component_order=['active-active','two active-exterior orientations','exterior-exterior'],
                rows=result)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
