"""Measure archive versus reconstructed daughter occupations on the same fixed packets."""
import argparse
import gzip
import json
from math import fsum

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.special import expit

from early_pole_decay_mode_clock_reference import G,KAPPA
from early_pole_daughter_local_fraction_curvature_attribution import initial
from early_pole_frozen_dynamic_vlasov_reference import load
from early_pole_live_pauli_gauge_moment_decay_reference import Tool


def compare(bundle):
    C,n,(y0,ym,y1),_,depth=initial(bundle)
    M=C.M
    with gzip.open(bundle,'rt') as f:
        grid=json.load(f)['branches']['steps']['grid']
    archive=load(f'results/early-pole-spectrum-{grid}-m06-2026-11-26.json.gz')
    p=np.asarray(archive['endpoint_fermion_momentum_edges'])
    f=np.asarray(archive['endpoint_fermion_occupations'])
    y=p**3/3.
    P=PchipInterpolator(y,np.r_[0.,np.cumsum(f*np.diff(y))],extrapolate=False).derivative()
    rows=[]
    e,F0=M.weights(y0,True)[:2]
    F1=M.weights(y1,True)[1]
    FM=M.weights(ym,True)[1]
    E,_,_,fs,m=C.scalar(n)
    x=.05*np.exp(getattr(M,'eta',0.))
    scale=C.sw*G*KAPPA*x*((1.06*x)**2-m*m)/(C.q*E)
    dg=M.r[2]
    cells=[i for i in range(dg.start,dg.stop)
           if M.cells[i].lo>=1.4138888888888888-2e-14 and M.cells[i].hi<=1.4666666666666666+2e-14]
    cells_by_q=[]
    sample_counter=0
    for j,es in enumerate(E):
        terms=np.zeros((len(cells),3))
        for i,k,u,v,w in Tool.events(C.lambdas[j],dict(coeff=e[M.start:],es=es)):
            i+=M.start
            k+=M.start
            if k not in cells:
                continue
            A=M.cells[k]
            q=A.pfrombare(A.ec+A.de*v)
            a=np.searchsorted(p,q)-1
            if np.any((a<0)|(a>=len(f))):
                raise ValueError('Out of the exact archived daughter support.')
            U=np.maximum(P(q**3/3.),0.)
            if np.any(~np.isfinite(U)) or np.any(U>1.):
                raise ValueError('Invalid interpolation on accepted archive shell.')
            Q=f[a]
            ref,mid,to=[expit(h[k,0]+h[k,1]*v) for h in (F0,FM,F1)]
            parent=expit(F0[i,0]+F0[i,1]*u)
            weight=-w*(parent+fs[j])*scale[j]
            terms[cells.index(k)]+=[float(weight@d) for d in (U-Q,to-ref,mid-ref-.5*(to-ref))]
            sample_counter+=len(q)
        cells_by_q.append(terms)
    if not sample_counter:
        raise ValueError('Expected active cold packets over the archived daughter aperture.')
    v=np.asarray(cells_by_q).sum(axis=1).T
    direct,fit,curve=map(fsum,v)
    return dict(scope='Keep native packet energy and parent occupation; evaluate exact archive and fitted daughter shapes on the same sampled daughter momenta.',
                hard_depth=depth,cold_modes=len(n),daughter_cells=cells,sample_counter=sample_counter,
                pointwise_archive_paired_total=direct,reconstructed_paired_total=fit,
                reconstructed_minus_pointwise_archive=fit-direct,
                reconstructed_curvature_fixed_packet_geometry=curve,
                columns=['pointwise_archive_pair','fitted_moment_pair','fitted_midpoint_nonlinearity'],
                aperture_terms_by_mode_and_cell=np.asarray(cells_by_q).tolist())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('archives',nargs='+')
    args=p.parse_args()
    print(json.dumps([compare(a) for a in args.archives],indent=2,allow_nan=False))
