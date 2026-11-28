"""Frozen fixed-k photon source at selected q on matched scalar-only histories."""
import argparse
import gzip
import json
from pathlib import Path

import numpy as np

from early_pole_cell_endpoint_transverse_audit import charge_field,continuum_moments,model_callback,piecewise_logit,root
from early_pole_endpoint_cold_aperture_attribution import K,coefficients
from early_pole_endpoint_scalar_subshell_pilot import ScalarSubshellChart
from early_pole_endpoint_stencil_remap_pilot import SPECS,X,load


Q=(2.6e-6,3.5e-6,4.0e-6,4.7e-6,8.2e-6)


def read(p):
    op=gzip.open if p.suffix=='.gz' else open
    with op(p,'rt',encoding='utf-8') as f:
        return json.load(f)


def card(c):
    g=c['grid']
    a=load(g)[0]
    d=ScalarSubshellChart(a,SPECS[g],c['charge_div'],c['scalar_div'],c['cut'])
    x=c['x_end']
    f=np.asarray(c['endpoint_occupation'])
    if (c['counts']!=[d.ns,d.nf,d.np] or c['reference_x']!=X or f.shape!=(d.n,) or not np.all(np.isfinite(f))
            or min(f)<0 or max(f[d.ns:])>1):
        raise ArithmeticError('Invalid state or chart')
    m=d.state(x,f)[0]
    fit=[]; ch=[]; cb=[]
    for mass,edges,N in ((x,d.p_edges,d.wp*f[d.iP:]),(x*d.ratio,d.d_edges,d.wd*f[d.ns:d.iP])):
        p,F,_,_,a,b=continuum_moments(edges,mass,N,36)
        fit.append((p,F,mass))
        ch.append(piecewise_logit(mass,edges,a,b,N,1))
        cb.append(model_callback(mass,edges,a,b))
    _,w,R=root(lambda z:sum(charge_field(p,F,m,z,36) for p,F,m in fit),K)
    measures=[]
    for q in Q:
        A,err=coefficients(ch,m,w,R,q,cb)
        i=np.searchsorted(d.ds.ed,q,side='right')-1
        s=f[i]
        B=np.array([(1+s)*A[0],s*A[1],A[0]+s*A[2]])
        if min(B)<=0 or abs(B[0]-B[1]-B[2])>7e-13*B[0]:
            raise ArithmeticError('Affine conditional source ledger failed')
        measures.append(dict(q=q,occupation=float(s),scalar_shell_edges=d.ds.ed[i:i+2].tolist(),
                             raw_affine=A.tolist(),angular_mesh_discrepancy=float(err),
                             conditional_affine=B.tolist(),conditional_photon_number_balance=float(B[1]/B[2])))
    return dict(scalar_div=c['scalar_div'],m=float(m),transverse_frequency=w,transverse_residue=R,points=measures)


def run(states):
    a,b=states
    if (a['grid'],a['x_end'],a['charge_div'],a['scalar_div'])!=('fine',.05,2,2) or (
        b['grid'],b['x_end'],b['charge_div'],b['scalar_div'])!=('fine',.05,2,4) or any(a[t]!=b[t] for t in ('steps','cut')):
        raise ArithmeticError('Need matched fine d=2 charge charts, under scalar S2→S4 only')
    cards=[card(c) for c in states]
    ratios=[]
    for r,s in zip(*(c['points'] for c in cards),strict=True):
        source=s['conditional_affine'][0]/r['conditional_affine'][0]
        scalar_factor=(1+s['occupation'])/(1+r['occupation'])
        raw_factor=s['raw_affine'][0]/r['raw_affine'][0]
        if abs(source/(scalar_factor*raw_factor)-1)>3e-15:
            raise ArithmeticError('Raw and scalar-occupation factors do not compose')
        ratios.append(dict(q=r['q'],forward_total_ratio=source,scalar_occupation_factor=scalar_factor,
                           raw_kernel_factor=raw_factor,
                           local_photon_balance_ratio=s['conditional_photon_number_balance']/r['conditional_photon_number_balance']))
    return dict(grid='fine',fixed_k=K,x=.05,cards=cards,ratios=ratios,
                scope='Frozen selected q points; neither a q integral nor actual photon populations.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('states',type=Path,nargs=2)
    args=p.parse_args()
    print(json.dumps(run([read(p) for p in args.states]),indent=2,allow_nan=False))
