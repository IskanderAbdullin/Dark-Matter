"""Replay the signed momentum-band receipts and bound the omitted pole-energy moments."""
import argparse
import gzip
import json
from math import fsum
from pathlib import Path

import numpy as np

from finite_mass_qed_frozen_transverse_plasmon_gate_reference import Q
from finite_mass_qed_timelike_pole_batched_absorption_reference import load, pole


def open_rows(path):
    if path.suffix=='.gz':
        with gzip.open(path,'rt') as r:
            return [json.loads(s) for s in r]
    return [json.loads(s) for s in path.read_text().splitlines()]


def pair(root,g,p,idx,gzip_input=False):
    end='.jsonl.gz' if gzip_input else '.jsonl'
    L=open_rows(root/f'pole-band-{g}-{p}{idx}{end}')
    n=L[0]['total']
    if len(L)!=n or [x['i'] for x in L]!=list(range(n)):
        raise ValueError(f'{g} {p} mesh {idx}: only {len(L)} of {n} contiguous records')
    if any(x['grid']!=g or x['pol']!=p or x['total']!=n for x in L):
        raise ArithmeticError('A stored row has the wrong scope')
    bounds=np.r_[[r['range'][0] for r in L],L[-1]['range'][1]]
    if any(x['range']!=bounds[j:j+2].tolist() for j,x in enumerate(L)) or np.any(np.diff(bounds)<=0):
        raise ArithmeticError('A finite panel is missing or duplicated')
    V=np.array([x['raw_signed_coefficients'] for x in L])
    C=np.array([fsum(V[:,i]) for i in range(3)])
    if np.any(V[:,:2]<0) or abs(C[0]-C[1]-C[2])>2e-13*max(C[:2]):
        raise ArithmeticError('Signed Pauli ledger inconsistent')
    return bounds,C,V,L


def run(root,gzip_input=False):
    out=[]
    for g in ('fine','third'):
        d,ch,es=load(g)
        n=d['endpoint_scalar_occupations'][np.searchsorted(d['scalar_momentum_edges'],Q)-1]
        dS=1+1/n
        for p in ('T','L'):
            b0,C0,_,_=pair(root,g,p,0,gzip_input)
            b1,C1,V,L=pair(root,g,p,1,gzip_input)
            if not np.array_equal(b0,b1):
                raise ArithmeticError('Higher quadrature changed exact momentum panels')
            rel=float(max(abs(C0/C1-1.)))
            if rel>5e-5:
                raise ArithmeticError(('Finite band mesh did not resolve both positive arms',g,p,rel))
            om=np.array([pole(ch,p,x)[0] for x in b1])
            if np.any(np.diff(om)<=0):
                raise ArithmeticError('Monotone energy bracket does not apply')
            lo,hi=(s @ V[:,:2] for s in (om[:-1],om[1:]))
            if any(lo<=0.) or dS*lo[0]<=hi[1]:
                raise ArithmeticError('Heat threshold has no checked positive slope')
            number=C1[1]/(C1[2]+C1[0]/n)
            low=lo[1]/(dS*hi[0]-lo[1])
            high=hi[1]/(dS*lo[0]-hi[1])
            margin=low/number-1.
            significant=margin>max(1e-4,10*rel)
            norm=1/(32*np.pi**3*es)
            out.append(dict(grid=g,mode=p,momentum_interval=b1[[0,-1]].tolist(),panels=len(L),
                            panel_mesh_relative_defect=rel,phase_prefactor_excluding_couplings_and_n_s=norm,
                            raw_integrated_number_coefficients=C1.tolist(),number_balancing_occupation=number,
                            lower_upper_heat_threshold=[low,high],heat_above_number_resolved_by_bounds=bool(significant),
                            full_gain_photon_energy_bounds=[lo[0]/C1[0],hi[0]/C1[0]],
                            full_inverse_photon_energy_bounds=[lo[1]/C1[1],hi[1]/C1[1]],
                            forward_only_angle_nodes=sum(x['forward_only_angle_nodes'] for x in L)))
    return out


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('root',type=Path)
    p.add_argument('--gzip-input',action='store_true')
    a=p.parse_args()
    print(json.dumps(run(a.root,a.gzip_input),indent=2,allow_nan=False))
