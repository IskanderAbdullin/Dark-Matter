"""Audit within-packet Gauss convergence without changing collision support."""
import argparse
import json
from pathlib import Path

import numpy as np


def compare(base, pair):
    if any(len(a)!=4 for a in (base,*pair)):
        raise ValueError('Expected all four hard/cold fixtures in every receipt.')
    values=[]
    names=('pointwise_archive_paired_total','reconstructed_paired_total',
           'reconstructed_curvature_fixed_packet_geometry')
    sizes=[]
    for j,k,l,n in zip(base,*pair,strict=True):
        a,b,c=(q['quadrature_order'] for q in (k,l,n))
        if (a,b,c)!=(12,24,48):
            raise ValueError('Expected specified strictly nested order values.')
        for q in (j,k,l,n):
            if any(q[f]!=j[f] for f in ('hard_depth','cold_modes','daughter_cells','columns')):
                raise ValueError('Different physical support or source diagnostics.')
        packet=np.array([k['sample_counter']//a,l['sample_counter']//b,n['sample_counter']//c])
        if np.any(packet!=packet[0]) or any(q['sample_counter']%q['quadrature_order']!=0 for q in (k,l,n)):
            raise ValueError('Quadrature changed collision support.')
        if j['sample_counter']%packet[0]!=0:
            raise ValueError('Old default quadrature samples do not share the same packet support.')
        a0,a1,a2,a3=(np.asarray(q['aperture_terms_by_mode_and_cell']) for q in (j,k,l,n))
        if any(z.shape!=a0.shape for z in (a1,a2,a3)):
            raise ValueError('Incomplete mode/cell contents.')
        if not np.all(np.isfinite((a0,a1,a2,a3))):
            raise ValueError('Nonfinite data.')
        val={}
        for tag,d in (('default_to_48',a3-a0),('12_to_24',a2-a1),('24_to_48',a3-a2)):
            val[tag]=dict(per_mode_cell_max_by_diagnostic=np.max(abs(d),axis=(0,1)).tolist(),
                          max_after_summing_daughter_cells=np.max(abs(d.sum(axis=1)),axis=0).tolist())
            sizes.append(val[tag]['per_mode_cell_max_by_diagnostic'][0])
        values.append(dict(hard_depth=j['hard_depth'],cold_modes=j['cold_modes'],
                           default_order=int(j['sample_counter']//packet[0]),packet_count=int(packet[0]),
                           reference_totals={name:n[name] for name in names},differences=val))
    return dict(scope='Numerical within-packet Gauss refinement at identical frozen momentum bounds; not aperture-cut or physical convergence.',
                cases=values,maximum_pointwise_mode_cell_change=float(max(sizes)))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('base')
    p.add_argument('refined',nargs=3)
    a=p.parse_args()
    r=lambda f:json.loads(Path(f).read_text())
    print(json.dumps(compare(r(a.base),[r(f) for f in a.refined]),indent=2,allow_nan=False))
