"""Exact independent oscillators as a no-transport phase-mixing comparator."""
import argparse
import gzip
import json
from pathlib import Path

import numpy as np
from scipy.special import ellipj


def one(s):
    modulation=s['initial_low_mode_amplitude']
    assert modulation>0.
    initial=1.+modulation*np.cos(np.linspace(0,2*np.pi,16384,endpoint=False))
    energy=np.mean(initial**4)/4
    t=np.array([v['t_in_Omega_units'] for v in s['rows']])
    al=[]
    for time,row in zip(t,s['rows'],strict=True):
        sn,cn,dn,_=ellipj(initial*time,.5)
        m=np.mean(initial*cn)
        v=np.mean(-initial**2*sn*dn)
        no_grad=1-(.5*v*v+.25*m**4)/energy
        al.append(dict(time=float(time),no_gradient_zero_mode_energy_loss=float(no_grad),
                       lattice_zero_mode_apparent_loss=row['inhomogeneous_energy_fraction'],
                       lattice_quadratic_hard_fraction=row['quadratic_k_ge_point6_energy_fraction']))
    both=[x for x in al if x['time']<=45]
    err=max(abs(x['no_gradient_zero_mode_energy_loss']-x['lattice_zero_mode_apparent_loss']) for x in both)
    return dict(initial_low_modulation=modulation,early_max_proxy_disagreement=float(err),
                actual_lattice_max_K_above_cut=max(r['lattice_quadratic_hard_fraction'] for r in al),
                zero_mode_proxy_lattice_max=max(r['lattice_zero_mode_apparent_loss'] for r in al),
                no_gradient_zero_mode_proxy_max=max(r['no_gradient_zero_mode_energy_loss'] for r in al),rows=al)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('lattices',nargs='+')
    a=p.parse_args()
    d={}
    for name in a.lattices:
        f=Path(name)
        raw=gzip.decompress(f.read_bytes()) if f.suffix=='.gz' else f.read_bytes()
        d[f.name]=one(json.loads(raw))
    print(json.dumps(d,indent=2,allow_nan=False))
