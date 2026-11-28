"""Reconstruct comoving scalar and stable entropy on the stored two-temperature path."""
import argparse
import gzip
import json
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from three_species_decay_onsager_reference import m, mp
from scalar_isolated_entropy_saturation_reference import thermodynamic_entropy


def fermi(r,phase):
    def occ(u):
        z=u*u/(np.hypot(u,r)+r)-np.log(phase)
        e=np.exp(-z)
        return z,e/(1+e),e
    edge=200+40*np.sqrt(max(1.,2*r))
    def w(u,key):
        z,f,e=occ(u)
        if key=='entropy':
            c=np.log1p(e)+z*f
        elif key=='n':
            c=f
        elif key=='energy':
            c=np.hypot(u,r)*f
        else:
            c=u*u*f/(3*np.hypot(u,r))
        return u*u*c
    v={k:quad(lambda u:w(u,k),0.,edge,epsabs=1e-14,
              epsrel=3e-10)[0] for k in ('entropy','n','energy','pressure')}
    entropy=v['entropy']/v['n']
    check=(v['energy']+v['pressure'])/v['n']-r-np.log(phase)
    assert abs(entropy/check-1)<3e-10
    return entropy


def run(hist):
    with gzip.open(hist,'rt') as h:
        arms=json.load(h)['arms']
    res=[]
    for a in arms:
        rows=[]
        for b in a['trajectory']:
            l=fermi(m/b['light_scalar_T'][0],b['phase_factors'][0])
            s=thermodynamic_entropy(mp/b['light_scalar_T'][1],b['phase_factors'][1])
            y1,y2=b['yields']
            rows.append(dict(T_visible_GeV=b['T_visible'],
                             stable_S_over_N=l,scalar_S_over_N=s,
                             stable_entropy_per_visible_entropy=y1*l,
                             scalar_entropy_per_visible_entropy=y2*s,
                             total_dark_entropy_per_visible_entropy=y1*l+y2*s))
        fld=('stable_entropy_per_visible_entropy','scalar_entropy_per_visible_entropy',
             'total_dark_entropy_per_visible_entropy','scalar_S_over_N')
        changes={k:rows[-1][k]/rows[0][k]-1 for k in fld}
        diffs=np.diff([q[fld[2]] for q in rows])
        assert min(diffs)>0
        res.append(dict(nodes=a['nodes'],fractional_endpoint_changes=changes,
                        minimum_adjacent_total_entropy_increment=float(min(diffs)),
                        sampled_rows=rows))
    return dict(scope="Reconstructed ideal FD/BE thermodynamic entropy at the archived points "
                      "of a reduced approximate collision model. Positive sampled total "
                      "entropy drift is not an H-theorem proof for the MB closure or the "
                      "underlying unresolved momentum distributions.",results=res)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('history',type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.history),indent=2,allow_nan=False))
