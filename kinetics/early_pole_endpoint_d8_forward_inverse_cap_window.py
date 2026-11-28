"""Locate the geometric forward-only window from the actual retained D8 supports."""
import argparse
from functools import lru_cache
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from early_pole_cell_endpoint_transverse_audit import charge_field,continuum_moments,root
from early_pole_endpoint_past_k_cap_aperture import aperture
from early_pole_endpoint_peer_charge_shape_transport import unpack
from early_pole_endpoint_scalar_subshell_chain_receipt import read

QS=(3.5e-6,1.5e-5)


def run(card,sha):
    g=card['grid']
    d,z,m=unpack(card,g,8)
    x=.05
    mf=x*d.ratio
    charged=((x,d.p_edges,d.wp*z[d.iP:]),(mf,d.d_edges,d.wd*z[d.ns:d.iP]))
    fit=[]
    caps=[]
    catalogue=[]
    for mass,edges,count in charged:
        use=np.flatnonzero(count.reshape(-1,2).sum(1)>0)
        if len(use)==0 or np.any(count<0):
            raise ArithmeticError('Need finite occupied charge support')
        caps.append(float(np.hypot(mass,edges[use[-1]+1])))
        catalogue.append(float(np.hypot(mass,edges[-1])))
        p,w,_,_,_,_=continuum_moments(edges,mass,count,36)
        fit.append((p,w,mass))
    @lru_cache(maxsize=180)
    def pole(k):
        _,w,R=root(lambda beta:sum(charge_field(*r,beta,36) for r in fit),k)
        return w
    @lru_cache(maxsize=512)
    def get(k,q,j):
        w=pole(k)
        E=np.hypot(m,q)
        a,b=caps
        effective=(a,b+E-w)[j]
        if effective<=x:
            raise ArithmeticError('Parent-equivalent daughter support is invalid')
        return aperture(k,q,w,m,x,mf,effective)
    data=[]
    K=np.linspace(.00135,.00165,25)
    for q in QS:
        rows=[]
        brackets=[]
        for k in K:
            a,b=[get(k,q,j) for j in range(2)]
            rows.append(dict(k=k,forward_fraction=a['allowed_cosine_measure_fraction'],
                             inverse_fraction=b['allowed_cosine_measure_fraction']))
        for j in range(2):
            for target in (-1.,1.):
                fn=lambda k: get(k,q,j)['cosine_upper_before_clipping']-target
                v=np.array([fn(k) for k in K])
                pairs=np.where(v[:-1]*v[1:]<0)[0]
                for h in pairs:
                    k=brentq(fn,K[h],K[h+1],xtol=5e-16)
                    if abs(fn(k))>3e-8:
                        raise ArithmeticError('Detected support crossing did not replay')
                    brackets.append(dict(channel=('forward_parent','inverse_daughter')[j],cosine=target,
                                         sample_pair=K[h:h+2].tolist(),k_crossing=k))
        data.append(dict(q=q,samples=rows,sign_change_crossings=brackets))
    return dict(grid=g,x=x,charge_subcells=8,state_sha256_uncompressed=sha,
                parent_daughter_archive_energy_caps=catalogue,parent_daughter_nonzero_energy_caps=caps,
                measurements=data,scope='Geometric caps of positive archived charge populations; classify occupied support, not amplitude or rate. No uniqueness assumed for bracketed roots.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('state',type=Path,nargs='+')
    a=p.parse_args()
    print(json.dumps([run(*read(p)) for p in a.state],indent=2,allow_nan=False))
