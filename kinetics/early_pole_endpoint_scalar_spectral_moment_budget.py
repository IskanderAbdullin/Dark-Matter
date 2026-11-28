"""Audit which momentum bands carry the scalar count, pressure and kinetic moments."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre

from early_pole_endpoint_scalar_subshell_spectrum_audit import read
from early_pole_endpoint_scalar_subshell_pilot import ScalarSubshellChart
from early_pole_endpoint_stencil_remap_pilot import SPECS,X,load


CUTS=(1e-5,3e-4,.01,.1)
LABELS=('count','q','q2','kinetic_energy','pressure_numerator')


def load_card(path,n):
    c=read(path)
    if (c.get('scalar_div',1),c.get('charge_div',c.get('subcells_per_archived_cell')),
            c['grid'],c['x_end'],c['steps'],c['reference_x'],c.get('cut',1e-5))!=(n,2,'fine',.05,2048,X,1e-5):
        raise ArithmeticError('Need a matched fine scalar ladder endpoint at fixed charge d=2')
    d=ScalarSubshellChart(load('fine')[0],SPECS['fine'],2,n,1e-5)
    a=np.asarray(c['endpoint_occupation'])
    if c['counts']!=[d.ns,d.nf,d.np] or a.shape!=(d.n,) or not np.all(np.isfinite(a)) or (
            min(a)<0 or max(a[d.ns:])>1):
        raise ArithmeticError('The archive state did not match its chart')
    return d,a,d.state(.05,a)[0]


def integrate(d,a,m,n):
    q=d.ds.ed
    split=np.union1d(q,CUTS)
    if split[0]!=q[0] or split[-1]!=q[-1]:
        raise ArithmeticError('Cuts must lie in the retained scalar domain')
    L,H=split[:-1],split[1:]
    idx=np.searchsorted(q,(L+H)/2,side='right')-1
    cell=a[idx]
    idxband=np.searchsorted(CUTS,(L+H)/2,side='right')
    k,w=roots_legendre(n)
    nodes=(L+H)[:,None]/2+(H-L)[:,None]*k/2
    E=np.hypot(nodes,m)
    # Stable even at ultracold q: E-m = q²/(E+m).
    X=np.stack((nodes**0,nodes,nodes**2,nodes**2/(E+m),nodes**2/E),axis=-1)
    W=cell[:,None]*nodes**2*w[None,:]*(H-L)[:,None]/(4*np.pi**2)
    v=np.sum(W[:,:,None]*X,axis=1)
    # Analytic q-polynomial control for the first three columns.
    exact=np.stack([cell*(H-L)*sum(H**(p-1-j)*L**j for j in range(p))/(p*2*np.pi**2) for p in (3,4,5)],axis=-1)
    if np.max(abs(v[:,:3]-exact)/np.maximum(abs(exact),1e-105))>7e-14:
        raise ArithmeticError('The first three moment integrals did not match their exact polynomial forms')
    ans=np.array([[fsum(v[idxband==i,k]) for k in range(len(LABELS))] for i in range(len(CUTS)+1)])
    N=fsum(d.ds.w*a[:d.ns])
    if abs(fsum(ans[:,0])/N-1)>2e-14:
        raise ArithmeticError('The partition changed the scalar number')
    return ans


def run(paths):
    levels=(1,2,4,8,16)[:len(paths)]
    ds=[load_card(p,n) for p,n in zip(paths,levels,strict=True)]
    m=ds[1][-1]  # hold the mass fixed to isolate this spectral comparison.
    out=[]
    for n,(d,a,_) in zip(levels,ds,strict=True):
        u,v=[integrate(d,a,m,i) for i in (24,48)]
        z=np.array([fsum(v[:,j]) for j in range(len(LABELS))])
        defect=np.max(abs(u-v)/np.maximum(abs(v),1e-105))
        if defect>2e-12 or min(z)<=0:
            raise ArithmeticError('A reported band moment is unresolved')
        fractions=v/z
        if np.max(abs(np.sum(fractions,axis=0)-1))>2e-14:
            raise ArithmeticError('Moment fractions do not add to one')
        out.append(dict(scalar_div=n,measurements=v.tolist(),total=z.tolist(),
                        per_particle_mean=(z[1:]/z[0]).tolist(),band_fractions=fractions.tolist(),
                        quadrature_relative_defect=float(defect)))
    return dict(scope='Conditional moments of three numerical piecewise-constant histograms at one fixed mass, not cosmological predictions.',
                x=.05,fitted_mass_fixed=float(m),cut_upper_edges=CUTS,columns=LABELS,cards=out)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('state',type=Path,nargs='+')
    args=p.parse_args()
    if len(args.state) not in (3,4,5):
        p.error('Supply native,S2,S4[,S8[,S16]]')
    print(json.dumps(run(args.state),indent=2,allow_nan=False))
