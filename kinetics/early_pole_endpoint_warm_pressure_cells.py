"""Localize the W2 pressure change on native cells; separate count and shape."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_endpoint_scalar_subshell_spectrum_audit import read
from early_pole_endpoint_scalar_warm_evolution_audit import unpack
from early_pole_endpoint_warm_observable_quantiles import parts


def partition(before,after):
    if (before.shape!=after.shape or before.ndim!=2 or before.shape[1]!=2
            or not np.all(np.isfinite([before,after])) or np.min([before,after])<0
            or np.any((before[:,0]==0)&(before[:,1]!=0))
            or np.any((after[:,0]==0)&(after[:,1]!=0))):
        raise ValueError('Need finite nonnegative [count,pressure] pairs')
    n,p=before.T
    m,q=after.T
    b0=np.divide(p,n,out=np.zeros_like(p),where=n>0)
    a0=np.divide(q,m,out=np.zeros_like(q),where=m>0)
    # When a bin is absent in one histogram no intrabin shape comparison
    # is defined. Assign the change entirely to count, using the present mean.
    b=np.where(n>0,b0,a0)
    a=np.where(m>0,a0,b0)
    d=np.stack((q-p,(m-n)*(.5*(a+b)),(.5*(m+n))*(a-b)),axis=-1)
    if np.any(abs(d[:,0]-d[:,1]-d[:,2])>2e-14*np.maximum(p,q)):
        raise ArithmeticError('Symmetric decomposition did not close on a native cell')
    return d


def integrate(a,y,b,z,m,order):
    old=a.ds.ed
    new=b.ds.ed
    ix=np.searchsorted(new,old)
    if new[0]!=old[0] or new[-1]!=old[-1] or np.any(ix>=len(new)) or not np.array_equal(new[ix],old):
        raise ArithmeticError('The refined spectrum is not strictly nested')
    raw=parts(new[:-1],new[1:],z[:b.ns],m,order)[:,(0,2)]
    v=np.array([[fsum(raw[i:j,k]) for k in (0,1)] for i,j in zip(ix[:-1],ix[1:],strict=True)])
    u=parts(old[:-1],old[1:],y[:a.ns],m,order)[:,(0,2)]
    for grid,occ,pairs in ((a,y,u),(b,z,v)):
        if abs(fsum(pairs[:,0])/fsum(grid.ds.w*occ[:grid.ns])-1)>2e-14:
            raise ArithmeticError('Coarsening lost scalar number')
    return u,v,partition(u,v)


def run(base,warm,x):
    (a,y),(b,z)=unpack(base,False,x),unpack(warm,True,x)
    if any(not np.array_equal(getattr(a,k),getattr(b,k)) for k in ('qd','qp','wd','wp')):
        raise ArithmeticError('The comparison must keep the charge chart fixed')
    m=a.state(x,y)[0]
    lo,hi=[integrate(a,y,b,z,m,n) for n in (24,48)]
    u,v,d=hi
    baseP=fsum(u[:,1])
    scales=(u,v,np.maximum(u[:,1],v[:,1])[:,None])
    defect=max(np.max(abs(c-f)/np.maximum(s,1e-105)) for c,f,s in zip(lo,hi,scales,strict=True))
    if defect>2e-12:
        raise ArithmeticError('A cell difference did not resolve')
    totals=lambda v:[fsum(c)/baseP for c in v.T]
    edges=a.ds.ed
    untouched=edges[:-1]<.01
    if np.max(abs(d[untouched,2])/np.maximum(u[untouched,1],1e-105))>1e-13:
        raise ArithmeticError('Unexpected shape change on an unrefined cell')
    intervals=[]
    for i in np.argsort(-abs(d[:,0]))[:8]:
        intervals.append(dict(q=edges[i:i+2].tolist(),delta_pressure_components_relative_to_native_total=(d[i]/baseP).tolist(),
                              old_new_count=[u[i,0],v[i,0]], old_new_pressure=[u[i,1],v[i,1]]))
    return dict(x=x,mass_fixed=float(m),labels=('pressure_change','count_at_mean_pressure_per_particle',
                   'changed_pressure_per_particle_at_mean_count'),
                relative_total=totals(d),relative_from_whole_native_cells_starting_above_point01=totals(d[~untouched]),
                sum_absolute_cellwise_relative=totals(abs(d)),top_eight_cells=intervals,
                quadrature_relative_to_gross_cell_scale=defect,
                scope='Symmetric algebraic partition of fixed-mass spectral moments on nested numerical histograms, not a decomposition of the collision history.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('native',type=Path)
    p.add_argument('refined',type=Path)
    p.add_argument('--x',type=float,choices=(.01,.02,.04,.05),default=.05)
    a=p.parse_args()
    print(json.dumps(run(read(a.native),read(a.refined),a.x),indent=2,allow_nan=False))
