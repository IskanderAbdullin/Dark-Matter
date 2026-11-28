"""Initial cold-source curvature along a convex local hard-moment replacement."""
import argparse
import gzip
import json
from math import fsum

import numpy as np

from early_pole_archived_shared_cold_window_reference import Cold


def run(bundle):
    with gzip.open(bundle,'rt') as f:
        branches=json.load(f)['branches']
    a,b=(branches[k] for k in ('steps','cell-conservative'))
    depth=a.get('aperture_subdivision',1)
    if depth!=b.get('aperture_subdivision',1) or depth not in (1,4):
        raise ValueError('Unmatched aperture depth.')
    z,y=(np.asarray(x['snapshots'][0]['hard_moments']) for x in (a,b))
    n0=np.asarray(a['snapshots'][0]['cold_bins'])
    C=Cold(a['grid'],1,0,len(n0),True,exposure=a['duration'],
           reference_aperture_grid=a['aperture_grid'],archive_projection='steps',
           aperture_subdivision=depth)
    M=C.M
    np.testing.assert_array_equal(M.init,z)
    np.testing.assert_array_equal(C.ns,n0)
    np.testing.assert_array_equal(n0,b['snapshots'][0]['cold_bins'])
    for j,x,k in zip(M.edge_sets,a['hard_edge_sets'],b['hard_edge_sets'],strict=True):
        np.testing.assert_array_equal(j,x)
        np.testing.assert_array_equal(j,k)
    edges=np.asarray(M.edge_sets[2])
    limits=(1.4138888888888888,1.4666666666666666)
    cells=np.flatnonzero((edges[:-1]>=limits[0]-2e-14)&(edges[1:]<=limits[1]+2e-14))
    if len(cells)!=2*depth or np.any(np.diff(cells)!=1):
        raise ValueError('The fixed interval was not reproduced.')
    np.testing.assert_allclose(edges[[cells[0],cells[-1]+1]],limits,rtol=0.,atol=2e-14)
    j=M.r[2].start+cells
    moments=z[j]
    diff=y[j]-moments
    samples=np.linspace(0.,1.,9)
    r=[]
    for x in samples:
        state=z.copy()
        state[j]=y[j] if x==1. else moments+x*diff
        r.append(Cold.rhs(C,np.r_[state.ravel(),n0,0.,0.],False)[C.j0:-2])
    r=np.asarray(r)
    np.testing.assert_array_equal(r[0],a['snapshots'][0]['instantaneous_source_by_q'])
    if not np.all(np.isfinite(r)):
        raise ValueError('Non-finite source rate.')
    shift=r-r[0]
    curved=shift-samples[:,None]*shift[-1]
    mean=np.array([fsum(x) for x in shift])
    slopes=np.diff(mean)/(samples[1]-samples[0])
    return dict(scope='Nine convex local hard moment states over one fixed interval; frozen initialized RHS only.',
                depth=depth,cold_modes=len(n0),indices=j.tolist(),fractions=samples.tolist(),
                mean_shift=mean.tolist(),local_rate_by_mode=r.tolist(),shift_by_mode=shift.tolist(),
                deviation_from_linear_by_mode=curved.tolist(),
                every_mode_monotone_nonincreasing=bool(np.all(np.diff(r,axis=0)<=0.)),
                half_fraction_midpoint_relative_nonlinearity=float((mean[4]-.5*mean[-1])/abs(.5*mean[-1])),
                minimum_fraction_slope=float(min(slopes)),maximum_fraction_slope=float(max(slopes)))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('archives',nargs='+')
    a=p.parse_args()
    print(json.dumps([run(f) for f in a.archives],indent=2,allow_nan=False))
