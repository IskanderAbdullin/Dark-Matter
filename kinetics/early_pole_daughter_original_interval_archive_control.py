"""Validate exact local hard-moment replacement on both nested hard grids."""
import argparse
import gzip
import json

import numpy as np

from early_pole_archived_shared_cold_window_reference import Cold


def run(path):
    with gzip.open(path,'rt') as f:
        b=json.load(f)['branches']
    n,c=(b[k] for k in ('steps','cell-conservative'))
    dep=n.get('aperture_subdivision',1)
    if dep!=c.get('aperture_subdivision',1) or dep not in (1,4):
        raise ValueError('Expected a validated nested aperture depth.')
    y0,y1=(np.asarray(x['snapshots'][0]['hard_moments']) for x in (n,c))
    cold=np.asarray(n['snapshots'][0]['cold_bins'])
    np.testing.assert_array_equal(cold,c['snapshots'][0]['cold_bins'])
    C=Cold(n['grid'],1,0,len(cold),True,exposure=n['duration'],
           reference_aperture_grid=n['aperture_grid'],archive_projection='steps',
           aperture_subdivision=dep)
    M=C.M
    np.testing.assert_array_equal(M.init,y0)
    np.testing.assert_array_equal(C.ns,cold)
    for x,f,g in zip(M.edge_sets,n['hard_edge_sets'],c['hard_edge_sets'],strict=True):
        np.testing.assert_array_equal(x,f)
        np.testing.assert_array_equal(x,g)
    lo,hi=1.4138888888888888,1.4666666666666666
    edges=np.asarray(M.edge_sets[2])
    ids=np.flatnonzero((edges[:-1]>=lo-2e-14)&(edges[1:]<=hi+2e-14))
    if len(ids)!=2*dep or np.any(np.diff(ids)!=1):
        raise ValueError('Original physical support not partitioned as expected.')
    np.testing.assert_allclose([edges[ids[0]],edges[ids[-1]+1]],[lo,hi],rtol=0.,atol=2e-14)
    ix=M.r[2].start+ids
    y=y0.copy()
    y[ix]=y1[ix]
    untouched=np.setdiff1d(np.arange(len(y)),ix)
    np.testing.assert_array_equal(y[untouched],y0[untouched])
    rate=lambda z:Cold.rhs(C,np.r_[z.ravel(),cold,0.,0.],False)[C.j0:-2]
    base,local=rate(y0),rate(y)
    np.testing.assert_array_equal(base,n['snapshots'][0]['instantaneous_source_by_q'])
    return dict(aperture_subdivision=dep,indices=ix.tolist(),
                momentum_interval=[float(edges[ids[0]]),float(edges[ids[-1]+1])],
                native_initialized_cold_rate=float(base.sum()),
                local_initialized_cold_rate=float(local.sum()),
                local_minus_native_initialized_rate=float((local-base).sum()),
                local_minus_native_initialized_rate_by_mode=(local-base).tolist())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--original',required=True)
    p.add_argument('--refined',required=True)
    a=p.parse_args()
    print(json.dumps({k:run(v) for k,v in vars(a).items()},indent=2,allow_nan=False))
