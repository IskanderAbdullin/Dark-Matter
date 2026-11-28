"""Cross independent hard-aperture and cold-output grids at the initial local swap."""
import argparse
import gzip
import json
from math import fsum

import numpy as np

from early_pole_archived_shared_cold_window_reference import Cold


def grids(path):
    with gzip.open(path,'rt') as f:
        d=json.load(f)['branches']
    a,b=[d[k] for k in ('steps','cell-conservative')]
    depth=a.get('aperture_subdivision',1)
    if depth!=b.get('aperture_subdivision',1) or depth not in (1,4):
        raise ValueError('Incorrect aperture hierarchy.')
    base,source=[np.asarray(v['snapshots'][0]['hard_moments']) for v in (a,b)]
    result={}
    archived=np.asarray(a['snapshots'][0]['cold_bins'])
    np.testing.assert_array_equal(archived,b['snapshots'][0]['cold_bins'])
    hold=None
    for n in (8,32):
        C=Cold(a['grid'],1,0,n,True,exposure=a['duration'],
               reference_aperture_grid=a['aperture_grid'],archive_projection='steps',
               aperture_subdivision=depth)
        M=C.M
        np.testing.assert_array_equal(M.init,base)
        for j,x,y in zip(M.edge_sets,a['hard_edge_sets'],b['hard_edge_sets'],strict=True):
            np.testing.assert_array_equal(j,x)
            np.testing.assert_array_equal(j,y)
        edges=np.asarray(M.edge_sets[2])
        selected=np.flatnonzero((edges[:-1]>=1.4138888888888888-2e-14)
                               &(edges[1:]<=1.4666666666666666+2e-14))
        if len(selected)!=2*depth or np.any(np.diff(selected)!=1):
            raise ValueError('Original daughter interval is not retained.')
        np.testing.assert_allclose(edges[[selected[0],selected[-1]+1]],(1.4138888888888888,1.4666666666666666),
                                   rtol=0.,atol=2e-14)
        ids=M.r[2].start+selected
        local=base.copy()
        local[ids]=source[ids]
        if n==8:
            np.testing.assert_array_equal(C.ns,archived)
            hold=local.copy()
        else:
            np.testing.assert_array_equal(hold,local)
        model_edges=np.geomspace(1e-6,1e-5,n+1)
        n0=C.ns
        grouped=n0.reshape(8,n//8).sum(axis=1)
        np.testing.assert_allclose(grouped,archived,rtol=0.,atol=4e-18)
        e0,_,_,_,m0=C.scalar(n0)
        if not np.isclose(m0,C.m0,rtol=3e-10,atol=1e-16):
            raise ValueError('The new seed did not initialize the same scalar pole.')
        d0,d1=[Cold.rhs(C,np.r_[v.ravel(),n0,0.,0.],False)[C.j0:-2] for v in (base,local)]
        if n==8:
            np.testing.assert_array_equal(d0,a['snapshots'][0]['instantaneous_source_by_q'])
        shift=d1-d0
        result[str(n)]=dict(native_rate=fsum(d0),local_rate=fsum(d1),paired_rate=fsum(shift),
                           mode_pair=shift.tolist(),common_group_pair=shift.reshape(8,n//8).sum(axis=1).tolist(),
                           grouped_seed_peak_difference=float(max(abs(grouped-archived))),
                           pole_over_T=float(m0),minimum_scalar_energy=float(min(e0)),edges=model_edges.tolist())
    return result


def cross(original,refined):
    if not np.array_equal(np.asarray(original['32']['edges'])[::4],original['8']['edges']):
        raise ValueError('Cold supports must form two exactly nested edge grids.')
    a,b,c,d=[v[k]['paired_rate'] for v in (original,refined) for k in ('8','32')]
    return dict(scope='Static initialized cold RHS only; cross chosen original-interval hard depth and retained cold mesh.',
                cases=dict(original=original,refined=refined),
                cold32_minus_8_original=b-a,cold32_minus_8_refined=d-c,
                hard4_minus_1_cold8=c-a,hard4_minus_1_cold32=d-b,
                crossed_difference=(d-c)-(b-a),
                cross_over_hard_change_on_cold8=(d-c-b+a)/(c-a))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--original',required=True)
    p.add_argument('--refined',required=True)
    o=p.parse_args()
    print(json.dumps(cross(grids(o.original),grids(o.refined)),indent=2,allow_nan=False))
