"""Build exact 32-cold-mode initial fixtures without changing archived hard inputs."""
import argparse
import gzip
import json
from pathlib import Path

import numpy as np

from early_pole_archived_shared_cold_window_reference import Cold


def rebin(path):
    with gzip.open(path,'rt') as f:
        data=json.load(f)['branches']
    output={}
    checks={}
    keys=('steps','cell-conservative')
    for label in keys:
        v=data[label]
        shot=v['snapshots'][0]
        x=np.asarray(shot['hard_moments'])
        cold=np.asarray(shot['cold_bins'])
        dep=v.get('aperture_subdivision',1)
        if dep not in (1,2,4) or cold.shape!=(8,) or v['projection']!=label:
            raise ValueError('Use an accepted 8-mode paired initializer on a validated hard mesh.')
        def make(n):
            model=Cold(v['grid'],1,0,n,True,exposure=v['duration'],
                       reference_aperture_grid=v['aperture_grid'],archive_projection=label,
                       aperture_subdivision=dep)
            np.testing.assert_array_equal(model.M.init,x)
            for e,g in zip(model.M.edge_sets,v['hard_edge_sets'],strict=True):
                np.testing.assert_array_equal(e,g)
            d=Cold.rhs(model,np.r_[x.ravel(),model.ns,0.,0.],False)[model.j0:-2]
            return model,d
        c0,d0=make(8)
        np.testing.assert_array_equal(c0.ns,cold)
        np.testing.assert_array_equal(d0,shot['instantaneous_source_by_q'])
        c,d=make(32)
        grouped=c.ns.reshape(8,4).sum(axis=-1)
        np.testing.assert_allclose(grouped,cold,atol=4e-18,rtol=0.)
        record={k:v[k] for k in ('duration','grid','hard_edge_sets','aperture_grid','projection','steps')}
        if dep!=1:
            record['aperture_subdivision']=dep
        initial=dict(shot,cold_bins=c.ns.tolist(),instantaneous_source_by_q=d.tolist())
        record['snapshots']=[initial]
        output[label]=record
        checks[label]=dict(hard_depth=dep,mode_count=len(d),largest_grouped_seed_difference=float(max(abs(grouped-cold))),
                           original_initial_rhs=float(d0.sum()),rebinned_initial_rhs=float(d.sum()),
                           provenance='Only the initial cold support and corresponding cold RHS are reconstructed.')
    u,v=(output[k] for k in keys)
    if u.get('aperture_subdivision',1)!=v.get('aperture_subdivision',1) or u['hard_edge_sets']!=v['hard_edge_sets']:
        raise ValueError('Mismatched paired hard mesh.')
    np.testing.assert_array_equal(u['snapshots'][0]['cold_bins'],v['snapshots'][0]['cold_bins'])
    return dict(branches=output,cold_rebin_checks=checks)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('input')
    p.add_argument('output')
    a=p.parse_args()
    v=rebin(a.input)
    Path(a.output).write_bytes(gzip.compress((json.dumps(v,separators=(',',':'),allow_nan=False)+'\n').encode(),mtime=0))
