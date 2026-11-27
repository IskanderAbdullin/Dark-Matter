"""Build static held fixtures at nested hard depths and compare raw/fitted rates."""
import argparse
import gzip
import json
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from early_pole_archived_shared_cold_window_reference import Cold
from early_pole_daughter_archive_pointwise_frozen_packets import compare


def build(template,depth,cold):
    base=template['branches']['steps']
    if depth not in (1,2,4,8,16) or cold not in (8,32):
        raise ValueError('Use the declared nested hard and cold supports.')
    d={}
    for k in ('steps','cell-conservative'):
        C=Cold(base['grid'],1,0,cold,True,exposure=base['duration'],
               reference_aperture_grid=base['aperture_grid'],archive_projection=k,
               aperture_subdivision=depth)
        p=C.M
        q=Cold.rhs(C,np.r_[p.init.ravel(),C.ns,0.,0.],False)[C.j0:-2]
        d[k]=dict(grid=base['grid'],aperture_grid=base['aperture_grid'],projection=k,duration=base['duration'],
                  hard_edge_sets=[list(v) for v in p.edge_sets],aperture_subdivision=depth,
                  snapshots=[dict(hard_moments=p.init.tolist(),cold_bins=C.ns.tolist(),
                                  instantaneous_source_by_q=q.tolist())])
    return dict(branches=d)


def run(original,refined):
    def load(path):
        with gzip.open(path,'rt') as f:
            return json.load(f)
    a,b=[load(f) for f in (original,refined)]
    result=[]
    green=[]
    # Validate exact old fixtures before constructing unarchived finer initial states.
    for dep,ref in ((1,a),(4,b)):
        z=build(a,dep,8)
        for k,c in z['branches'].items():
            q=ref['branches'][k]
            for f in ('grid','aperture_grid','projection','hard_edge_sets'):
                if c[f]!=q[f]:
                    raise ValueError('Did not regenerate the accepted geometry.')
            for f in ('hard_moments','cold_bins','instantaneous_source_by_q'):
                np.testing.assert_array_equal(c['snapshots'][0][f],q['snapshots'][0][f])
        green.append(dep)
    with TemporaryDirectory() as temp:
        p=Path(temp)/'initialized.json.gz'
        for dep in (1,2,4,8,16):
            for n in (8,32):
                with gzip.open(p,'wt') as f:
                    json.dump(build(a,dep,n),f,allow_nan=False)
                t=compare(p)
                terms=np.asarray(t['aperture_terms_by_mode_and_cell'])
                diff=terms[:,:,1]-terms[:,:,0]
                r={k:v for k,v in t.items() if k not in ('scope','aperture_terms_by_mode_and_cell','columns')}
                result.append(dict(r,abs_mode_cell_sum=float(abs(diff).sum()),
                                   max_absolute_mode_delta=float(max(abs(diff.sum(axis=1))))))
    return dict(scope='Initial frozen pointwise and fitted daughter responses only; positive depth refinements do not evolve new histories.',
                regenerated_accepted_depths_exact=green,ladder=result)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('original')
    p.add_argument('refined')
    q=p.parse_args()
    print(json.dumps(run(q.original,q.refined),indent=2,allow_nan=False))
