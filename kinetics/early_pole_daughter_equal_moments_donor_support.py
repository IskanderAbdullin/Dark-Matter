"""Vary where the same two conserved bare moments are restored."""
import argparse
import gzip
import json

import numpy as np

from early_pole_archived_shared_cold_window_reference import Cold
from early_pole_daughter_complement_retune_reference import prepare


def choose_donors(full, fixed):
    below=full[full<fixed[0]]
    above=full[full>fixed[-1]]
    return dict(all_outside=np.r_[below,above],below=below,above=above,
                no_neighbors=np.r_[below[:-1],above[1:]])


def run(bundle):
    with gzip.open(bundle,'rt') as f:
        branches=json.load(f)['branches']
    a,b=(branches[k] for k in ('steps','cell-conservative'))
    u,v=(q['snapshots'][0] for q in (a,b))
    n=np.asarray(u['cold_bins'])
    y0,y1=(np.asarray(q['hard_moments']) for q in (u,v))
    np.testing.assert_array_equal(n,v['cold_bins'])
    C=Cold(a['grid'],1,0,len(n),True,exposure=a['duration'],
           reference_aperture_grid=a['aperture_grid'],archive_projection='steps')
    M=C.M
    np.testing.assert_array_equal(M.init,y0)
    for e,f,g in zip(M.edge_sets,a['hard_edge_sets'],b['hard_edge_sets'],strict=True):
        np.testing.assert_array_equal(e,f)
        np.testing.assert_array_equal(e,g)
    d=M.r[2]
    edge=np.asarray(M.edge_sets[2])
    ix=d.start+np.flatnonzero((edge[:-1]>=1.4138888888888888-2e-14)&(edge[1:]<=1.4666666666666666+2e-14))
    if len(ix)!=2:
        raise ValueError('Expected exactly two local daughter cells.')
    source=lambda z:Cold.rhs(C,np.r_[z.ravel(),n,0.,0.],False)[C.j0:-2]
    z0=source(y0)
    np.testing.assert_array_equal(z0,u['instantaneous_source_by_q'])
    out={}
    for k,c in choose_donors(np.arange(d.start,d.stop),ix).items():
        try:
            _,y,metadata=prepare(M,y0,y1,ix,donor_ids=None if k=='all_outside' else c)
        except ValueError as exc:
            out[k]=dict(rejected=str(exc),donor_cells=c.tolist())
            continue
        counts=source(y)-z0
        out[k]=dict(metadata,initial_mean_effect=float(counts.sum()),
                    initial_effect_by_mode=counts.tolist(),
                    max_cell_moment_shift=float(np.max(abs(y[c]-y0[c]))))
    return dict(cases=out,scope='Static, identical full daughter bare moments and local swap; donor supports vary.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('bundle')
    q=p.parse_args()
    print(json.dumps(run(q.bundle),indent=2,allow_nan=False))
