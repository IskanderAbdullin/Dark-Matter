"""Match one total daughter moment while preserving the other at its raw value."""
import argparse
import gzip
import json

import numpy as np

from early_pole_archived_shared_cold_window_reference import Cold
from early_pole_daughter_complement_retune_reference import prepare


def run(bundle):
    with gzip.open(bundle,'rt') as f:
        branches=json.load(f)['branches']
    a,b=(branches[k] for k in ('steps','cell-conservative'))
    if a['projection']!='steps' or b['projection']!='cell-conservative':
        raise ValueError('Expected matched archived branches.')
    Y0,Y1=(np.asarray(x['snapshots'][0]['hard_moments']) for x in (a,b))
    N=np.asarray(a['snapshots'][0]['cold_bins'])
    np.testing.assert_array_equal(N,b['snapshots'][0]['cold_bins'])
    C=Cold(a['grid'],1,0,len(N),True,exposure=a['duration'],
           reference_aperture_grid=a['aperture_grid'],archive_projection='steps')
    M=C.M
    np.testing.assert_array_equal(M.init,Y0)
    for e,x,y in zip(M.edge_sets,a['hard_edge_sets'],b['hard_edge_sets'],strict=True):
        np.testing.assert_array_equal(e,x)
        np.testing.assert_array_equal(e,y)
    edges=np.asarray(M.edge_sets[2])
    loc=np.flatnonzero((edges[:-1]>=1.4138888888888888-2e-14)&(edges[1:]<=1.4666666666666666+2e-14))
    if len(loc)!=2:
        raise ValueError('Require the two original daughter aperture cells.')
    ix=M.r[2].start+loc
    source=lambda y:Cold.rhs(C,np.r_[y.ravel(),N,0.,0.],False)[C.j0:-2]
    ref=source(Y0)
    np.testing.assert_array_equal(ref,a['snapshots'][0]['instantaneous_source_by_q'])
    out={}
    for k,match in (('both',(True,True)),('number',(True,False)),('bare_energy',(False,True))):
        try:
            raw,fit,check=prepare(M,Y0,Y1,ix,match)
        except ValueError as exc:
            out[k]={'rejected':str(exc)}
            continue
        raw_shift=float((source(raw)-ref).sum())
        shifted=float((source(fit)-ref).sum())
        out[k]=dict(check,raw_shift=raw_shift,retuned_shift=shifted,
                   relative_change=float(shifted/raw_shift-1))
    return dict(cases=out,scope='Static controlled complement retunes, no time integration.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('bundle')
    x=p.parse_args()
    print(json.dumps(run(x.bundle),indent=2,allow_nan=False))
