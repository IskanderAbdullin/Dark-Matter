"""Check operator identity under actual held channel selectors at the same state."""
import argparse
import gzip
import json

import numpy as np

from early_pole_archived_cold_screened_elastic_reference import Screened
from early_pole_archived_shared_cold_window_reference import Cold
from early_pole_daughter_complement_retune_reference import prepare
from early_pole_screened_elastic_channel_sum_reference import ElasticChannels


def run(bundle):
    with gzip.open(bundle,'rt') as f:
        saved=json.load(f)['branches']
    a,b=(saved[k] for k in ('steps','cell-conservative'))
    q0,q1=(r['snapshots'][0] for r in (a,b))
    y0,y1=(np.asarray(x['hard_moments']) for x in (q0,q1))
    n=np.asarray(q0['cold_bins'])
    C=Screened(a['grid'],1,len(n),a['duration'],(1,1,2),.08,1.,a['aperture_grid'],'steps')
    M=C.M
    np.testing.assert_array_equal(y0,M.init)
    np.testing.assert_array_equal(n,C.ns)
    for e,f,g in zip(M.edge_sets,a['hard_edge_sets'],b['hard_edge_sets'],strict=True):
        np.testing.assert_array_equal(e,f)
        np.testing.assert_array_equal(e,g)
    edge=np.asarray(M.edge_sets[2])
    ix=M.r[2].start+np.flatnonzero((edge[:-1]>=1.4138888888888888-2e-14)&(edge[1:]<=1.4666666666666666+2e-14))
    if len(ix)!=2:
        raise ValueError('Expected original two-cell daughter preparation.')
    local,comp,_=prepare(M,y0,y1,ix)
    comp[ix]=y0[ix]
    both=ElasticChannels(C,a['grid'],1.)
    out={}
    for key,y in (('native',y0),('local',local),('compensation_only',comp)):
        _,F,G,*_=M.weights(y,True)
        z=np.r_[y.ravel(),n,0.,0.]
        C.strength=1.
        full=both.collision(G,F,G)
        unlike=C.collision(G,F,G)
        normal=Cold.rhs(C,z,True)
        try:
            C.strength=0.
            same=both.collision(G,F,G)
            switched=Cold.rhs(C,z,True)
        finally:
            C.strength=1.
        np.testing.assert_array_equal(normal,switched)
        e=full[0]-unlike[0]-same[0]
        scale=np.max(abs(full[0]))
        energy=float(np.sum(G*full[0]))
        result=dict(max_absolute_array_additivity_error=float(np.max(abs(e))),
                    relative_additivity_error=float(np.max(abs(e))/scale),
                    energy_pairing=energy,
                    entropy_additivity_error=float(full[1]-unlike[1]-same[1]),
                    work_additivity_error=float(full[2]-unlike[2]-same[2]),
                    nonelastic_bit_equal=True)
        if result['relative_additivity_error']>5e-15:
            raise ValueError('The chosen elastic vector fields are not additive at identical coordinates.')
        out[key]=result
    return dict(cases=out,scope='Static algebraic closure of channel switch, not a finite time superposition claim.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('bundle')
    q=p.parse_args()
    print(json.dumps(run(q.bundle),indent=2,allow_nan=False))
