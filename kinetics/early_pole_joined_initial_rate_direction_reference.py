"""Attribute the initial change of the cold source under the separately evaluated channels."""
import argparse
import gzip
import json

import numpy as np

from early_pole_archived_cold_screened_elastic_reference import Screened
from early_pole_archived_live_compton_reference import ComptonKernel
from early_pole_archived_live_flavor_conversion_reference import ConversionKernel
from early_pole_screened_active_response_update_reference import ActiveResponse
from early_pole_screened_elastic_channel_sum_reference import ElasticChannels
from early_pole_archived_shared_cold_window_reference import Cold


def run(bundle,projection):
    with gzip.open(bundle,'rt') as f:
        old=json.load(f)['branches'][projection]
    grid=old['grid']
    shot=old['snapshots'][0]
    N=np.array(shot['cold_bins'])
    C=Screened(grid,1,len(N),old['duration'],(1,1,2),.08,1.,old['aperture_grid'],projection)
    M=C.M
    Y=np.array(shot['hard_moments'])
    np.testing.assert_array_equal(Y,M.init)
    np.testing.assert_array_equal(N,C.ns)
    for a,b in zip(M.edge_sets,old['hard_edge_sets'],strict=True):
        np.testing.assert_array_equal(a,b)
    active=ActiveResponse(C,.12)
    active.refresh(Y,'full')
    el=ElasticChannels(C,grid,1.)
    E,F,G,*_=M.weights(Y,True)
    unlike=el.unlike.collision(G,F,G)[0]
    same=el.soft.collision(G,F,G)[0]
    D=Cold.rhs(C,np.r_[Y.ravel(),N,0.,0.],False)
    df=D[:C.j0].reshape(Y.shape)
    dn=D[C.j0:-2]
    P=M.pair_rhs(Y,True)[:-2].reshape(Y.shape)
    Co=ComptonKernel(M).collision(E,F,G)[0]
    Cv=ConversionKernel(M).collision(E,F,G)[0]
    Z=np.zeros_like(N)
    modes={'elastic_unlike':(unlike,Z),'elastic_same':(same,Z),
           'decay_reciprocal':(df,dn),'pair':(P,Z),'compton':(Co,Z),'conversion':(Cv,Z)}
    modes['total']=(sum(v[0] for v in modes.values()),dn)

    def source(y,n):
        return Cold.rhs(C,np.r_[y.ravel(),n,0.,0.],False)[C.j0:-2]

    b=source(Y,N)
    np.testing.assert_allclose(b,np.array(shot['instantaneous_source_by_q']),rtol=1e-12,atol=1e-14)
    agg=[]
    for t in (1e-14,3e-15,1e-15):
        item={}
        for k,(u,v) in modes.items():
            am=source(Y-t*u,N-t*v)
            ap=source(Y+t*u,N+t*v)
            s=(ap-am)/(2*t)
            item[k]=dict(mode_directional_derivative=s.tolist(),total_directional_derivative=float(s.sum()),
                         max_abs_hard_update=float(np.max(abs(t*u))),max_abs_cold_update=float(np.max(abs(t*v))))
        ref=np.array(item['total']['mode_directional_derivative'])
        terms=sum(np.array(item[k]['mode_directional_derivative']) for k in item if k!='total')
        agg.append(dict(probe_exposure=t,channels=item,
                        relative_additivity_defect=float(np.max(abs(ref-terms))/np.max(abs(ref)))))
    return dict(projection=projection,initial_exposure=shot['exposure'],source=float(b.sum()),directions=agg,
                scope='Tangent from the prepared input at time zero, not a finite trajectory or fitted rate.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('bundle')
    p.add_argument('--projection',choices=('steps','cell-conservative'),default='steps')
    a=p.parse_args()
    print(json.dumps(run(**vars(a)),indent=2,allow_nan=False))
