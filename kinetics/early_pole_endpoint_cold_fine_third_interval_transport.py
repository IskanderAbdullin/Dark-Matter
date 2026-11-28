"""Audit the minimum q interval cost under both native D2 panel ledgers."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_endpoint_cold_non_nested_cdf_contrast import audit,enclosure
from early_pole_endpoint_cold_profile_interval_transport import fixed,transport


def run(fine,third,states,archive,fixture):
    a,b=[audit(p,s,archive,g) for p,s,g in zip((fine,third),states,('fine','third'),strict=True)]
    if [a[1],b[1]]!=fixture['states_sha256_uncompressed'] or [a[2],b[2]]!=fixture['gates']:
        raise ArithmeticError('The two charts do not bind the accepted nonnested response')
    by_mesh=[]
    for (e,A),(f,B) in zip(a[0],b[0],strict=True):
        k=np.union1d(e,f)
        mid=k[:-1]+np.diff(k)/2
        rows={}
        for label,i,j in (('forward_number',0,0),('inverse_number',0,1),('forward_excess',2,0)):
            u,v=A[:,i,j],B[:,i,j]
            lo,hi=[enclosure(q,x,mid) for q,x in ((e,u),(f,v))]
            gap=np.maximum(0.,np.maximum(lo[:,0]-hi[:,1],hi[:,0]-lo[:,1]))
            independent_cdf_bound=fsum(np.diff(k)*gap)
            P,Q=[fixed(dict(base=z/fsum(z))) for z in (u,v)]
            z=transport(P,Q,e,f,'fixed',independent_bracket=True)
            l,h=z['independent_cost_lower_upper_q']
            if h+2e-14<independent_cdf_bound or l<0.:
                raise ArithmeticError('Interval result contradicts CDF support implication')
            rows[label]=dict(continuous_CDF_pointwise_enclosure_lower_bound=independent_cdf_bound,
                            full_independent_interval_coupling=z,
                            full_normalized_finite_TV_lower_bound=fixture['observables'][label]['sharp_compatible_TV']['minimum_TV_bounds'][0])
        by_mesh.append(rows)
    a,b=by_mesh
    error=max(abs(a[k]['full_independent_interval_coupling']['minimum_interval_transport_cost_q']-
                  v['full_independent_interval_coupling']['minimum_interval_transport_cost_q']) for k,v in b.items())
    if error>6e-13:
        raise ArithmeticError('Pair did not independently agree on interval transport')
    return dict(grids=['fine','third'],x=.05,charge_subcells=2,full_band_q=[1e-6,1.5e-5],
                states_sha256_uncompressed=fixture['states_sha256_uncompressed'],
                observables=b,max_paired_mesh_cost_q_difference=error,
                scope='Arbitrary within-panel positions at unchanged native gross mass; not a physical momentum transfer.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('fine','third','states'):
        p.add_argument('--'+k,nargs=2,type=Path,required=True)
    for k in ('archive','fixture'):
        p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args()
    print(json.dumps(run(a.fine,a.third,a.states,a.archive,json.loads(a.fixture.read_text())),indent=2,allow_nan=False))
