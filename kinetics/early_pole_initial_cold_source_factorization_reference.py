"""Separate static cold occupancy and dispersion changes at a common scalar seed."""
import argparse
import gzip
import json

import numpy as np

from early_pole_archived_shared_cold_window_reference import Cold
from early_pole_joined_scalar_exterior_source_reference import band_kernel


def run(bundle):
    with gzip.open(bundle,'rt') as f:
        raw=json.load(f)['branches']
    a,b=(raw[p] for p in ('steps','cell-conservative'))
    sa,sb=[x['snapshots'][0] for x in (a,b)]
    n=np.array(sa['cold_bins'])
    np.testing.assert_array_equal(n,sb['cold_bins'])
    for x,y in zip(a['hard_edge_sets'],b['hard_edge_sets'],strict=True):
        np.testing.assert_array_equal(x,y)
    C=Cold(a['grid'],1,0,len(n),True,exposure=a['duration'],
           reference_aperture_grid=a['aperture_grid'],archive_projection='steps')
    for x,y in zip(a['hard_edge_sets'],C.M.edge_sets,strict=True):
        np.testing.assert_array_equal(x,y)
    _,_,_,fs,pole=C.scalar(n)
    models=[C.M.weights(np.asarray(shot['hard_moments']),True)[:2] for shot in (sa,sb)]
    def evaluate(E,F):
        e=np.array([w*band_kernel(C,E,F,q,pole,v,12)
                    for q,w,v in zip(C.q,C.sw,fs,strict=True)])
        free=np.array([w*band_kernel(C,E,F,q,pole,0.,12)
                       for q,w in zip(C.q,C.sw,strict=True)])
        np.testing.assert_allclose(e[:,0]-e[:,1],e[:,2],rtol=2e-11,atol=3e-14)
        induced=e[:,2]-free[:,2]
        item=dict(total_net_rate=float(e[:,2].sum()),net_rate_each_mode=e[:,2].tolist(),
                  spontaneous=float(free[:,2].sum()),signed_induced_contrast=float(induced.sum()),
                  gross_gain=float(e[:,0].sum()),gross_inverse=float(e[:,1].sum()),
                  scalar_energy_source=float(e[:,3].sum()))
        return item
    results={}
    for k,j in [(0,0),(1,0),(0,1),(1,1)]:
        item=evaluate(models[k][0],models[j][1])
        tag=('continuous' if k else 'native')+'_dispersion_'+('continuous' if j else 'native')+'_occupations'
        results[tag]=item
    for tag,shot in (('native_dispersion_native_occupations',sa),('continuous_dispersion_continuous_occupations',sb)):
        np.testing.assert_allclose(results[tag]['net_rate_each_mode'],shot['instantaneous_source_by_q'],
                                   rtol=1e-8,atol=3e-12)
    swaps={}
    for name,i in [('parent_only',1),('daughter_only',2)]:
        F=models[0][1].copy()
        F[C.M.r[i]]=models[1][1][C.M.r[i]]
        swaps[name]=evaluate(models[0][0],F)
    y=[v['total_net_rate'] for v in results.values()]
    return dict(scope='Static diagnostic with fixed scalar seed and scale; unmatched hybrid hard energies/occupations are not evolved.',
                cold_grid=C.q.tolist(),scalar_pole=pole,scalar_occupation_range=[float(fs.min()),float(fs.max())],
                retained_aperture=a['aperture_grid'],cases=results,single_occupation_swaps_at_native_dispersion=swaps,
                held_occupation_dispersion_difference=y[1]-y[0],held_dispersion_occupation_difference=y[2]-y[0],
                mixed_rate_difference=y[3]-y[1]-y[2]+y[0],total_continuous_minus_native=y[3]-y[0])


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('bundle')
    a=p.parse_args()
    print(json.dumps(run(a.bundle),indent=2,allow_nan=False))
