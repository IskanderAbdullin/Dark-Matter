"""Compare occupation substitutions along the accepted first-window snapshots."""
import argparse
import gzip
import json

import numpy as np

from early_pole_archived_shared_cold_window_reference import Cold
from early_pole_joined_scalar_exterior_source_reference import band_kernel


def run(bundle):
    with gzip.open(bundle,'rt') as f:
        raw=json.load(f)['accepted_64']
    a,b=[raw[p] for p in ('steps','cell-conservative')]
    for field in ('grid','duration','steps','method','strength','screening_update','pair','same_soft_strength',
                  'compton_strength','conversion_strength','aperture_grid','hard_edge_sets'):
        if a[field]!=b[field]:
            raise ValueError(f'Unexpected branch mismatch: {field}')
    C=Cold(a['grid'],1,0,len(a['snapshots'][0]['cold_bins']),True,exposure=a['duration'],
           reference_aperture_grid=a['aperture_grid'],archive_projection='steps')
    for x,y in zip(C.M.edge_sets,a['hard_edge_sets'],strict=True):
        np.testing.assert_array_equal(x,y)

    def net(E,F,pole,fs):
        mode=np.array([v*band_kernel(C,E,F,q,pole,f,12)[2]
                       for q,v,f in zip(C.q,C.sw,fs,strict=True)])
        return mode

    rows=[]
    for sa,sb in zip(a['snapshots'],b['snapshots'],strict=True):
        assert sa['step']==sb['step']
        assert sa['exposure']==sb['exposure']
        ea,fa,*_=C.M.weights(np.asarray(sa['hard_moments']),True)
        eb,fb,*_=C.M.weights(np.asarray(sb['hard_moments']),True)
        _,_,_,f0,m0=C.scalar(np.asarray(sa['cold_bins']))
        _,_,_,f1,m1=C.scalar(np.asarray(sb['cold_bins']))
        base=net(ea,fa,m0,f0)
        end=net(eb,fb,m1,f1)
        for x,z in ((base,sa),(end,sb)):
            np.testing.assert_allclose(x,z['instantaneous_source_by_q'],rtol=1e-8,atol=3e-12)
        ret=np.zeros((2,3))
        highest=[]
        titles=['parent','daughter']
        for l,span in enumerate(C.M.r[1:]):
            x=fa.copy();x[span]=fb[span]
            only=net(ea,x,m0,f0)-base
            ret[l,0]=only.sum()
            edge=np.asarray(C.M.edge_sets[l+1])
            ids=np.flatnonzero((edge[:-1]>=1.4138888888888888-2e-14)&(edge[1:]<=1.4666666666666666+2e-14))
            assert len(ids)==2
            j=span.start+ids
            y=fa.copy();y[j]=fb[j]
            ret[l,1]=(net(ea,y,m0,f0)-base).sum()
            ret[l,2]=ret[l,0]-ret[l,1]
            highest.append(dict(species=titles[l],all_retained_swap=ret[l,0],original_two_cells_swap=ret[l,1],
                                outside_original_two_cell_difference=ret[l,2]))
        both=net(ea,fb,m0,f0)
        dispersion=net(eb,fb,m0,f0)-both
        seed=end-net(eb,fb,m0,f0)
        d=(end-base).sum()
        cross=(both-base).sum()-ret[:,0].sum()
        assert abs(ret[:,0].sum()+cross+dispersion.sum()+seed.sum()-d)<2e-15
        rows.append(dict(step=sa['step'],exposure=sa['exposure'],
                         native_rate=float(base.sum()),continuous_rate=float(end.sum()),
                         total_continuous_minus_native=float(d),occupation_swaps=highest,
                         occupation_cross_difference=float(cross),subsequent_dispersion_change=float(dispersion.sum()),
                         subsequent_cold_seed_change=float(seed.sum())))
    return dict(scope='Snapshot counterfactuals only; substitution order explicit, not separate evolved paths.',
                fine_steps=a['steps'],duration=a['duration'],rows=rows)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('bundle')
    a=p.parse_args()
    print(json.dumps(run(a.bundle),indent=2,allow_nan=False))
