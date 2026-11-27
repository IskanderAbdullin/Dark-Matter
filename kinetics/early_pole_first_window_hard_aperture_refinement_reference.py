"""Compare nested hard aperture supports and separately refine the initial packet rate."""
import argparse
import gzip
import json

import numpy as np

from early_pole_archived_shared_cold_window_reference import Cold
from early_pole_joined_temporal_refinement_reference import conservation,select


def static(bundle,projection,widths):
    with gzip.open(bundle,'rt') as f:
        old=json.load(f)['branches'][projection]
    n=len(old['snapshots'][0]['cold_bins'])
    out={}
    for sub in widths:
        c=Cold(old['grid'],1,0,n,True,exposure=old['duration'],
               reference_aperture_grid=old['aperture_grid'],archive_projection=projection,
               aperture_subdivision=sub)
        y=c.M.init
        np.testing.assert_array_equal(c.ns,old['snapshots'][0]['cold_bins'])
        if sub==1:
            np.testing.assert_array_equal(y,old['snapshots'][0]['hard_moments'])
        z=np.r_[y.ravel(),c.ns,0.,0.]
        d=Cold.rhs(c,z,False)
        cr=d[c.j0:-2]
        out[str(sub)]=dict(total_rate=float(cr.sum()),by_mode=cr.tolist(),
            retained_species_counts=[float(y[q,0].sum()) for q in c.M.r],
            retained_species_bare_energies=[float(np.sum(c.M.shape[q]*y[q])) for q in c.M.r],
            hard_cell_sizes=[len(e)-1 for e in c.M.edge_sets])
    return out


def compare(base,test):
    keys=('grid','projection','duration','steps','pair','method','strength','same_soft_strength',
          'compton_strength','conversion_strength','screening_update','response_smear','aperture_grid',
          'elastic_energy_budget')
    if any(base[k]!=test[k] for k in keys) or test.get('aperture_subdivision')!=2 or base.get('aperture_subdivision',1)!=1:
        raise ValueError('Compare the original and twice-subdivided aperture runs only.')
    a,b=(x['snapshots'][0] for x in (base,test))
    np.testing.assert_array_equal(a['cold_bins'],b['cold_bins'])
    np.testing.assert_allclose(a['charged_and_photon_counts'],b['charged_and_photon_counts'],atol=2e-14,rtol=0.)
    x,y=(p['hard_edge_sets'] for p in (base,test))
    for c,d in zip(x,y,strict=True):
        assert c[0]==d[0] and c[-1]==d[-1]
        assert set(c)<=set(d)
    new=[sorted(set(d)-set(c)) for c,d in zip(x,y,strict=True)]
    n0,n1=(r['scalar_gain'] for r in (base,test))
    return dict(extra_hard_anchors=new,mean_rates=[n0/base['duration'],n1/base['duration']],
                gain_change_fraction=n1/n0-1.,terminal_rates=[
                    sum(x['snapshots'][-1]['instantaneous_source_by_q']) for x in (base,test)],
                conservation=[conservation(x) for x in (base,test)])


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('bundle')
    p.add_argument('base')
    p.add_argument('refined')
    p.add_argument('--projection',choices=('steps','cell-conservative'),default='steps')
    a=p.parse_args()
    c,d=(select(p,a.projection) for p in (a.base,a.refined))
    k=compare(c,d)
    initial=static(a.bundle,a.projection,(1,2,4,8))
    for t,s in [('1',c),('2',d)]:
        np.testing.assert_allclose(initial[t]['by_mode'],s['snapshots'][0]['instantaneous_source_by_q'],
                                   atol=2e-14,rtol=1e-12)
    print(json.dumps(dict(grid=c['grid'],projection=a.projection,duration=c['duration'],
               finite=k,initial_only_depth=initial,
               scope='Additional hard cells only in the original aperture. Depths 4 and 8 are static, not trajectories.'),
               indent=2,allow_nan=False))
