"""Combine finite sign witnesses with one-arm-fixed numerical dual limits."""
import argparse
import gzip
import json
from math import fsum
from pathlib import Path


def compare(lo,hi,env):
    out=[]
    for a,b,g in zip(lo,hi,env['raw'][-1],strict=True):
        key=(a['hard_depth'],a['cold_modes'])
        if any((v['hard_depth'],v['cold_modes'])!=key for v in (b,g)):
            raise ValueError('Different physical fixture.')
        if (a['subdivision'],b['subdivision'])!=(12,24) or not g['fixed_gauge_moment']:
            raise ValueError('Expected retained moment and refinement controls.')
        if abs(g['independent_fitted_rate']-b['fitted_pair'])>1e-11:
            raise ValueError('Fitted initial-rate reference did not reproduce.')
        bands={}
        for moving,fixed,sign in (('native','local',1.),('local','native',-1.)):
            ds=[v['native' if moving=='native' else 'archive_local'] for v in g['cells']]
            anchor=b[fixed+'_fixed_influence']
            limits=sorted(sign*(fsum(v[z] for v in ds)-anchor) for z in ('dual_min','dual_max'))
            bands[moving+'_only']=dict(numerical_continuum_dual_range=limits,
                                       max_per_cell_dual_order_delta=max(v['max_dual_order_delta'] for v in ds))
        values={}
        for v,w in ((a,'12'),(b,'24')):
            values[w]={}
            for lab in ('both_arms','native_only','local_only'):
                p=v[lab]
                if p['status']=='verified_positive_witness':
                    if (p['positive_paired_rate']<=0. or p['max_normalized_moment_residual']>5e-13
                            or abs(p['verified_pointwise_amplitude']-p['solver_amplitude'])>1e-5):
                        raise ValueError('An admissibility witness did not pass.')
                elif p['status']!='infeasible_in_this_finite_approximation':
                    raise ValueError('Unexpected LP status.')
                values[w][lab]=p
        out.append(dict(hard_depth=key[0],cold_modes=key[1],one_arm_numerical_duals=bands,
                        finite_amplitude_witnesses=values))
    return dict(scope='Finite step-profile sign witnesses and non-rigorous continuous dual quadratures in the same frozen geometry.',
                fixtures=out)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('sub12')
    p.add_argument('sub24')
    p.add_argument('envelope')
    q=p.parse_args()
    with gzip.open(q.envelope,'rt') as f:
        print(json.dumps(compare(json.loads(Path(q.sub12).read_text()),
                                 json.loads(Path(q.sub24).read_text()),json.load(f)),indent=2,allow_nan=False))
