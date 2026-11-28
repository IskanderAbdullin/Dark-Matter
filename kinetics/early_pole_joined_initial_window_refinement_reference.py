"""Compare nested widths on the same early window of accepted joined paths."""
import argparse
import json

import numpy as np

from early_pole_joined_temporal_refinement_reference import select,conservation


def scores(runs):
    ledgers=[conservation(r) for r in runs]
    base=runs[0]
    width=base['duration']/base['steps']
    min_time=min(r['duration'] for r in runs)
    n=round(min_time/width)
    if n<1 or abs(min_time-n*width)>8*np.spacing(min_time):
        raise ValueError('Expected an integer number of coarse intervals on the common window.')
    initial={}
    paths=[]
    divisors=[]
    measures=('grid','projection','pair','method','strength','hard_edge_sets','aperture_grid',
              'compton_strength','same_soft_strength','conversion_strength','response_smear',
              'screening_update')
    canonical=('hard_moments','cold_bins','instantaneous_source_by_q')
    for r in runs:
        keys=('scalar_bins','aperture_subdivision','witness_mix','response_radial_order')
        defaults=(8,1,0.,None)
        for f,v in zip(keys,defaults,strict=True):
            if base.get(f,v)!=r.get(f,v):
                raise ValueError(f'Changed physical/discretization control: {f}')
        if any(base[m]!=r[m] for m in measures):
            raise ValueError('Cannot compare mismatched physical controls.')
        ratio=width/(r['duration']/r['steps'])
        div=round(ratio)
        if div<1 or not np.isclose(ratio,div,rtol=2e-14,atol=0.) or r['steps']<div*n:
            raise ValueError('Expected nested time widths and a complete common window.')
        if r['snapshots'][0]['exposure']!=0.:
            raise ValueError('Initial snapshot not identified.')
        for f in canonical:
            y=np.asarray(r['snapshots'][0][f])
            if f in initial:
                np.testing.assert_array_equal(initial[f],y)
            else:
                initial[f]=y
        seq=np.asarray(r['state_history'])
        if seq.size!=r['steps']+1 or seq[-1]-seq[0]!=r['scalar_gain']:
            raise ValueError('Incomplete integration history.')
        if seq[0]!=base['state_history'][0] or not np.all(np.isfinite(seq)):
            raise ValueError('Incorrect initial count or nonfinite evolution.')
        paths.append(seq[:n*div+1:div]-seq[0])
        divisors.append(div)
    if divisors!=sorted(set(divisors)):
        raise ValueError('Order histories by strictly increasing temporal resolution.')
    x=np.array(paths)
    rates=np.diff(x,axis=1)/width
    assert rates.shape==(len(runs),n)
    ref=rates[-1]
    return dict(grid=base['grid'],projection=base['projection'],common_window=min_time,
                resolution_divisors=divisors,initial_instantaneous_source=float(initial['instantaneous_source_by_q'].sum()),
                common_interval=width,first_interval_averages=rates[:,0].tolist(),
                scalar_gain_at_end=x[:,-1].tolist(),absolute_first_interval_offsets_from_fine=
                (rates[:,0]-ref[0]).tolist(),final_interval_averages=rates[:,-1].tolist(),
                maximum_interval_rate_offsets_from_fine=np.max(abs(rates-ref),axis=1).tolist(),
                conservation=ledgers)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('paths',nargs='+')
    p.add_argument('--projection',choices=('steps','cell-conservative'),default='steps')
    a=p.parse_args()
    print(json.dumps(scores([select(p,a.projection) for p in a.paths]),indent=2,allow_nan=False))
