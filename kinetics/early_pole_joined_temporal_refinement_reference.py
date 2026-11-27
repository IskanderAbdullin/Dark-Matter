"""Compare nested macro grids on the same physical five-channel trajectory."""
import argparse
import json
import math

import numpy as np

from early_pole_joined_history_memory_reference import load


def conservation(run):
    energy,entropy=run['energy_drift'],run['entropy_closure']
    if not (abs(energy)<3e-12 and abs(entropy)<3e-11):
        raise ValueError('Only runs that pass the original terminal gates may be compared.')
    budget=run.get('elastic_energy_budget')
    drifts=[c['rho_drift'] for c in run['elastic']]
    peak=max(map(abs,drifts),default=0.)
    if budget is not None and (not np.isfinite(budget) or budget<=0. or peak>=budget):
        raise ValueError('Recorded step exceeds the stated elastic energy budget.')
    signed=math.fsum(drifts)
    return dict(terminal_energy_drift=energy,terminal_entropy_closure=entropy,
                elastic_energy_budget=budget,largest_step_energy_defect=peak,
                signed_elastic_energy_defect=signed,nonelastic_energy_remainder=energy-signed)


def compare(coarse, fine):
    controls = ('grid','pair','method','strength','projection','aperture_grid','duration',
                'hard_edge_sets','compton_strength','same_soft_strength','conversion_strength',
                'screening_update','response_smear')
    for k in controls:
        if coarse[k] != fine[k]:
            raise ValueError(f'Cannot compare different {k}.')
    for k,v in (('scalar_bins',8),('aperture_subdivision',1),('witness_mix',0.),
                ('response_radial_order',None)):
        if coarse.get(k,v)!=fine.get(k,v):
            raise ValueError(f'Cannot compare different {k}.')
    ledgers=[conservation(r) for r in (coarse,fine)]
    n,m = coarse['steps'],fine['steps']
    if not 0<n<m or m%n:
        raise ValueError('A strictly nested macro grid at the same duration is required.')
    ratio = m//n
    a,b = (np.array(r['state_history']) for r in (coarse,fine))
    if len(a)!=n+1 or len(b)!=m+1 or a[0]!=b[0]:
        raise ValueError('Expected exact same initial count and a complete history in both arms.')
    c,f = a-a[0], b[::ratio]-b[0]
    duration = coarse['duration']
    width = duration/n
    if not np.all(np.isfinite(c)) or not np.all(np.isfinite(f)):
        raise ValueError('Nonfinite count history.')
    if c[-1]!=coarse['scalar_gain'] or f[-1]!=fine['scalar_gain']:
        raise ValueError('Gain does not match the recorded complete history.')
    rc,rf = np.diff(c)/width,np.diff(f)/width
    diff = f-c
    sample = np.unique(np.r_[1,2,4,8,n//4,n//2,3*n//4,n])
    sample = sample[(sample>0)&(sample<=n)]
    rows = [dict(exposure=float(i*width),
                 coarse_gain=float(c[i]),fine_gain=float(f[i]),gain_change=float(diff[i]),
                 coarse_rate=float(rc[i-1]),fine_rate=float(rf[i-1]),
                 rate_change=float(rf[i-1]-rc[i-1])) for i in sample]
    sc = {r['step']:r for r in coarse.get('snapshots',())}
    sf = {r['step']:r for r in fine.get('snapshots',())}
    for fn in ('cold_bins','hard_moments'):
        if sc and sf:
            np.testing.assert_array_equal(sc[0][fn],sf[0][fn])
    snapshots = []
    for j in sorted(sc):
        if j*ratio not in sf:
            continue
        x,y = (np.asarray(s['instantaneous_source_by_q']) for s in (sc[j],sf[j*ratio]))
        np.testing.assert_allclose(sc[j]['exposure'],sf[ratio*j]['exposure'],rtol=2e-14,atol=0.)
        assert x.shape==y.shape
        gap = y-x
        snapshots.append(dict(exposure=sc[j]['exposure'],mode_sum_coarse=float(x.sum()),
                              mode_sum_fine=float(y.sum()),signed_gap=float(gap.sum()),
                              absolute_mode_sum_gap=float(np.abs(gap).sum())))
    return dict(grid=coarse['grid'],projection=coarse['projection'],duration=duration,coarse_steps=n,fine_steps=m,
                final_gain_coarse=float(c[-1]),final_gain_fine=float(f[-1]),final_gain_difference=float(diff[-1]),
                final_gain_relative_to_coarse=float(diff[-1]/abs(c[-1])) if c[-1] else None,
                maximum_absolute_history_difference=float(np.max(np.abs(diff))),
                maximum_absolute_aligned_interval_difference=float(np.max(np.abs(rf-rc))),
                minimum_coarse_and_aggregated_fine_interval_rates=[float(rc.min()),float(rf.min())],
                conservation=ledgers, aligned_intervals=rows,snapshot_instantaneous_differences=snapshots)


def fixture():
    key = dict(grid='fine',pair=True,method='backward',strength=1.,projection='steps',aperture_grid='fine',
               hard_edge_sets=[],duration=.25,compton_strength=1.,same_soft_strength=1.,
               conversion_strength=1.,screening_update='full',response_smear=.12,
               energy_drift=0.,entropy_closure=0.,elastic=[])
    coarse = dict(key,steps=2,state_history=[.25,.5,.75],scalar_gain=.5)
    fine = dict(key,steps=4,state_history=[.25,.375,.5625,.75,1.],scalar_gain=.75)
    d=compare(coarse,fine)
    assert d['final_gain_relative_to_coarse']==.5
    assert d['maximum_absolute_aligned_interval_difference']==1.5
    tests=(dict(projection='cell-conservative'),dict(energy_drift=3.1e-12),
           dict(elastic_energy_budget=5e-15,elastic=[dict(rho_drift=6e-15)]),
           dict(scalar_gain=0.))
    for t in tests:
        try:
            compare(coarse,dict(fine,**t))
        except ValueError:
            pass
        else:
            raise AssertionError('Mismatched control, energy or gain did not fail.')


def select(path,projection):
    d=load(path)
    if 'branches' in d:
        return d['branches'][projection]
    if 'evolution' in d:
        a,= (r for r in d['evolution'] if r['projection']==projection)
        return a
    if projection in d and isinstance(d[projection],dict):
        return d[projection]
    if d['projection']!=projection:
        raise ValueError('Raw trajectory projection does not match the requested one.')
    return d


if __name__=='__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('coarse')
    p.add_argument('fine')
    p.add_argument('--projection',choices=('steps','cell-conservative'),default='steps')
    a = p.parse_args()
    fixture()
    print(json.dumps(compare(select(a.coarse,a.projection),select(a.fine,a.projection)),indent=2,allow_nan=False))
