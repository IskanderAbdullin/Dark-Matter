"""Locate when the paired cold-count gap accrues at fixed inner width."""
import argparse
import json
from pathlib import Path

import numpy as np

from early_pole_daughter_subcycle_width_compare import small_ledger
from early_pole_joined_temporal_refinement_reference import conservation


def read_profile(pair, n, h, ref):
    if len(pair)!=2:
        raise ValueError('Expected native then local.')
    a,b=pair
    controls=('projection','duration','grid','aperture_grid','hard_edge_sets',
              'strength','same_soft_strength','compton_strength','conversion_strength',
              'screening_update','method','pair','response_smear')
    if a['screening_update']!='frozen' or b['prepared_charged_swap']!='local-daughter':
        raise ValueError('Wrong screen control or preparation.')
    for v,w in zip(pair,ref,strict=True):
        if v['steps']!=n or any(v[k]!=w[k] for k in controls):
            raise ValueError('Mismatched held controls.')
        for k in ('hard_moments','cold_bins','instantaneous_source_by_q'):
            np.testing.assert_array_equal(v['snapshots'][0][k],w['snapshots'][0][k])
        if n<64:
            small_ledger(v)
        else:
            conservation(v)
            np.testing.assert_allclose([x['h'] for x in v['elastic']],h,
                                       rtol=1e-15,atol=0.)
    x,y=(np.asarray(v['state_history']) for v in pair)
    if x.shape!=y.shape or x.shape!=(n+1,) or x[0]!=y[0]:
        raise ValueError('Missing a count or cold start changed.')
    z={}
    for i,j in zip(a['snapshots'],b['snapshots'],strict=True):
        if i['exposure']!=j['exposure']:
            raise ValueError('Snapshot times differ.')
        z[i['exposure']]=(
            np.asarray(j['instantaneous_source_by_q'])-i['instantaneous_source_by_q'],
            np.asarray(j['hard_moments'])-i['hard_moments'])
    return ((y-y[0])-(x-x[0]))[::n//16],z


def compare(outer16,outer32,outer64):
    duration=outer64[0]['duration']
    t=np.linspace(0,duration,17)
    h=duration/64
    data={n:read_profile(p,n,h,outer64)
          for n,p in ((16,outer16),(32,outer32),(64,outer64))}
    shared=sorted(set.intersection(*(set(p[1]) for p in data.values())))
    if not shared or shared[0]!=0 or shared[-1]!=duration:
        raise ValueError('Missing common initial or final state sample.')
    out={}
    for n in (16,32):
        x,snap=data[n]
        y,reference=data[64]
        gap=x-y
        if not gap[-1]:
            raise ValueError('The terminal difference must be resolved.')
        rates=np.asarray([(snap[u][0]-reference[u][0]).sum() for u in shared])
        hard=np.asarray([np.max(abs(snap[u][1]-reference[u][1])) for u in shared])
        out[n]=dict(final_gap=float(gap[-1]),
                    first_coarse_interval_fraction=float(gap[1]/gap[-1]),
                    first_two_coarse_intervals_fraction=float(gap[2]/gap[-1]),
                    first_four_coarse_intervals_fraction=float(gap[4]/gap[-1]),
                    aligned_cumulative_gap=gap.tolist(),shared_snapshot_times=shared,
                    max_absolute_gap_of_paired_snapshot_rates=float(np.max(abs(rates))),
                    maximum_of_paired_hard_moment_component_gap=float(hard.max()),
                    paired_snapshot_rate_gap=rates.tolist(),paired_hard_moment_component_gap=hard.tolist())
    return dict(cases=out,aligned_exposure=t.tolist(),inner_width=h,
                scope='Recorded samples and integrated counts; no interpolation or claim about unobserved states.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('outer16','outer32','outer64'):
        p.add_argument('--'+k,nargs=2,required=True,help='Native then local.')
    a=p.parse_args()
    print(json.dumps(compare(*([json.loads(Path(f).read_text()) for f in paths]
                                 for paths in vars(a).values())),indent=2,allow_nan=False))
