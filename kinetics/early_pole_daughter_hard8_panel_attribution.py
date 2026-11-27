"""Resolve the temporal location of the depth-four to depth-eight effect."""
import argparse
import gzip,json
from math import fsum
from pathlib import Path

import numpy as np


def raw(path):
    if path.endswith('.gz'):
        with gzip.open(path,'rt') as f:return json.load(f)
    return json.loads(Path(path).read_text())


def load_pair(name,depth,t,q,base_seed):
    s=[]
    for file,arm in zip(name,('native','local-daughter'),strict=True):
        x=raw(file)
        flags=(x.get('aperture_subdivision',1),x.get('prepared_charged_swap','native'),
               x.get('selected_elastic_channel','both'))
        if flags!=(depth,arm,'both') or not x.get('elastic_channels',True) or x['report_panels']!=16:
            raise ValueError('Unmatched arm or unsupported depth.')
        np.testing.assert_array_equal(x['history_exposure'],t)
        np.testing.assert_array_equal(x['cold_momentum_nodes_over_T'],q)
        a=np.array(x['cold_mode_count_history'])
        if a.shape!=(17,32):
            raise ValueError('Expected complete mode history.')
        np.testing.assert_array_equal(a[0],base_seed)
        if abs(x['energy_defect'])>=3e-12 or abs(x['entropy_defect'])>=3e-11 or abs(x['hard_number_defect'])>=3e-12:
            raise ValueError('Source ledger failed.')
        s.append(a-a[0])
    return s[1]-s[0]


def run(h4,h8,ladder):
    v=raw(h4[0])
    t=np.array(v['history_exposure'])
    q=np.array(v['cold_momentum_nodes_over_T'])
    anchor=np.array(v['cold_mode_count_history'])[0]
    M,N=[load_pair(s,d,t,q,anchor) for s,d in ((h4,4),(h8,8))]
    P=N-M
    panels=np.diff(P,axis=0)
    change=np.array([fsum(z) for z in P])
    x=np.diff(change)
    a,b=[np.array(next(s['paired_counts'] for s in ladder['levels'] if s['depth']==d)) for d in (4,8)]
    if max(abs(change-(b-a)))>1e-22 or max(abs(x-np.array([fsum(z) for z in panels])))>1e-22:
        raise ArithmeticError('Separately recomputed modal panel contributions do not close.')
    secs=[change[4*i+4]-change[4*i] for i in range(4)]
    turns=np.flatnonzero((x[:-1]<0)&(x[1:]>0))
    def split_part(k):
        return np.array([fsum(z) for z in k.reshape(4,8)])
    return dict(scope='Modal attribution of actual matched hard4 and hard8 first window only.',
                exposures=t.tolist(),cumulative_h8_minus_h4=change.tolist(),
                adjacent_panel_mean_rates=(x/np.diff(t)).tolist(),
                increasing_mode_count_by_panel=np.count_nonzero(panels>0,axis=-1).tolist(),
                sign_change_adjacent_panel_midpoints=[[(t[i]+t[i+1])/2,(t[i+1]+t[i+2])/2] for i in turns],
                first_quarter_fraction_of_terminal_difference=float(secs[0]/change[-1]),
                total_h8_minus_h4_by_quarter=secs,
                quartile_by_momentum=[split_part(P[i]).tolist() for i in (4,8,12,16)],
                mode_ratios_end_range=[float(np.min(N[-1]/M[-1])),float(np.max(N[-1]/M[-1]))])


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--hard4',nargs=2,required=True)
    p.add_argument('--hard8',nargs=2,required=True)
    p.add_argument('--ladder',required=True)
    a=p.parse_args()
    print(json.dumps(run(a.hard4,a.hard8,raw(a.ladder)),indent=2,allow_nan=False))
