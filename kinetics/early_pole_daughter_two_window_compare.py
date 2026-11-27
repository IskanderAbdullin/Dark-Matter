"""Compare two consecutive held windows without resetting occupations."""
import argparse
import gzip
import json
from pathlib import Path
from math import fsum

import numpy as np


def check(short,long,depth=1):
    if depth not in (1,2,4):
        raise ValueError("Use the accepted hard refinement depths.")
    t=np.asarray(short[0]['history_exposure'])
    n=len(t)-1
    finish=2*t[-1]
    t2=np.asarray(long[0]['history_exposure'])
    if len(t2)!=2*n+1 or t2[-1]!=finish or not np.array_equal(t,t2[:n+1]):
        raise ValueError('Two runs do not share the same nominal reported grid.')
    res=[]
    for a,b,k in zip(short,long,('native','local-daughter'),strict=True):
        scope=('projection','rtol','base_atol','cold_momentum_nodes_over_T')
        if any(a[f]!=b[f] for f in scope):
            raise ValueError('The time extension changed another control.')
        for v,panels,expected_times in ((a,n,t),(b,2*n,t2)):
            if (not np.array_equal(v['history_exposure'],expected_times) or v['report_panels']!=panels or v.get('prepared_charged_swap','native')!=k
                    or v.get('aperture_subdivision',1)!=depth or not v.get('elastic_channels',True)
                    or v.get('selected_elastic_channel','both')!='both'
                    or abs(v['energy_defect'])>=3e-12 or abs(v['entropy_defect'])>=3e-11
                    or abs(v['hard_number_defect'])>=3e-12):
                raise ValueError('Time data failed the original solver or scope gates.')
        x,y=[np.asarray(v['cold_mode_count_history']) for v in (a,b)]
        if x.shape!=(n+1,32) or y.shape!=(2*n+1,32) or not np.all(np.isfinite(y)&(y>0.)):
            raise ValueError('Invalid cold-mode history.')
        np.testing.assert_array_equal(x[0],y[0])
        if not np.all(np.isfinite(x)&(x>0.)):
            raise ValueError('Invalid short-control scalar history.')
        for data,raw in ((x,a),(y,b)):
            np.testing.assert_array_equal(data.T.copy(order='C').sum(axis=0),raw['cold_count_history'])
        err=float(np.max(abs(y[:n+1]-x)))
        if err>2e-18:
            raise ValueError('Extension did not recover the first short trajectory.')
        res.append(dict(max_absolute_early_reproduction_per_mode=err,gain=y-y[0]))
    modes=res[1]['gain']-res[0]['gain']
    d=np.array([fsum(row) for row in modes])
    u, v=(np.asarray(arm["cold_count_history"]) for arm in long)
    if np.max(abs(d-((v-v[0])-(u-u[0]))))>3e-17:
        raise ValueError("Paired compensated modes did not match recorded totals.")
    changes=(modes[n],modes[-1]-modes[n],modes[-1])
    return dict(scope='Two contiguous windows of the same held benchmark, without reinitialization or expansion.',
                exposure=t2.tolist(),paired_count=d.tolist(),
                negative_to_nonnegative_brackets=[t2[i:i+2].tolist() for i in np.flatnonzero((d[:-1]<0.)&(d[1:]>=0.))],
                mean_shifts=dict(first=float(d[n]/t[-1]),second=float((d[-1]-d[n])/t[-1]),whole=float(d[-1]/finish)),
                positive_modes=[int(np.count_nonzero(v>0.)) for v in changes],
                positive_paired_count_in_all_second_window_panels=bool(np.all(d[n+1:]>0.)),
                positive_panel_paired_gain_count_second_window=int(np.count_nonzero(np.diff(d[n:])>0.)),
                second_to_first_paired_rate=float((d[-1]-d[n])/d[n]),
                separately_reconstructed_energy_ledgers=[v['energy_defect'] for v in long],
                paired_mode_sum_by_interval=[v.tolist() for v in changes],
                shared_prefix_max_abs_error=[a['max_absolute_early_reproduction_per_mode'] for a in res])


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--short',nargs=2,required=True)
    p.add_argument('--long',nargs=2,required=True)
    p.add_argument('--hard-depth',type=int,choices=(1,2,4),default=1)
    a=p.parse_args()
    def r(f):
        if f.endswith('.gz'):
            with gzip.open(f,'rt') as h:
                return json.load(h)
        return json.loads(Path(f).read_text())
    print(json.dumps(check([r(f) for f in a.short],[r(f) for f in a.long],a.hard_depth),indent=2,allow_nan=False))
