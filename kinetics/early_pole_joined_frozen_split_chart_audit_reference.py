"""Compare the admissible simultaneous chart with strictly held-screening split data."""
import argparse
import json
from pathlib import Path

import numpy as np


def compare(chart,split):
    if (chart['projection']!=split['projection'] or chart['duration']!=split['duration']
            or split['screening_update']!='frozen' or split['strength']!=1.
            or split['same_soft_strength']!=1. or not split['pair']
            or split['compton_strength']!=1. or split['conversion_strength']!=1.):
        raise ValueError('Same prepared input and held screening required.')
    for a,b in [(chart['energy_defect'],chart['entropy_defect']),
                (split['energy_drift'],split['entropy_closure'])]:
        if not abs(a)<3e-12 or not abs(b)<3e-11:
            raise ValueError('Both paths must pass their terminal conservation gates.')
    a=np.asarray(chart['cold_count_history'])
    b=np.asarray(split['state_history'])
    steps=split['steps']
    spans=chart['report_panels']
    if steps%spans or len(b)!=steps+1 or len(a)!=spans+1:
        raise ValueError('Use nested complete count and sample histories.')
    np.testing.assert_allclose(chart['history_exposure'],np.arange(spans+1)*chart['duration']/spans,
                               rtol=2e-14,atol=0.)
    b=b[::(steps//spans)]
    assert abs(a[0]-b[0])<3e-17
    diff=(a-a[0])-(b-b[0])
    va,vb=(q[-1]-q[0] for q in (a,b))
    assert abs(vb-split['scalar_gain'])<3e-21 and abs(va/chart['duration']-chart['mean_cold_rate'])<1e-15
    return dict(projection=chart['projection'],split_macro_steps=steps,chart_report_panels=spans,
                chart_calls=chart['calls'],mean_rates=[va/chart['duration'],vb/chart['duration']],
                integrated_chart_vs_split_relative=float(va/vb-1),
                maximum_absolute_integrated_difference=float(np.max(abs(diff))),
                scope='Finite fixed-scale accepted held-screening comparison. Chart accuracy must be checked separately.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('chart')
    p.add_argument('split')
    a=p.parse_args()
    print(json.dumps(compare(json.loads(Path(a.chart).read_text()),json.loads(Path(a.split).read_text())),
                     indent=2,allow_nan=False))
