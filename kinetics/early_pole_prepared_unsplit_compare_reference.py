"""Isolate aperture-local and outside daughter effects in a common unsplit graph."""
import argparse
import json
from pathlib import Path

import numpy as np


def compare(native, local, outside):
    common = ('projection','duration','report_panels','rtol','base_atol','history_exposure')
    ref=np.asarray(native['cold_count_history'])
    out={}
    for key,x in (('native',native),('local-daughter',local),('outside-daughter',outside)):
        if any(native[k]!=x[k] for k in common):
            raise ValueError(f'Mismatched unsplit control or history grid: {key}.')
        if key!='native' and x.get('prepared_charged_swap')!=key:
            raise ValueError(f'Mismatched aperture swap: {key}.')
        history=np.asarray(x['cold_count_history'])
        if history.shape!=(x['report_panels']+1,) or history[0]!=ref[0]:
            raise ValueError(f'Mismatched seed count: {key}.')
        if not (abs(x['energy_defect'])<3e-12 and abs(x['entropy_defect'])<3e-11
                and abs(x['hard_number_defect'])<3e-12):
            raise ValueError(f'Failed closure: {key}.')
        rate=(history[-1]-history[0])/x['duration']
        np.testing.assert_allclose(x['mean_cold_rate'],rate,rtol=2e-15,atol=0.)
        shift=(history-ref)
        out[key]=dict(mean_cold_rate=float(rate),mean_effect_from_native=float(shift[-1]/x['duration']),
                      cumulative_effect=shift.tolist(),energy_defect=x['energy_defect'],
                      entropy_defect=x['entropy_defect'],calls=x['calls'])
    return dict(cases=out,history_exposure=native['history_exposure'],
                scope='Same frozen archived elastic graph, isolated prepared daughter input. No active screening feedback.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('native')
    p.add_argument('local_daughter')
    p.add_argument('outside_daughter')
    a=p.parse_args()
    print(json.dumps(compare(*(json.loads(Path(f).read_text()) for f in vars(a).values())),
                     indent=2,allow_nan=False))
