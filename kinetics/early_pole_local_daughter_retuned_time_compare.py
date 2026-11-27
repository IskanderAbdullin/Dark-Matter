"""Compare the aperture-local prepared chart with and without fixed bare moments."""
import argparse
import json
from pathlib import Path

import numpy as np


def compare(native,local,retuned,initial):
    items={'native':native,'local-daughter':local,'retuned-local-daughter':retuned}
    static=initial['cases']['local']
    check=retuned['complement_retune']
    if any(check[k]!=static[k] for k in ('fixed_cells','complement_cells',
            'original_swap_moment_defect','number_and_bare_energy_defect','affine_logit_multipliers')):
        raise ValueError('The static and evolved retunes are not identical.')
    if max(map(abs,check['number_and_bare_energy_defect']))>=5e-14:
        raise ValueError('Unexpected two-moment residual.')
    raw={}
    for key,d in items.items():
        if any(d[k]!=native[k] for k in ('duration','projection','history_exposure','rtol','base_atol')):
            raise ValueError(f'Mismatched chart controls for {key}.')
        if d.get('prepared_charged_swap','native')!=key:
            raise ValueError(f'Mismatched prepared chart for {key}.')
        if not (abs(d['energy_defect'])<3e-12 and abs(d['entropy_defect'])<3e-11
                and abs(d['hard_number_defect'])<3e-12):
            raise ValueError(f'Terminal closure rejected for {key}.')
        h=np.array(d['cold_count_history'])
        if len(h)!=len(native['history_exposure']) or h[0]!=native['cold_count_history'][0]:
            raise ValueError('Inconsistent initial count or history width.')
        raw[key]=h-h[0]
    tm=native['duration']
    a,b,c=(raw[k] for k in items)
    old=(b[-1]-a[-1])/tm
    new=(c[-1]-a[-1])/tm
    reference=static['original_shift']
    return dict(initial_shifts=[reference,static['retuned_shift']],
                initial_relative_change=float(static['retuned_shift']/reference-1),
                finite_shifts=[float(old),float(new)],
                finite_mean_change=float(new-old),
                finite_relative_change=float(new/old-1) if old else None,
                cumulative_differential=(c-b).tolist(),exposure=native['history_exposure'],
                scope='Two bare daughter moments only, retuned outside the identical local aperture cells.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('native','local','retuned','initial'):
        p.add_argument(k)
    a=p.parse_args()
    print(json.dumps(compare(*(json.loads(Path(f).read_text()) for f in vars(a).values())),
                     indent=2,allow_nan=False))
