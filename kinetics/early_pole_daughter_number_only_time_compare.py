"""Quantify the finite effect of the accepted number-only complement retune."""
import argparse
import json
from pathlib import Path

import numpy as np


def compare(native, raw, number, both, static):
    labels=('native','local-daughter','number-only-retuned-local-daughter','retuned-local-daughter')
    models=(native,raw,number,both)
    s=static['cases']['number']
    c=number['complement_retune']
    for k in ('fixed_cells','complement_cells','match_axes','affine_logit_multipliers',
              'target_moments_relative_to_native','actual_moments_relative_to_native'):
        if s[k]!=c[k]:
            raise ValueError('Different static and evolved number-only retune.')
    h={}
    for key,x in zip(labels,models,strict=True):
        if x.get('prepared_charged_swap','native')!=key:
            raise ValueError('Unexpected preparation label.')
        for k in ('projection','duration','rtol','base_atol','history_exposure'):
            if x[k]!=native[k]:
                raise ValueError(f'Mismatched integration setup: {k}.')
        if not (abs(x['energy_defect'])<3e-12 and abs(x['entropy_defect'])<3e-11
                and abs(x['hard_number_defect'])<3e-12):
            raise ValueError('Original terminal gate failed.')
        a=np.array(x['cold_count_history'])
        if a.shape!=(x['report_panels']+1,) or a[0]!=native['cold_count_history'][0]:
            raise ValueError('Unexpected seed count or reporting grid.')
        h[key]=a-a[0]
    t=native['duration']
    relative={k:float((v[-1]-h['native'][-1])/t) for k,v in h.items()}
    return dict(mean_effect_from_native=relative,
                number_retune_minus_raw=float((h[labels[2]][-1]-h[labels[1]][-1])/t),
                two_moment_minus_number_retune=float((h[labels[3]][-1]-h[labels[2]][-1])/t),
                exposure=native['history_exposure'],
                cumulative_number_retune_minus_raw=(h[labels[2]]-h[labels[1]]).tolist(),
                static_number_fit=c,
                scope='Accepted manufactured preparations. The energy-only recipe failed its amplitude gate.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('native','raw','number','both','static'):
        p.add_argument(k)
    a=p.parse_args()
    print(json.dumps(compare(*(json.loads(Path(n).read_text()) for n in vars(a).values())),
                     indent=2,allow_nan=False))
