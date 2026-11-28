"""Separate local substitution and the matching outside retune on one held chart."""
import argparse
import json
from pathlib import Path

import numpy as np


def compare(native, local, outside, both):
    labels=('native','local-daughter','local-complement-only','retuned-local-daughter')
    arms=(native,local,outside,both)
    values={}
    controls=('projection','duration','report_panels','rtol','base_atol','history_exposure')
    common=outside['complement_retune']
    full=both['complement_retune']
    for k in ('fixed_cells','complement_cells','affine_logit_multipliers',
              'original_swap_moment_defect','number_and_bare_energy_defect'):
        if common[k]!=full[k]:
            raise ValueError('Require precisely the same retuning operation in both arms.')
    d0=np.asarray(full['original_swap_moment_defect'])
    d1=np.asarray(common['actual_complement_only_number_and_bare_energy_difference'])
    np.testing.assert_allclose(d0+d1,full['number_and_bare_energy_defect'],rtol=0.,atol=5e-14)
    for x,k in zip(arms,labels,strict=True):
        if x.get('prepared_charged_swap','native')!=k:
            raise ValueError(f'Swapped or unknown treatment: {k}.')
        if any(native[z]!=x[z] for z in controls):
            raise ValueError('Mismatched solver controls.')
        if not (abs(x['energy_defect'])<3e-12 and abs(x['entropy_defect'])<3e-11
                and abs(x['hard_number_defect'])<3e-12):
            raise ValueError('Prepared case failed the original terminal gates.')
        c=np.asarray(x['cold_count_history'])
        if c.shape!=(x['report_panels']+1,) or c[0]!=native['cold_count_history'][0]:
            raise ValueError('Mismatched seed or sample grid.')
        values[k]=c-c[0]
    a,b,c,d=(values[k] for k in labels)
    tm=native['duration']
    interaction=d-b-c+a
    return dict(mean={k:float(values[k][-1]/tm) for k in labels},
                local_effect=float((b[-1]-a[-1])/tm),
                compensation_only_effect=float((c[-1]-a[-1])/tm),
                combined_effect=float((d[-1]-a[-1])/tm),
                signed_interaction=float(interaction[-1]/tm),
                signed_interaction_history=interaction.tolist(),
                exposure=native['history_exposure'],
                scope='Finite interaction under one manufactured held-response retune; bare number and energy only.')


def compare_tolerances(nominal,tighter):
    for a,b in zip(nominal,tighter,strict=True):
        if a['rtol']!=3e-9 or b['rtol']!=3e-10 or a['base_atol']!=2e-14 or b['base_atol']!=2e-15:
            raise ValueError('Expected a tenfold tightening of both tolerances.')
        if any(a[k]!=b[k] for k in ('projection','duration','report_panels','history_exposure')):
            raise ValueError('The paired strict run changes more than tolerances.')
        if a.get('complement_retune')!=b.get('complement_retune'):
            raise ValueError('The complement retune changed with tolerance.')
    low,high=compare(*nominal),compare(*tighter)
    fields=('local_effect','compensation_only_effect','combined_effect','signed_interaction')
    return dict(nominal=low,tighter=high,
                signed_change={k:high[k]-low[k] for k in fields},
                largest_absolute_change_in_count_history=max(
                    float(np.max(abs(np.asarray(b['cold_count_history'])-a['cold_count_history'])))
                    for a,b in zip(nominal,tighter,strict=True)),
                scope='Tolerance sensitivity of the same held compensation, not a general error bound.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('native','local','outside','both'):
        p.add_argument(k)
    p.add_argument('--tighter',nargs=4,
                   help='Tighter native, local, compensation-only and combined raw receipts.')
    q=p.parse_args()
    def read(paths):
        return [json.loads(Path(f).read_text()) for f in paths]
    nominal=read([q.native,q.local,q.outside,q.both])
    print(json.dumps(compare_tolerances(nominal,read(q.tighter)) if q.tighter else compare(*nominal),
                     indent=2,allow_nan=False))
