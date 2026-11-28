"""Audit numerical stability of the matched unsplit prepared-input contrasts."""
import argparse
import json
from pathlib import Path

import numpy as np


def compare(nominal, tighter):
    keys=('native','local-daughter','outside-daughter')
    names=keys+('daughter',) if 'daughter' in nominal else keys
    if set(nominal)!=set(tighter) or set(nominal)!=set(names):
        raise ValueError('Require matched native, local and outside data; full daughter is optional.')
    metrics={}
    for key in names:
        a,b=nominal[key],tighter[key]
        if any(a[k]!=b[k] for k in ('projection','duration','report_panels','history_exposure')):
            raise ValueError(f'Mismatched tolerance inputs for {key}.')
        if a['rtol']!=3e-9 or b['rtol']!=3e-10 or a['base_atol']!=2e-14 or b['base_atol']!=2e-15:
            raise ValueError(f'Mismatched tolerance settings for {key}.')
        h=[]
        for d in (a,b):
            if d.get('prepared_charged_swap','native')!=key or not (
                    abs(d['energy_defect'])<3e-12 and abs(d['entropy_defect'])<3e-11
                    and abs(d['hard_number_defect'])<3e-12):
                raise ValueError(f'Preparation or conservation check failed: {key}.')
            h.append(np.asarray(d['cold_count_history']))
        np.testing.assert_equal(h[0][0],h[1][0])
        rates=np.array([(v[-1]-v[0])/a['duration'] for v in h])
        np.testing.assert_allclose(rates,[a['mean_cold_rate'],b['mean_cold_rate']],
                                   atol=0.,rtol=2e-15)
        metrics[key]=dict(means=rates.tolist(),
                          maximum_absolute_count_history_change=float(np.max(abs(h[1]-h[0]))),
                          rhs_calls=[a['calls'],b['calls']],
                          energy_defect=[a['energy_defect'],b['energy_defect']],
                          entropy_defect=[a['entropy_defect'],b['entropy_defect']])
    contrasts={}
    for key in names[1:]:
        ra=np.asarray(metrics['native']['means'])
        rb=np.asarray(metrics[key]['means'])
        d=rb-ra
        contrasts[key]=dict(contrast=d.tolist(),absolute_shift=float(d[1]-d[0]),
                            relative_change=float(d[1]/d[0]-1))
    factorial=None
    if 'daughter' in metrics:
        a,b,c,d=(np.asarray(metrics[k]['means']) for k in names)
        interactions=(d-a)-(b-a)-(c-a)
        factorial=dict(interaction=interactions.tolist(),
                       absolute_shift=float(interactions[1]-interactions[0]),
                       relative_change=float(interactions[1]/interactions[0]-1) if interactions[0] else None)
    return dict(cases=metrics,contrasts=contrasts,
                **(dict(factorial=factorial) if factorial is not None else {}),
                scope='Tolerance sensitivity of one fixed prepared-input benchmark, not a general accuracy bound.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('nominal','tighter'):
        p.add_argument('--'+k,nargs='+',required=True,help='Native, local, outside, optionally full (in order).')
    q=p.parse_args()
    def read(paths):
        keys=('native','local-daughter','outside-daughter','daughter')
        if len(paths) not in (3,4):
            raise ValueError('Require three or four inputs per profile.')
        return {k:json.loads(Path(f).read_text()) for k,f in zip(keys[:len(paths)],paths,strict=True)}
    print(json.dumps(compare(read(q.nominal),read(q.tighter)),indent=2,allow_nan=False))
