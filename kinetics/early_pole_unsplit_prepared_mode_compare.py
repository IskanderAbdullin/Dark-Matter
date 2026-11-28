"""Audit cold mode histories in the same nominal unsplit prepared benchmark."""
import argparse
import json
from pathlib import Path

import numpy as np


KEYS=('native','local-daughter','local-complement-only')


def examine(mode_cases,scalar_cases,split_mode_receipt):
    base=mode_cases[0]
    samples=np.asarray(base['history_exposure'])
    momentum=np.asarray(base['cold_momentum_nodes_over_T'])
    if samples[0]!=0 or len(momentum)!=8 or np.any(np.diff(momentum)<=0):
        raise ValueError('Unexpected cold aperture history.')
    y={}
    replay={}
    for k,x,ref in zip(KEYS,mode_cases,scalar_cases,strict=True):
        if x.get('prepared_charged_swap','native')!=k or ref.get('prepared_charged_swap','native')!=k:
            raise ValueError('Wrong prepared input order.')
        for q in ('duration','projection','base_atol','rtol','report_panels','history_exposure'):
            if x[q]!=base[q] or ref[q]!=base[q]:
                raise ValueError('Mismatched numerical settings.')
        if x.get('complement_retune')!=ref.get('complement_retune'):
            raise ValueError('Compensation metadata changed.')
        if not (abs(x['energy_defect'])<3e-12 and abs(x['entropy_defect'])<3e-11
                and abs(x['hard_number_defect'])<3e-12):
            raise ValueError('Terminal chart ledger not accepted.')
        c=np.asarray(x['cold_mode_count_history'])
        if c.shape!=(len(samples),len(momentum)):
            raise ValueError('Missing a reported cold mode.')
        np.testing.assert_array_equal(x['cold_momentum_nodes_over_T'],momentum)
        np.testing.assert_allclose(c.sum(axis=1),x['cold_count_history'],rtol=0.,atol=1e-16)
        np.testing.assert_array_equal(x['cold_count_history'],ref['cold_count_history'])
        y[k]=c
        replay[k]=dict(energy_defect=x['energy_defect'],entropy_defect=x['entropy_defect'],
                       scalar_history_exact=True)
    n=y['native']
    out={}
    for k in KEYS[1:]:
        q=y[k]
        np.testing.assert_array_equal(q[0],n[0])
        v=q-n
        terminal=v[-1]/samples[-1]
        out[k]=dict(terminal_mean_effect_by_mode=terminal.tolist(),
                    terminal_negative_mode_indices=np.flatnonzero(terminal<0).tolist(),
                    terminal_positive_mode_indices=np.flatnonzero(terminal>0).tolist(),
                    cumulative_effect_by_mode=v.tolist(),
                    cumulative_negative_to_nonnegative_brackets=[
                      [[float(samples[i]),float(samples[i+1])] for i in range(1,len(v)-1)
                        if col[i]<0<=col[i+1]] for col in v.T])
    spl=np.asarray(split_mode_receipt['terminal_mode_effects']['outer64'])
    chart=np.asarray(out['local-daughter']['terminal_mean_effect_by_mode'])
    return dict(momentum_nodes_over_T=momentum.tolist(),exposure=samples.tolist(),cases=out,
                held64_minus_unsplit_local_terminal_by_mode=(spl-chart).tolist(),
                replay=replay,
                scope='Original 8-mode held chart; neither physical spectrum nor momentum-grid convergence.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('modes','scalar'):
        p.add_argument('--'+k,nargs=3,required=True,help='Native, local, then complement-only.')
    p.add_argument('--split-modes',required=True)
    a=p.parse_args()
    def read(q):
        return json.loads(Path(q).read_text())
    print(json.dumps(examine(list(map(read,a.modes)),list(map(read,a.scalar)),read(a.split_modes)),
                     indent=2,allow_nan=False))
