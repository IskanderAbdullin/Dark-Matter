"""Validate a checked joined-model continuation against its uninterrupted control."""
import argparse
import gzip,json
from hashlib import sha256
from pathlib import Path

import numpy as np


def inspect(d):
    x,y,z=[d[k] for k in ('local-first','local-second','local-whole')]
    a,b,c=[t['checked_end_state'] for t in (x,y,z)]
    h=a['provenance']
    if (h!=b['provenance'] or h!=c['provenance'] or h.get('state_schema')!='joined-logit-f64-v1'
            or h.get('charged_swap')!='local-daughter' or len(h['module_sources'])<4):
        raise ValueError('Recorded joined-model checkpoints are not comparable.')
    for t in (x,y,z):
        if max(abs(t['energy_defect']),abs(t['hard_number_defect']))>3e-12 or abs(t['entropy_defect'])>3e-11:
            raise ArithmeticError('At least one measured history did not pass conservation.')
        if t['cold_momentum_nodes_over_T']!=x['cold_momentum_nodes_over_T']:
            raise ArithmeticError('The cold grid changed.')
    if not (a['time']==x['history_exposure'][-1]==y['history_exposure'][0]==z['history_exposure'][1]
            and b['time']==c['time']==y['history_exposure'][-1]==z['history_exposure'][-1]):
        raise ArithmeticError('The retained exposure did not resume exactly.')
    A,B,C=[np.asarray(t['cold_mode_count_history']) for t in (x,y,z)]
    np.testing.assert_array_equal(A[-1],B[0])
    np.testing.assert_array_equal(A[0],C[0])
    err=max(abs(B[-1]-C[-1]))
    native=max(abs(np.asarray(b['state'])-c['state']))
    if err>2e-17 or native>2e-8:
        raise ArithmeticError('The numerical segment restart did not recover the continuous control.')
    program=Path(__file__).with_name('early_pole_joined_simultaneous_chart_reference.py')
    digest=sha256(program.read_bytes()).hexdigest()
    return dict(scope='A small verified H1 local-daughter checkpoint comparison. No claim for long H8 restarts.',
                source_program_matches=bool(h['program']==digest),source_manifest_sha256=h['program'],
                fingerprinted_model_modules=len(h['module_sources']),fingerprinted_spectra=len(h['spectra']),
                continued_boundaries=[x['history_exposure'],y['history_exposure']],
                continuous_boundaries=z['history_exposure'],
                exact_cold_join=True,terminal_mode_defect=float(err),terminal_native_chart_defect=float(native),
                prefix_direct_discrepancy=float(max(abs(A[-1]-C[1]))),
                max_energy_number_entropy=[max(abs(t[k]) for t in (x,y,z))
                                          for k in ('energy_defect','hard_number_defect','entropy_defect')])


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('archive')
    a=p.parse_args()
    with gzip.open(a.archive,'rt') as f:d=json.load(f)
    print(json.dumps(inspect(d),indent=2,allow_nan=False))
