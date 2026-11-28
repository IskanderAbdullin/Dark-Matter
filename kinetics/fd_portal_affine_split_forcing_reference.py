"""Distinguish exact decay from maintenance of a constant affine forcing."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


def case(z):
    decay=np.exp(-z)
    cayley=(1.-z/2.)/(1.+z/2.)
    exact_source=-np.expm1(-z)
    exact_split_source=z*(1.+decay)/2.
    cayley_split_source=z*(1.+cayley)/2.
    e=exact_split_source/exact_source
    c=cayley_split_source/(1.-cayley)
    u=0.5*z/np.tanh(0.5*z)
    if abs(c-1.)>9e-15 or abs(e/u-1.)>4e-15 or e<=1.:
        raise ArithmeticError('Split affine forcing disagrees with independently regrouped formulas')
    return dict(z=z,homogeneous_exact=decay,homogeneous_midpoint=cayley,
                homogeneous_exact_split_steady_forcing_ratio=e,midpoint_split_steady_forcing_ratio=c)


def run(lock):
    b=lock.read_bytes()
    d=json.loads(b)
    result=[]
    for r in d['rows']:
        a,c=[x for x in r['steps'] if x['count']==8]
        if (a['method'],c['method'])!=('unfiltered','filtered') or a['initial_fast_rate_times_step']!=c['initial_fast_rate_times_step']:
            raise ArithmeticError('Require the exact same archived trajectory conditions')
        result.append(dict(T_SM=r['initial_visible_T'],linear_constant_forcing_model=case(a['initial_fast_rate_times_step']),
                           observed_nonlinear_eight_step_error_ratio=c['maximum_scaled_state_error']/a['maximum_scaled_state_error']))
    return dict(manufactured_samples=[case(z) for z in (.05,.25,1.,2.,4.,8.,20.,80.)],compared_probe=result,
                probe_sha256=hashlib.sha256(b).hexdigest(),
                scope='Exact scalar constant forcing identity suggests a mechanism; does not predict the nonlinear expanding endpoint errors quantitatively.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('probe',type=Path)
    print(json.dumps(run(p.parse_args().probe),indent=2,allow_nan=False))
