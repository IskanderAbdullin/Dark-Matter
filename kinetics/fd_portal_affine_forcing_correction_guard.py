"""Exact affine correction cannot be an entropy-decreasing isolated kick."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


def case(z,steady):
    r=np.exp(-z)
    exact=-np.expm1(-z)*steady
    split=.5*z*(1+r)*steady
    c=exact-split
    amplitude=np.array([-2.,-.3,0.,.5,2.])
    add=.5*z*steady
    expected=r*amplitude+exact
    corrected=r*(amplitude+add)+c+add
    err=np.max(abs(corrected-expected))
    if err>2e-14*max(1.,abs(c)) or c*steady>=0.:
        raise ArithmeticError('Exact affine compensation failed')
    # If this offset is reassigned to the isolated middle collision, its input
    # can be exactly equilibrium. Then passive quadratic free energy increases.
    gain=.5*c*c
    if not gain>0. or abs(r*(steady+add)+c+add-steady)>2e-14*max(1.,abs(steady)):
        raise ArithmeticError('Isolated kick did not realize the predicted obstruction')
    return dict(z=z,steady=steady,compensation_offset=c,maximum_full_step_affine_error=err,
                isolated_equilibrium_free_energy_increase=gain)


def run(decay):
    b=decay.read_bytes()
    x=json.loads(b)
    grids=[(T['T_SM'],T['linear_constant_forcing_model']['z']) for T in x['compared_probe']]
    return dict(manufactured=[case(z,b) for z in (.05,.25,1.,2.,4.,8.,20.,80.) for b in (-.4,.7)],
                archived_step_parameters=[dict(T_SM=t,unit_steady_case=case(z,1.)) for t,z in grids],
                forcing_comparison_sha256=hashlib.sha256(b).hexdigest(),
                scope='Scalar affine numerical obstruction only. Do not reassign this offset to the isolated physical FD collision step.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('affine',type=Path)
    print(json.dumps(run(p.parse_args().affine),indent=2,allow_nan=False))
