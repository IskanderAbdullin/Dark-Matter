"""Compare coarse/fine pressure contrasts at one fixed final measurement mass."""
import argparse
import json
from pathlib import Path

import numpy as np

from early_pole_endpoint_scalar_subshell_spectrum_audit import read


def values(card,kind):
    if (card['x_start'],card['x_end'],card['steps'])!=(.04,.05,[1024,2048]):
        raise ArithmeticError('Expected paired last-segment clocks')
    if kind=='Gauss':
        x=np.asarray([card['coarse_count_and_pressures'],card['fine_count_and_pressures']])
        declared=np.asarray(card['coarse_over_fine_minus_one'])
        trace=card['maximum_trace_absolute_difference']
    else:
        if card['warm_div']!=kind:
            raise ArithmeticError('Wrong geometric warm chart')
        x=np.asarray(card['coarse_fine_count_pressure_and_warm_pressure'])
        declared=np.asarray(card['coarse_fine_relative_change'])
        trace=card['max_archived_trace_difference']
    if (x.shape!=(2,3) or declared.shape!=(3,) or min(x.ravel())<=0 or not np.all(np.isfinite(x))
            or not np.all(np.isfinite(declared)) or not 0<=trace<=5e-13):
        raise ArithmeticError('Invalid paired moments or trace replay')
    if np.max(abs(x[0]/x[1]-1-declared))>1e-15:
        raise ArithmeticError('Declared coarse-to-fine shift disagrees with saved rates')
    return x


def run(w1,w2,g):
    mass=w1['fixed_pressure_mass']
    if not np.isfinite(mass) or mass<=0 or any(card['fixed_pressure_mass']!=mass for card in (w2,g)):
        raise ArithmeticError('Not all moments share one fixed pressure mass')
    a,b,c=[values(r,k) for r,k in ((w1,1),(w2,2),(g,'Gauss'))]
    pairs=dict(W2_over_W1=b/a-1,Gauss_over_W1=c/a-1,Gauss_over_W2=c/b-1)
    return dict(fixed_mass=mass,columns=('number','pressure','warm_pressure'),steps=(1024,2048),
                individual_coarse_over_fine_minus_one=[(r[0]/r[1]-1).tolist() for r in (a,b,c)],
                same_clock_relative_contrasts={k:v.tolist() for k,v in pairs.items()},
                fine_minus_coarse_contrast={k:(v[1]-v[0]).tolist() for k,v in pairs.items()},
                scope='Three same-restart last-segment sensitivities, with each chart using its own archived .04 state. No full-history timestep error bound.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('w1',type=Path)
    p.add_argument('w2',type=Path)
    p.add_argument('gauss',type=Path)
    a=p.parse_args()
    print(json.dumps(run(*map(read,(a.w1,a.w2,a.gauss))),indent=2,allow_nan=False))
