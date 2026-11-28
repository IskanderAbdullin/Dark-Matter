"""Close the signed cold witness mesh and independent integration gates."""
import argparse
import json
from pathlib import Path

import numpy as np


def load(path):
    with open(path) as stream:
        return [json.loads(line) for line in stream if line.strip()]


def assemble(fine,third,reverse):
    rows=load(fine)+load(third)
    back=sorted(load(reverse),key=lambda r:r['k'])
    if len(rows)!=12 or len(back)!=2 or [x['k'] for x in back]!=[-1e-5,1e-5]:
        raise ArithmeticError('All signed witness and reversal controls are required')
    expected={(g,k,m) for g in ('fine','third') for k in (3e-6,1e-5) for m in (0,1,2)}
    maprow={(x['grid'],x['omega'],x['mesh']):x for x in rows}
    if set(maprow)!=expected:
        raise ArithmeticError('Unexpected or duplicate grid, frequency or mesh')
    if any(r['mix']!=[0.,.85,1.] or r['affine_defect']>1e-13 or r['thermal_pointwise_defect']>3e-13
           for r in rows):
        raise ArithmeticError('Bounded mixture linearity or thermal balance failed')
    maxfinal=0.
    for g,k in {(t[0],t[1]) for t in expected}:
        a,b=(np.array(maprow[g,k,m]['signed_vacuum_stimulus_total']) for m in (1,2))
        if np.any(np.sign(a)!=np.sign(b)):
            raise ArithmeticError(('Refinement did not fix the sign',g,k))
        defect=float(np.max(abs(a-b)))
        maxfinal=max(maxfinal,defect)
        if defect>1e-9 or np.any(abs(a[:,2]-b[:,2])>8e-5*np.abs(b[:,2])):
            raise ArithmeticError(('Pointwise frequency mesh check failed',g,k,defect))
    reverse_max=max(x['absolute_defect'] for x in back)
    if reverse_max>2e-13 or any(r['mix']!=.85 for r in back):
        raise ArithmeticError('Independent fixed-parent signed components not verified')
    return dict(scope='Pointwise signed vacuum soft data before dω/ω and gauge factor; no integrated rate.',
                final_angular_parent_mesh_absolute_defect=maxfinal,
                mixture_linearity_absolute_defect=max(r['affine_defect'] for r in rows),
                common_thermal_directional_relative_defect=max(r['thermal_pointwise_defect'] for r in rows),
                parent_first_stimulus_absolute_defect=reverse_max,reverse=back,
                rows=[maprow[z] for z in sorted(maprow)])


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('fine','third','reverse'):
        p.add_argument(n,type=Path)
    o=p.parse_args()
    print(json.dumps(assemble(o.fine,o.third,o.reverse),indent=2,allow_nan=False))
