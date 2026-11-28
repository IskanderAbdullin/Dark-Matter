"""Frozen native-occupation fixture for scalar-only nested shells."""
import json
from math import fsum

import numpy as np

from early_pole_cell_endpoint_collision_pilot import initial_checks,CellChart
from early_pole_endpoint_scalar_subshell_pilot import ScalarSubshellChart
from early_pole_cell_endpoint_startup_pilot import drive
from early_pole_endpoint_stencil_remap_pilot import SPECS,X,load


def run(g):
    a=load(g)[0]
    baseline=CellChart(a,SPECS[g],2,reference_x=X)
    unit=ScalarSubshellChart(a,SPECS[g],2,1,1e-5)
    frozen=lambda d: drive(d,.05,32,x_start=.049997,initial=d.initial)
    ref,y=frozen(baseline)
    alt,z=frozen(unit)
    if ref!=alt or not np.array_equal(y,z):
        raise ArithmeticError('Unit scalar refinement is not identical')
    rows=[]
    for n in (1,2,4):
        d=ScalarSubshellChart(a,SPECS[g],2,n,1e-5)
        if any(not np.array_equal(getattr(baseline,k),getattr(d,k))
               for k in ('d_edges','p_edges','wd','wp','qd','qp')):
            raise ArithmeticError('Scalar refinement changed the charge chart')
        diag=initial_checks(d,.05)
        m=d.state(.05,d.initial)[0]
        rows.append(dict(scalar_div=n,scalar_count=d.ns,scalar_number=fsum(d.w[:d.ns]*d.initial[:d.ns]),
                         scalar_mass=m,**diag))
    if any(r['scalar_number']!=rows[0]['scalar_number'] for r in rows):
        raise ArithmeticError('Scalar refinement lost total number')
    return dict(grid=g,charge_div=2,unit_short_trajectory_matches=True,x=.05,
                exact_native_cut=1e-5,cards=rows)


if __name__=='__main__':
    print(json.dumps([run(g) for g in ('fine','third')],indent=2,allow_nan=False))
