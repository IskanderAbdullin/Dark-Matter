"""Compare final scalar shape across charged refinements with the scalar grid fixed."""
import argparse
import json
from pathlib import Path

import numpy as np

from early_pole_cell_endpoint_collision_pilot import CellChart
from early_pole_endpoint_scalar_spectral_moment_budget import CUTS,LABELS,integrate
from early_pole_endpoint_scalar_shape_transport import CUT,cdf,pair
from early_pole_endpoint_scalar_subshell_spectrum_audit import read
from early_pole_endpoint_stencil_remap_pilot import SPECS,X,load


def run(paths):
    cards=[read(p) for p in paths]
    g=cards[0]['grid']
    native=load(g)[0]
    levels=[c['subcells_per_archived_cell'] for c in cards]
    if any(b!=2*a for a,b in zip(levels[:-1],levels[1:],strict=True)):
        raise ArithmeticError('Expected nested consecutive charge charts')
    batches=[]
    for c,n in zip(cards,levels,strict=True):
        d=CellChart(native,SPECS[g],n,reference_x=X)
        y=np.asarray(c['endpoint_occupation'])
        if (c['grid'],c['reference_x'],c['x_end'],c['steps'],c['counts'])!=(g,X,.05,2048,[d.ns,d.nf,d.np]) or (
            y.shape!=(d.n,) or not np.all(np.isfinite(y)) or min(y)<0 or max(y[d.ns:])>1):
            raise ArithmeticError('Expected same archived branch, scalar grid and terminal clock')
        batches.append((d,y,d.state(.05,y)[0]))
    m=batches[0][-1]
    out=[]
    for (a,v,u),(b,z,w),i,j in zip(batches[:-1],batches[1:],levels[:-1],levels[1:],strict=True):
        if not np.array_equal(a.ds.ed,b.ds.ed):
            raise ArithmeticError('Scalar support changed in the charge-only ladder')
        A,B=[integrate(d,y,m,48) for d,y in ((a,v),(b,z))]
        totals=[d.sum(axis=0) for d in (A,B)]
        out.append(dict(subcells=[i,j],cold=pair(cdf(a,v,CUT),cdf(b,z,CUT)),
                        full=pair(cdf(a,v,a.ds.ed[-1]),cdf(b,z,b.ds.ed[-1])),
                        total_moment_relative_changes=(totals[1]/totals[0]-1).tolist(),
                        above_point01_moment_relative_changes=(B[-2:].sum(axis=0)/A[-2:].sum(axis=0)-1).tolist(),
                        gap_relative_change=float(w/u-1)))
    return dict(grid=g,x=.05,fixed_moment_mass=float(m),moments=LABELS,cuts=CUTS,rows=out,
                scope='Charge refinement sensitivity on archived scalar histograms; not a continuum or cosmological prediction.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('state',type=Path,nargs='+')
    a=p.parse_args()
    if len(a.state)<2:
        p.error('Need two or more consecutive charge resolution states')
    print(json.dumps(run(a.state),indent=2,allow_nan=False))
