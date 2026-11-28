"""Compare scalar-only refinements on the same native shell support."""
import argparse
import gzip
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_endpoint_scalar_subshell_pilot import ScalarSubshellChart
from early_pole_endpoint_stencil_remap_pilot import SPECS,X,load
from early_pole_resolved_parent_reference import S


def read(p):
    op=gzip.open if p.suffix=='.gz' else open
    with op(p,'rt',encoding='utf-8') as f:
        return json.load(f)


def compare(cards):
    old=cards[0]
    g=old['grid']
    qcut=cards[1]['cut']
    x=old['x_end']
    if old['subcells_per_archived_cell']!=2 or old['reference_x']!=X:
        raise ArithmeticError('Expected the two-charge-subcell native control')
    a=load(g)[0]
    native_edges=np.asarray(a['scalar_momentum_edges'])
    levels=[1]+[c['scalar_div'] for c in cards[1:]]
    if any(q!=2*p for p,q in zip(levels[:-1],levels[1:],strict=True)):
        raise ArithmeticError('Expected successive nested scalar doubling')
    records=[]
    for n,c in enumerate(cards):
        s=1 if not n else c['scalar_div']
        if c['grid']!=g or c['x_end']!=x or c['reference_x']!=X or (
            n and (c['charge_div']!=2 or c['cut']!=qcut or c['steps']!=old['steps'])):
            raise ArithmeticError('Mismatched common archive clock, charge chart, or time order')
        d=ScalarSubshellChart(a,SPECS[g],2,s,qcut)
        if c['counts']!=[d.ns,d.nf,d.np]:
            raise ArithmeticError('Wrong shell discretization')
        y=np.asarray(c['endpoint_occupation'])
        if y.shape!=(d.n,) or not np.all(np.isfinite(y)) or min(y)<0 or max(y[d.ns:])>1:
            raise ArithmeticError('Invalid occupations')
        m,E,*_=d.state(x,y)
        S0=d.w[:d.ns]*y[:d.ns]
        M=np.bincount(d.native_scalar_labels,weights=S0)
        N=fsum(S0)
        if len(M)!=len(native_edges)-1 or not N>0:
            raise ArithmeticError('Native support or positive counts are missing')
        cut=[1e-6,1e-5,.0003]
        frac=lambda u:[fsum(S0[d.ds.ed[1:]<=v])/N for v in u]
        unitE=fsum(S0*E[:d.ns])/N
        excess=fsum(S0*d.ds.p*d.ds.p/(E[:d.ns]+m))/N
        if abs((m+excess)/unitE-1)>2e-14 or abs(fsum(M)/N-1)>2e-14:
            raise ArithmeticError('Scalar measure did not close')
        mass_cdf=np.cumsum(S0/N)
        qtile=[]
        # Compare both the common native shell counts and the explicitly
        # piecewise-constant fitted distribution inside refined shells.
        for p in (.1,.5,.9,.999):
            i=int(np.searchsorted(mass_cdf,p))
            p0=0. if i==0 else mass_cdf[i-1]
            f=(p-p0)/(mass_cdf[i]-p0)
            ql3,qh3=d.ds.ed[i:i+2]**3
            qtile.append(float(np.cbrt(ql3+f*(qh3-ql3))))
        records.append((dict(scalar_div=s,scalar_bins=d.ns,total_yield=N/S,scalar_mass=float(m),
             node_rms_q=float(np.sqrt(fsum(S0*d.ds.p**2)/N)),mean_excess=excess,
             cut_fractions=frac(cut),q_quantiles_from_piecewise_constant_cells=qtile),M/N))
    out=[]
    for i,(r,p) in enumerate(records):
        if i:
            old,p0=records[i-1]
            cdf=np.cumsum(p-p0)
            j=np.argmax(abs(cdf))
            r.update(yield_change_from_previous=r['total_yield']/old['total_yield']-1,
                     pole_change_from_previous=r['scalar_mass']/old['scalar_mass']-1,
                     native_coarsened_total_variation=float(fsum(abs(p-p0))/2),
                     native_max_signed_cdf_difference=float(cdf[j]),
                     native_max_cdf_difference_upper_q=float(native_edges[j+1]))
        out.append(r)
    return dict(grid=g,x=x,cuts=cut,quantile_levels=[.1,.5,.9,.999],scope='One matched charge chart. Coarsening TV does not measure intrabin shape.',rows=out)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('control',type=Path)
    p.add_argument('refinement',type=Path,nargs='+')
    a=p.parse_args()
    print(json.dumps(compare([read(f) for f in [a.control,*a.refinement]]),indent=2,allow_nan=False))
