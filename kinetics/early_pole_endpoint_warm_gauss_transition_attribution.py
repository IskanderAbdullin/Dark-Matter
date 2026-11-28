"""Resolve paired Gauss pressure changes by conserved moments and parent cells."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_endpoint_gauss_warm_comparison import chart_for, gauss
from early_pole_endpoint_scalar_subshell_spectrum_audit import read
from early_pole_endpoint_scalar_warm_evolution_audit import unpack
from early_pole_endpoint_warm_pressure_cells import partition


def nodal_warm(card,m):
    order,d,repeats=chart_for(card)
    y=np.asarray(card['nodal_occupation'])
    w=d.ds.w*y[:d.ns]
    q=d.ds.p
    v=np.add.reduceat(w*q*q/np.hypot(q,m),np.r_[0,np.cumsum(repeats)[:-1]])
    return order,v,repeats>1


def run(base,cards,levels=(2,3,4)):
    levels=tuple(levels)
    if levels not in ((2,3,4),(4,5),(5,6)) or len(cards)!=len(levels):
        raise ArithmeticError('Need the accepted original or adjacent warm Gauss orders')
    x=base['x_end']
    d,y=unpack(base,False,x)
    m=d.state(x,y)[0]
    edges=d.ds.ed
    refined=edges[:-1]>=.01
    results=[]
    cell=[]
    for k,raw in zip(levels,cards,strict=True):
        z,env,ms=gauss(raw,d,m,x)
        j,mono,W=nodal_warm(raw,m)
        if j!=k or not np.array_equal(W,refined):
            raise ArithmeticError('Warm Gauss archive does not match proposed order and warm edges')
        direct=fsum(mono[W])+fsum(ms[~W,1])
        lo,hi=env['pressure_total_min_max']
        if not lo-9e-14*lo<=direct<=hi+9e-14*hi:
            raise ArithmeticError('Positive direct nodes leave their conditional pressure interval')
        results.append(dict(order=k,reconstructed=z.tolist(),hybrid_direct_warm_pressure=direct,
                            direct_over_reconstruction_minus_one=direct/z[1]-1))
        cell.append(ms)
    out=[]
    for i,(a,b) in enumerate(zip(cell[:-1],cell[1:],strict=True),levels[0]):
        P=fsum(a[:,1])
        Q=fsum(b[:,1])
        r=partition(a,b)/P
        agg=np.array([fsum(v) for v in r.T])
        if abs(agg[0]-(Q/P-1))>4e-15 or np.max(abs(r[~refined,2]))>5e-15:
            raise ArithmeticError('Paired pressure attribution does not close')
        net=fsum(r[:,0])
        gross=fsum(abs(r[:,0]))
        out.append(dict(lower_order=i,higher_order=i+1,total_count_change=fsum(b[:,0])/fsum(a[:,0])-1,
                    warm_count_change=fsum(b[refined,0])/fsum(a[refined,0])-1,
                    total_pressure_change_components=agg.tolist(),
                    warm_cells_pressure_change_components=[fsum(v) for v in r[refined].T],
                    sum_abs_cell_pressure_changes=gross, net_over_gross_abs=abs(net)/gross,
                    parent_q_edges=edges.tolist(),all_cell_pressure_change_components=r.tolist(),
                    largest_cell_changes=[dict(edges=edges[j:j+2].tolist(),parts=r[j].tolist())
                        for j in np.argsort(-abs(r[:,0]))[:5]]))
    return dict(x=x,fixed_mass=m,charts=results,transitions=out,
                components=['total','changed_cell_count','changed_pressure_per_particle'],
                scope='Differences between whole numerical trajectories, not independent causal interventions. Pressure in non-refined cells retains native constant convention, including any cell straddling .01.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('base',type=Path)
    p.add_argument('states',type=Path,nargs='+')
    p.add_argument('--levels',type=int,nargs='+',default=(2,3,4))
    a=p.parse_args()
    print(json.dumps(run(read(a.base),[read(p) for p in a.states],a.levels),indent=2,allow_nan=False))
