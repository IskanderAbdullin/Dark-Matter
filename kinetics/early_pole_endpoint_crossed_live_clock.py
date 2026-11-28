"""Cross held scalar spectra and charged states on exactly matched scalar grids."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_cell_endpoint_collision_pilot import CellChart
from early_pole_double_charge_stencil_pilot import event
from early_pole_endpoint_stencil_remap_pilot import SPECS,load


def case(d,y,x):
    m,E,_,H=d.state(x,y)
    q=[]
    for degree in (4,8):
        ch,defect=d.chart(x,m,degree)
        s,f,b=event(ch,y,d)
        left,mid,right=[fsum(r) for r in np.split(s,(d.ns,d.iP))]
        n=fsum(d.w[:d.ns]*y[:d.ns])
        g=[fsum(f),fsum(b)]
        if max(abs(left-mid),abs(left+right),abs(fsum(s*E)))>2e-12*max(g) or defect>2e-14:
            raise ArithmeticError('Crossed-state reaction violated paired moments')
        q.append(np.array([H*left/n,*[H*r/n for r in g]]))
    err=max(abs(q[0]/q[1]-1))
    if err>3e-6:
        raise ArithmeticError(('Crossed static clock still underresolved',float(err)))
    return dict(clock_net_forward_inverse=q[1].tolist(),scalar_mass=float(m),clock_prefactor=float(H),
                last_order_relative_drift=float(err))


def run(cards):
    x=cards[0]['x_end']; g=cards[0]['grid']; x0=cards[0]['reference_x']
    a=load(g)[0]
    Ds=[];states=[];labels=[]
    for c in cards:
        n=c['subcells_per_archived_cell']
        d=CellChart(a,SPECS[g],n,reference_x=x0)
        y=np.asarray(c['endpoint_occupation'])
        if c['grid']!=g or c['x_end']!=x or c['reference_x']!=x0 or c['counts']!=[d.ns,d.nf,d.np]:
            raise ArithmeticError('Cross only on the same grid and clock')
        if y.shape!=(d.n,) or not np.all(np.isfinite(y)) or min(y)<0 or max(y[d.ns:])>1:
            raise ArithmeticError('Invalid occupation source')
        Ds.append(d); states.append(y); labels.append(n)
    if len(Ds)<2 or labels!=sorted(set(labels)) or any(not np.array_equal(d.ds.ed,Ds[0].ds.ed) for d in Ds):
        raise ArithmeticError('Expected a strictly ordered charge comparison on one scalar grid')
    rows=[]
    for d,y in zip(Ds,states,strict=True):
        rows.append([case(d,np.r_[z[:d.ns],y[d.ns:]],x) for z in states])
    net=np.asarray([[z['clock_net_forward_inverse'][0] for z in r] for r in rows])
    comp=[]
    for i in range(1,len(rows)):
        j=i-1
        qe=net[i,j]-net[j,j]
        se=net[j,i]-net[j,j]
        ix=net[i,i]-net[i,j]-net[j,i]+net[j,j]
        comp.append(dict(from_to_subcells=[labels[j],labels[i]],total=float(net[i,i]-net[j,j]),
                         charge_chart_and_state_at_old_scalar=float(qe),scalar_at_old_charge=float(se),
                         interaction=float(ix)))
    return dict(grid=g,x=x,subcells=labels,rows_charge_chart_and_state_columns_scalar_state=rows,
                mixed_net_clock=net.tolist(),finite_difference_decompositions=comp,
                scope='Instantaneous frozen counterfactuals, not a coupled decomposition of the trajectory.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('state',nargs='+',type=Path)
    c=p.parse_args().state
    print(json.dumps(run([json.loads(f.read_text()) for f in c]),indent=2,allow_nan=False))
