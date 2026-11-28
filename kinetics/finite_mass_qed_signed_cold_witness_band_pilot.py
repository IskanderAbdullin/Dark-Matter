"""Restartable finite frequency band, same refined parents and inner kernel."""
import argparse
import gzip
import json
import os
from pathlib import Path
import numpy as np
from scipy.special import roots_legendre
from finite_mass_qed_signed_cold_witness_reference import build, net, thermal_gate, Q
from finite_mass_qed_parent_integrated_signed_soft_reference import cap_knots


def panel_run(g,mesh,path,offset,maximum):
    root=Path(__file__).parent/'results'
    with gzip.open(root/f'early-pole-spectrum-{g}-m06-2026-11-26.json.gz') as stream:
        d=json.load(stream)
    n0=d['endpoint_scalar_occupations'][np.searchsorted(d['scalar_momentum_edges'],Q)-1]
    src=[build(g,d,x) for x in (0.,1.)]
    lo,hi=3e-6,1e-5
    bounds=[lo,*[x for x in cap_knots(src[0],hi) if lo<x<hi],hi]
    orders=((36,12,24,2),(64,18,36,3))
    dm=orders[mesh]
    x,w=roots_legendre(dm[-1])
    bas=[]
    for s in src:
        A,_,X=s.signed(0.,*dm[:3])[0]
        bas.append(X+A/n0)
    eq=thermal_gate(float(src[0].Es))
    old={}
    if path.exists():
        with open(path) as rd:
            for line in rd:
                row=json.loads(line)
                i=row['i']
                if (row['mesh']!=mesh or row['grid']!=g or row['orders']!=list(dm)
                        or row['bounds']!=bounds[i:i+2] or i in old):
                    raise ArithmeticError('Mismatched or duplicate saved panel')
                old[i]=row
    with open(path,'a',buffering=1) as wr:
        for i,(a,b) in enumerate(zip(bounds[:-1],bounds[1:],strict=True)):
            if i<offset or (maximum is not None and i>=offset+maximum) or i in old:
                continue
            Z=((a+b)/2+(b-a)*x/2).astype(np.longdouble)
            r=np.array([[(net(z,.9,n0,s,dm[:3],C)/z) for s,C in zip(src,bas,strict=True)] for z in Z])
            vals=(np.tensordot(w,r,axes=(0,0))*(b-a)/2).astype(float).tolist()
            print(json.dumps(dict(grid=g,mesh=mesh,orders=dm,i=i,bounds=[a,b],total_panels=len(bounds)-1,
                                  parts=vals,common_thermal_defect=eq),allow_nan=False),file=wr,flush=True)
            os.fsync(wr.fileno())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('grid',choices=('fine','third'))
    p.add_argument('mesh',type=int,choices=(0,1))
    p.add_argument('output',type=Path)
    p.add_argument('--offset',type=int,default=0)
    p.add_argument('--max-panels',type=int)
    o=p.parse_args()
    panel_run(o.grid,o.mesh,o.output,o.offset,o.max_panels)
