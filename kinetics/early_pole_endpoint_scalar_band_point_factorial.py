"""Isolate the local cold response to scalar-bin refinements at one fixed charge chart."""
import argparse
import json
from pathlib import Path

import numpy as np

from early_pole_cell_endpoint_transverse_audit import charge_field,continuum_moments,model_callback,piecewise_logit,root
from early_pole_endpoint_cold_point_factorial import FACTORS,K,one,prepare
from early_pole_endpoint_scalar_subshell_chain_receipt import read
from early_pole_endpoint_scalar_subshell_pilot import ScalarSubshellChart
from early_pole_endpoint_stencil_remap_pilot import SPECS,X,load


def refined(c,n):
    if (c['grid'],c['x_end'],c['steps'],c['reference_x'],c['scalar_div'],c['charge_div'],c['cut'])!=(
            'fine',.05,2048,X,n,2,1e-5) or c.get('warm_floor') is not None:
        raise ArithmeticError('Expected a matched cold-only scalar state')
    d=ScalarSubshellChart(load('fine')[0],SPECS['fine'],2,n,1e-5)
    y=np.asarray(c['endpoint_occupation'])
    if c['counts']!=[d.ns,d.nf,d.np] or y.shape!=(d.n,) or not np.all(np.isfinite(y)) or min(y)<0 or max(y[d.ns:])>1:
        raise ArithmeticError('Invalid retained occupations')
    x=.05
    tbl=[]; ch=[]; cb=[]
    for mass,edges,occ in ((x,d.p_edges,d.wp*y[d.iP:]),(x*d.ratio,d.d_edges,d.wd*y[d.ns:d.iP])):
        p,f,_,_,a,b=continuum_moments(edges,mass,occ,36)
        tbl.append((p,f,mass))
        ch.append(piecewise_logit(mass,edges,a,b,occ,1))
        cb.append(model_callback(mass,edges,a,b))
    _,w,R=root(lambda beta:sum(charge_field(*v,beta,36) for v in tbl),K)
    return d.state(x,y)[0],(w,R),(ch,cb),d.ds.ed,y[:d.ns]


def run(data,levels=(1,2,4)):
    if len(data)!=3 or tuple(levels) not in ((1,2,4),(4,8,16)):
        raise ArithmeticError('Need lower or higher three consecutive scalar archives')
    (a,sa),(b,sb),(c,sc)=data
    if a['steps']!=2048:
        raise ArithmeticError('Unmatched native reference clock')
    chart=[prepare(a,2) if levels[0]==1 else refined(a,levels[0]),refined(b,levels[1]),refined(c,levels[2])]
    for (x,y) in zip(chart[:-1],chart[1:],strict=True):
        if any(not np.array_equal(p.E,q.E) for p,q in zip(x[2][0],y[2][0],strict=True)):
            raise ArithmeticError('The charged chart must remain fixed')
        if min(x[3])!=min(y[3]) or max(x[3])!=max(y[3]) or len(x[3])>=len(y[3]):
            raise ArithmeticError('Both scalar refinements must preserve the native support')
        j=np.searchsorted(y[3],x[3])
        lower=np.maximum(j-1,0)
        if np.any(np.minimum(abs(y[3][j]/x[3]-1),abs(y[3][lower]/x[3]-1))>5e-15):
            raise ArithmeticError('Scalar bins must form a strict nested refinement to roundoff')
    rows=[]
    for j,k in ((0,1),(1,2)):
        for q in (3.5e-6,4.7e-6):
            rows.append(dict(levels=[levels[j],levels[k]],**one([chart[j],chart[k]],q)))
    return dict(x=.05,k=K,grid='fine',charge_subcells=2,levels=levels,state_sha256=[sa,sb,sc],
                chart_factors=FACTORS,cols=('gross_forward','inverse'),
                scalar_gap_values=[v[0] for v in chart],photon_values=[v[1] for v in chart],rows=rows,
                scope='Four local cold factor swaps at fixed charge partition. Population, charged distribution and fitted poles can still vary on the evolved endpoints. No band-integrated attribution.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('states',type=Path,nargs=3)
    p.add_argument('--higher',action='store_true')
    a=p.parse_args()
    print(json.dumps(run([read(r) for r in a.states],(4,8,16) if a.higher else (1,2,4)),indent=2,allow_nan=False))
