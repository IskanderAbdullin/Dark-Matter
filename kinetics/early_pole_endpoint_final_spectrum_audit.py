"""Compare final endpoint and original scalar histograms on the same scalar grid."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_cell_endpoint_collision_pilot import CellChart
from early_pole_endpoint_stencil_remap_pilot import SPECS,load


def summary(e,p,v):
    t=fsum(v)
    ans=[]
    for c in (1e-6,1e-5,3e-4,1e-3,1e-2):
        f=np.clip((c**3-e[:-1]**3)/np.diff(e**3),0,1)
        ans.append(float(fsum(f*v)/t))
    return dict(count_over_T_cubed=t,mean_node_q=fsum(p*v)/t,
                rms_node_q=np.sqrt(fsum(p*p*v)/t),
                cut_fractions=ans)


def run(card):
    g=card['grid'];a=load(g)[0]
    if card['x_end']!=a['x_end']:
        raise ArithmeticError('Endpoint comparisons require the same archived x')
    x=card['x_end']
    d=CellChart(a,SPECS[g],card['subcells_per_archived_cell'],reference_x=card['reference_x'])
    y=np.asarray(card['endpoint_occupation'])
    if card['counts']!=[d.ns,d.nf,d.np] or len(y)!=d.n or min(y)<0 or max(y[d.ns:])>1:
        raise ArithmeticError('Malformed endpoint snapshot')
    n=d.ds.w*y[:d.ns]
    r=np.array(a['endpoint_scalar_signed_shell_number_over_T3'])
    A=summary(d.ds.ed,d.ds.p,r)
    B=summary(d.ds.ed,d.ds.p,n)
    m=d.state(x,y)[0]
    M=a['checkpoints'][-1]['scalar_pole_over_T']
    TV=fsum(abs(n/fsum(n)-r/fsum(r)))/2
    return dict(grid=g,x=x,cuts=[1e-6,1e-5,3e-4,1e-3,1e-2],native=A,cell_endpoint=B,
                relative_total_count=fsum(n)/fsum(r)-1,normalized_scalar_distribution_total_variation=TV,
                original_scalar_pole_mass=M,new_scalar_pole_mass=m,relative_scalar_mass=m/M-1)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('state',type=Path)
    o=p.parse_args()
    print(json.dumps(run(json.loads(o.state.read_text())),indent=2,allow_nan=False))
