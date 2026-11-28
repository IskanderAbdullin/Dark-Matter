"""Audit adjacent positive-node warm source changes on a fixed S4 state."""
import argparse
import json
from math import fsum
from pathlib import Path
import numpy as np

from early_pole_endpoint_warm_positive_gauss_ladder import WarmGauss
from early_pole_endpoint_scalar_subshell_spectrum_audit import read
from early_pole_endpoint_scalar_warm_evolution_audit import unpack
from early_pole_endpoint_stencil_remap_pilot import load
from early_pole_endpoint_warm_source_cells import local


def run(card,ref,through=5):
    if through not in (5,6):
        raise ArithmeticError('Use the original or one higher positive quadrature')
    x=card['x_end']
    a,y=unpack(card,False,x)
    m,_,_,H=a.state(x,y)
    ed=a.ds.ed
    floor=ref['warm_floor']
    warm=ed[:-1]>=floor
    e=ed[np.searchsorted(ed,floor):]
    if (x,m,H)!=(ref['x'],ref['mass_held'],ref['clock_held']) or e.tolist()!=ref['warm_native_edges']:
        raise ArithmeticError('Warm source does not match its frozen state/cells')
    V={r['rule']:np.asarray(r['warm_cell_forward_inverse_net']) for r in ref['cases']}
    R=V['Volume64']
    N=fsum(R[:,2])
    if N<=0 or abs(N/ref['ref_net_pressure']-1)>2e-14 or ref['mesh_reference']['global_L1_fraction']>4e-5:
        raise ArithmeticError('Expected positive paired volume reference')
    output=[]
    raw=[]
    for p in range(2,through+1):
        b=WarmGauss(load('fine')[0],p)
        z=b.transfer_step(y)
        q=b.ds.p
        W=np.where(q>=floor,q*q/np.hypot(q,m),0.)
        if any(not np.array_equal(getattr(a,k),getattr(b,k)) for k in ('wp','wd','qp','qd')):
            raise ArithmeticError('Changed charge grid')
        count=np.add.reduceat(b.ds.w*z[:b.ns],np.r_[0,np.cumsum(b.repeats)[:-1]])
        n=a.ds.w*y[:a.ns]
        if np.max(abs(count-n)/np.maximum(n,1e-105))>3e-14:
            raise ArithmeticError('Native count changed')
        lo,hi=[local(b,z,W,x,m,H,ed,k)[warm] for k in (4,8)]
        defect=float(np.max(abs(lo-hi)/np.maximum(hi[:,0,None]+hi[:,1,None],1e-100)))
        if defect>3e-6 or (p==2 and np.max(abs(V['Gauss2']-hi)/np.maximum(hi[:,0,None]+hi[:,1,None],1e-100))>3e-12):
            raise ArithmeticError('Frozen source lacks a resolved mesh or known order control')
        raw.append(hi)
        P=(hi[:,2]-R[:,2])/N
        output.append(dict(order=p, scalar_nodes=b.ns, aggregate_residual=fsum(P),
                           l1_signed_source_residual=fsum(abs(P)), paired_radial_defect=defect))
    changes=[]
    for i,(a,b) in enumerate(zip(raw[:-1],raw[1:],strict=True),2):
        h=(b-a)/N
        j=np.argmax(abs(h[:,2]))
        changes.append(dict(orders=[i,i+1],total_relative_components=[fsum(v) for v in h.T],
                            l1_local_net_change=fsum(abs(h[:,2])),
                            leading_parent=e[j:j+2].tolist(),components_at_leading_parent=h[j].tolist(),
                            all_warm_parent_changes=h.tolist()))
    return dict(x=x, archived_warm_edges=e.tolist(),reference_pressure_rate=N,
                reference_relative_resolution=ref['mesh_reference']['global_L1_fraction'],
                comparison_columns=['forward','inverse','net'],
                charts=output,changes=changes,
                scope='Held S4 scalar distribution, optics and charge: instantaneous finite source only, not a measured evolution or its causal decomposition.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('state',type=Path)
    p.add_argument('reference',type=Path)
    p.add_argument('--through',type=int,choices=(5,6),default=5)
    a=p.parse_args()
    print(json.dumps(run(read(a.state),read(a.reference),a.through),indent=2,allow_nan=False))
