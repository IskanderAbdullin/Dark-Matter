"""Extend the locked warm source quadrature from G8 to G9/G10/G12."""
import argparse
import hashlib
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_endpoint_scalar_subshell_chain_receipt import read
from early_pole_endpoint_scalar_warm_evolution_audit import unpack
from early_pole_endpoint_stencil_remap_pilot import load
from early_pole_endpoint_warm_positive_gauss_ladder import WarmGauss
from early_pole_endpoint_warm_source_cells import local


def one(state,ref,prior):
    card,digest=read(state)
    b,old=unpack(card,False,card['x_end'])
    x=card['x_end']
    m,_,_,H=b.state(x,old)
    rb=ref.read_bytes()
    pb=prior.read_bytes()
    r,p=json.loads(rb),json.loads(pb)
    ed=b.ds.ed
    floor=r['warm_floor']
    warm=ed[:-1]>=floor
    edges=ed[np.searchsorted(ed,floor):]
    A=np.asarray(next(d for d in r['cases'] if d['rule']=='Volume64')['warm_cell_forward_inverse_net'])
    N=fsum(A[:,2])
    if ((x,m,H)!=(r['x'],r['mass_held'],r['clock_held']) or x!=p['x']
            or edges.tolist()!=r['warm_native_edges'] or p['archived_warm_edges']!=edges.tolist()
            or abs(N/r['ref_net_pressure']-1.)>3e-14 or N!=p['reference_pressure_rate']
            or r['mesh_reference']['global_L1_fraction']>4e-5):
        raise ArithmeticError('Reference optics, chart, or source changed')
    data=[]
    keep=[]
    for k in (8,9,10,12):
        c=WarmGauss(load('fine')[0],k)
        y=c.transfer_step(old)
        z=c.ds.p
        W=np.where(z>=floor,z*z/np.hypot(z,m),0.)
        if any(not np.array_equal(getattr(b,g),getattr(c,g)) for g in ('wp','wd','qp','qd')):
            raise ArithmeticError('Held charge geometry changed')
        n=b.ds.w*old[:b.ns]
        nn=np.add.reduceat(c.ds.w*y[:c.ns],np.r_[0,np.cumsum(c.repeats)[:-1]])
        if max(abs(nn-n)/np.maximum(n,1e-105))>3e-14:
            raise ArithmeticError('Held native count changed')
        lo,hi=[local(c,y,W,x,m,H,ed,mesh)[warm] for mesh in (4,8)]
        scale=np.maximum(hi[:,0,None]+hi[:,1,None],1e-100)
        defect=np.max(abs(lo-hi)/scale)
        if defect>3e-6:
            raise ArithmeticError('Independent radial rule discrepancy')
        d=(hi[:,2]-A[:,2])/N
        q=float(fsum(d)); l=float(fsum(abs(d)))
        if k==8:
            exact=p['charts'][-1]
            if [(exact['order'],exact['scalar_nodes']),(q,l)]!=[(k,c.ns),(exact['aggregate_residual'],exact['l1_signed_source_residual'])]:
                raise ArithmeticError('Archived G8 was not exactly reproduced')
        data.append(dict(order=k,nodes=c.ns,polynomial_exactness_error=c.polynomial_error,
                         aggregate_residual=q,l1_net_cell_residual=l,paired_radial_defect=defect))
        keep.append(hi)
    tab=[]
    for i,j,y,z in zip((8,9,10),(9,10,12),keep[:-1],keep[1:],strict=True):
        a=(z-y)/N
        h=np.argmax(abs(a[:,2]))
        tab.append(dict(orders=[i,j],total_relative_components=[fsum(q) for q in a.T],
                        normalized_l1_local_net_drift=fsum(abs(a[:,2])),
                        largest_parent_q_bounds=edges[h:h+2].tolist(),largest_parent_drift=a[h].tolist()))
    return dict(x=x,state_sha256=digest,reference_and_prior_sha256=[hashlib.sha256(v).hexdigest() for v in (rb,pb)],
                orders=data,adjacent_changes=tab)


def run(states,references,priors):
    rows=[one(s,r,p) for s,r,p in zip(states,references,priors,strict=True)]
    if [a['x'] for a in rows]!=[.01,.05]:
        raise ArithmeticError('Require the two declared and separately archived early/terminal held states')
    return dict(rows=rows,scope='Frozen S4 warm pressure source on held states only. Neither the compared Volume64 nor G12 is assumed exact or an evolved history.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for x in ('states','references','priors'):
        p.add_argument('--'+x,type=Path,nargs=2,required=True)
    a=p.parse_args()
    print(json.dumps(run(a.states,a.references,a.priors),indent=2,allow_nan=False))
