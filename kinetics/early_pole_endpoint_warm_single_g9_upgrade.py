"""Directly test a one-node warm budget transfer on a separate held endpoint."""
import argparse
import hashlib
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_endpoint_scalar_subshell_chain_receipt import read
from early_pole_endpoint_scalar_warm_evolution_audit import unpack
from early_pole_endpoint_warm_budgeted_gauss_probe import Mixed
from early_pole_endpoint_warm_source_cells import local


def one(state,ref,accepted,idx,other):
    card,digest=read(state)
    if digest!=accepted['state_sha256']:
        raise ArithmeticError('Held baseline bytes do not match the archived order diagnostic')
    obj,y=unpack(card,False,card['x_end'])
    x=card['x_end']
    m,_,_,H=obj.state(x,y)
    raw=ref.read_bytes()
    r=json.loads(raw)
    ed=obj.ds.ed
    vol=r['warm_floor']
    v=ed[:-1]>=vol
    Q0=np.asarray(next(k for k in r['cases'] if k['rule']=='Volume64')['warm_cell_forward_inverse_net'])
    scale=r['ref_net_pressure']
    if (x,m,H)!=(r['x'],r['mass_held'],r['clock_held']) or hashlib.sha256(raw).hexdigest()!=accepted['reference_and_prior_sha256'][0]:
        raise ArithmeticError('Frozen physical and numerical parameters changed')

    def rule(orders):
        q=Mixed(orders)
        t=q.lift(y)
        if any(not np.array_equal(getattr(obj,k),getattr(q,k)) for k in ('wp','wd','qp','qd')):
            raise ArithmeticError('Charged layout changed in single-node test')
        n=obj.ds.w*y[:obj.ns]
        z=np.add.reduceat(q.ds.w*t[:q.ns],np.r_[0,np.cumsum(q.repeats)[:-1]])
        if np.max(abs(z-n)/np.maximum(n,1e-105))>3e-14:
            raise ArithmeticError('Number in the original scalar cells was not held fixed')
        w=np.where(q.ds.p>=vol,q.ds.p*q.ds.p/np.hypot(q.ds.p,m),0.)
        a,b=[local(q,t,w,x,m,H,ed,j)[v] for j in (4,8)]
        err=np.max(abs(a-b)/np.maximum(b[:,0,None]+b[:,1,None],1e-100))
        if err>3e-6:
            raise ArithmeticError('Assigned source rule did not resolve with two radial orders')
        return q.ns,b

    n=int(v.sum())
    n8,q8=rule([8]*n)
    n9,q9=rule([9]*n)
    d=(q9-q8)/scale
    if max(abs(np.array([fsum(z) for z in d.T])-accepted['adjacent_changes'][0]['total_relative_components']))>3e-12:
        raise ArithmeticError('Direct independently assembled G8→G9 differences did not close')
    ans=[]
    for j in (idx,other):
        p=[8]*n
        p[j]=9
        nq,q=rule(p)
        expect=q8.copy()
        expect[j]=q9[j]
        defect=np.max(abs(q-expect))/scale
        if nq!=n8+1 or defect>3e-12:
            raise ArithmeticError('Fresh one-node geometry did not reproduce the local marginal change')
        e=(q[:,2]-Q0[:,2])/scale
        ans.append(dict(upgraded_parent_index=j,total_nodes=nq,normalized_signed_net_residual=fsum(e),
                        normalized_l1_net_residual=fsum(abs(e)),direct_mask_model_defect=defect,
                        fraction_of_uniform_G8_to_G9_signed_change=d[j,2]/fsum(d[:,2])))
    return dict(x=x,baseline_nodes=[n8,n9],options=ans)


def run(of,refs,lock):
    z=lock.read_bytes()
    a=json.loads(z)
    edges=np.asarray(json.loads(refs[0].read_bytes())['warm_native_edges'])
    b=a['rows'][0]['adjacent_changes'][0]['largest_parent_q_bounds']
    c=a['rows'][0]['adjacent_changes'][1]['largest_parent_q_bounds']
    def ix(t):
        j=int(np.searchsorted(edges,t[0]))
        if edges[j:j+2].tolist()!=t:
            raise ArithmeticError('Trained upgrade does not bind the archived cell')
        return j
    j,k=ix(b),ix(c)
    u=[one(s,r,rec,j,k) for s,r,rec in zip(of,refs,a['rows'],strict=True)]
    return dict(training_x=.01,held_check_x=.05,assigned_and_control_parent_q=[b,c],rows=u,
                held_order_comparison_sha256=hashlib.sha256(z).hexdigest(),
                scope='Direct one-node intervention in a held collision source, not an evolved adaptive chart.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for x in ('states','refs'):
        p.add_argument('--'+x,type=Path,nargs=2,required=True)
    p.add_argument('--lock',type=Path,required=True)
    a=p.parse_args()
    print(json.dumps(run(a.states,a.refs,a.lock),indent=2,allow_nan=False))
