"""Bound the best fixed chart across four frozen warm sources, at a fixed budget."""
import argparse
import hashlib
import json
from math import fsum
from pathlib import Path

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, linprog, milp

from early_pole_endpoint_warm_budgeted_gauss_probe import choose, evaluate, table
from early_pole_endpoint_scalar_subshell_spectrum_audit import read
from early_pole_endpoint_scalar_warm_evolution_audit import unpack


def values(loss, orders, path):
    ix = [orders.index(o) for o in path]
    return np.array([fsum(q[np.arange(len(ix)), ix]) for q in loss])


def fit(loss, orders, budget):
    """Return a feasible integer witness and an independently valid lower bound."""
    U = np.asarray(loss, dtype=float)
    if U.ndim!=3 or U.shape[2]!=len(orders) or not np.all(np.isfinite(U)) or np.min(U)<0:
        raise ArithmeticError('Finite, nonnegative panel losses required')
    nt, n, p = U.shape
    scale = max(np.max(U), 1e-100)
    F = U/scale
    A = np.vstack((np.repeat(np.eye(n),p,axis=1),np.tile(orders,n),F.reshape(nt,-1)))
    rhs = np.r_[np.ones(n),budget,np.zeros(nt)]
    k = n*p
    c = np.r_[np.zeros(k),1.]
    m = milp(c,integrality=np.r_[np.ones(k),0],options={'mip_rel_gap':0.},
             bounds=Bounds(np.zeros(k+1),np.r_[np.ones(k),np.inf]),
             constraints=LinearConstraint(np.c_[A,np.r_[np.zeros(n+1),-np.ones(nt)]],
                                          np.r_[np.ones(n),np.full(1+nt,-np.inf)],rhs))
    if not m.success:
        raise ArithmeticError('No feasible integer minimax chart')
    W = np.rint(m.x[:k]).reshape(n,p)
    if np.max(abs(m.x[:k]-W.ravel()))>1e-7 or np.any(W.sum(axis=1)!=1):
        raise ArithmeticError('Invalid integer optimizer witness')
    path = tuple(orders[i] for i in np.argmax(W, axis=1))
    v = values(U,orders,path)
    if sum(path)>budget:
        raise ArithmeticError('The minimax witness exceeds its budget')
    # For every convex mixture of epochs, the per-cell Bellman optimum is
    # a lower bound on the worst epoch of every feasible fixed chart.
    pool=[values(F,orders,path)]
    lower=0.
    witnesses=[]
    w=np.ones(nt)/nt
    for it in range(128):
        _,(x,q)=choose(np.tensordot(w,F,axes=1),budget,orders)
        z=values(F,orders,q)
        bound=min(float(w@z),x)
        if bound>lower:
            lower=bound
            witnesses=w.tolist()
        if any(np.array_equal(z,a) for a in pool):
            break
        pool.append(z)
        a=linprog(np.r_[np.zeros(nt),-1.],
                  A_ub=np.c_[-np.array(pool),np.ones(len(pool))],b_ub=np.zeros(len(pool)),
                  A_eq=[np.r_[np.ones(nt),0.]],b_eq=[1.],
                  bounds=[(0.,1.)]*nt+[(0.,None)],method='highs')
        if not a.success:
            raise ArithmeticError('The weighted mixture search failed')
        w=np.maximum(a.x[:-1],0.)
        w/=fsum(w)
    guard=64*np.finfo(float).eps*n*max(F.ravel())
    interval=[max(0.,(lower-guard)*scale),float(max(v))]
    if interval[0]>interval[1]+guard*scale:
        raise ArithmeticError('Independent bound exceeds the feasible minimax witness')
    return path,v,interval,witnesses


def run(paths, selection, direct=False):
    states,refs,comparisons=[[read(f) for f in group] for group in paths]
    raw=selection.read_bytes()
    old=json.loads(raw)
    epochs=[.01,.02,.04,.05]
    if ([u['x_end'] for u in states]!=epochs or [u['x'] for u in refs]!=epochs
            or any(r['warm_native_edges']!=old['warm_edges'] for r in refs)):
        raise ArithmeticError('Need the four preregistered stored sources on one native S4 chart')
    orders=tuple(old['orders'])
    ts=[table(r,c) for r,c in zip(refs,comparisons,strict=True)]
    loss=np.asarray([np.array([abs(t[o][:,2]-r[:,2]) for o in orders]).T for t,r in ts])
    old_ns=unpack(states[0],False,.01)[0].ns
    cold=old_ns-loss.shape[1]
    out=[]
    for original in old['choices']:
        cap=original['maximum_total_scalar_nodes']-cold
        oldpath=original['warm_orders_by_cell']
        prev=values(loss,orders,oldpath)
        if not np.allclose(prev[[0,-1]],[x['reference_l1'] for x in original['paired_direct_results']],
                           atol=3e-12,rtol=0.):
            raise ArithmeticError('Frozen original mixed chart fails source receipt')
        new,v,bracket,w=fit(loss,orders,cap)
        oracles=[choose(q,cap,orders)[1][0] for q in loss]
        if min(v-np.array(oracles))<-3e-13 or max(v)>max(prev)+3e-13:
            raise ArithmeticError('Feasible minimax witness contradicts the exact snapshot optima or old chart')
        checks=[]
        if direct:
            checks=[evaluate(a,b,c,new) for a,b,c in zip(states,refs,comparisons,strict=True)]
            if max(abs(v - [q['reference_l1'] for q in checks]))>3e-12:
                raise ArithmeticError('Fresh paired integration disagrees with the minimax objective')
        out.append(dict(maximum_total_nodes=cap+cold,selected_nodes=sum(new)+cold,
                        fixed_warm_orders=new,previous_held_l1=prev.tolist(),minimax_l1=v.tolist(),
                        certified_worst_l1_lower_and_feasible_upper=bracket,lower_bound_epoch_weights=w,
                        exact_independent_snapshot_oracle_l1=oracles,direct_paired_rechecks=checks))
    return dict(epochs=epochs,original_selection_sha256=hashlib.sha256(raw).hexdigest(),budgets=out,
                scope='Retrospective four-snapshot optimization, not a new held test or an evolved rule. The lower bound is from weighted Bellman optima; any remaining gap is not a certificate of exact integer optimality.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for flag in ('states','reference','comparisons'):
        p.add_argument('--'+flag, type=Path, nargs=4, required=True)
    p.add_argument('--selection',type=Path,required=True)
    p.add_argument('--direct',action='store_true')
    a=p.parse_args()
    print(json.dumps(run((a.states,a.reference,a.comparisons),a.selection,a.direct),indent=2,allow_nan=False))
