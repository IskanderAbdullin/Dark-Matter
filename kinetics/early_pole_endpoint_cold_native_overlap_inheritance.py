"""Explain whether a coarser overlap certificate survives a real child split."""
import argparse
import hashlib
import json
from math import fsum
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

from early_pole_endpoint_cold_native_placement_refinement import CHANNELS,probability,project
from early_pole_endpoint_cold_non_nested_cdf_contrast import overlap_tv
from early_pole_endpoint_cold_normalized_profile_envelope import geometry


def idx(edges,cuts):
    if not np.all(np.isin(edges,cuts)) or tuple(edges[[0,-1]])!=tuple(cuts[[0,-1]]):
        raise ArithmeticError('Overlap support is not a common refinement')
    return np.searchsorted(edges,cuts[:-1],side='right')-1


def inherit(g,h,G,H,e,f,A,B):
    sa,ta,_=project(g,G,e,A)
    sb,tb,_=project(h,H,f,B)
    p=overlap_tv(g,h,G,H)
    q=np.union1d(e,f)
    cuts=np.union1d(g,h)
    k=len(q)-1
    i,j=idx(e,q),idx(f,q)
    parent=idx(cuts,q)
    cap=np.r_[probability(sa),probability(sb)]
    orig=np.asarray(p['common_measure_by_union_bin'])
    C=np.zeros((len(cap),k))
    ix=np.arange(k)
    C[i,ix]=1.
    C[len(sa)+j,ix]=1.
    E=np.zeros((len(orig),k))
    E[parent,ix]=1.
    # Scale before comparing the exact old allocation with new child
    # mass constraints. Feasible inherited mass may use any descendant
    # of its original overlap interval.
    z=linprog(-np.ones(k),A_ub=np.r_[C,E]*1e4,b_ub=np.r_[cap,orig]*1e4,
              bounds=(0.,None),method='highs',options={'primal_feasibility_tolerance':1e-9,
                                                      'dual_feasibility_tolerance':1e-9})
    if not z.success:
        raise ArithmeticError('Could not lift an inherited overlap allocation')
    x=np.maximum(z.x,0.)
    for row,c in zip(np.r_[C,E].astype(bool),np.r_[cap,orig],strict=True):
        w=fsum(x[row])
        if w>c:
            x[row]*=c/w
    guard=64*np.finfo(float).eps
    K=np.r_[C,E]
    Q=np.r_[cap,orig]
    cert=np.maximum(-z.ineqlin.marginals*1e4,0.)
    cert[:len(sa)]+=max(0.,float(np.max(1.-K.T@cert)))+guard
    if (min(K.T@cert)<1. or max(K@x-Q)>guard or not np.all(np.isfinite(cert))):
        raise ArithmeticError('Independent inherited-cap dual is not feasible')
    bound=fsum(cert*Q)+guard
    loss=[max(0.,fsum(orig)-bound),max(0.,fsum(orig)-fsum(x)+guard)]
    if loss[0]>loss[1] or loss[1]-loss[0]>2e-9:
        raise ArithmeticError('Inherited overlap deficit has no tight certification')
    prices=np.r_[np.asarray(p['native_panel_capacity_prices'][0])[idx(g,e)],
                 np.asarray(p['native_panel_capacity_prices'][1])[idx(h,f)]]
    dual=float(prices@cap)
    old_dual=np.concatenate(p['native_panel_capacity_prices'])@np.r_[probability(G),probability(H)]
    if (min(C.T@prices)<1. or abs(old_dual-dual)>1e-13 or
            np.max(C@x-cap)>1e-13 or np.max(E@x-orig)>1e-13):
        raise ArithmeticError('Inherited primal or lifted dual failed independent validation')
    opt=overlap_tv(e,f,sa,sb)
    c=np.asarray(opt['common_measure_by_union_bin'])
    brought=E@c
    best=fsum(c)
    if np.max(C@c-cap)>1e-13 or dual+1e-13<best:
        raise ArithmeticError('Independent refined overlap did not fit the lifted dual')
    return dict(alternative_projected_minimum_tv_bracket=opt['minimum_TV_bounds'],
                alternative_refined_common_mass=best,alternative_refined_common_mass_by_new_union=c.tolist(),
                alternative_reallocation_from_old_overlap=fsum(np.maximum(orig-brought,0.)),
                parent_compatible_minimum_tv_bracket=p['minimum_TV_bounds'],
                lifted_parent_dual_common_mass_upper=dual+64*np.finfo(float).eps,
                feasible_inherited_common_mass=fsum(x),lost_from_original_primal_bounds=loss,
                inherited_cap_dual_prices=cert.tolist(),
                true_empirical_minimum_TV_perturbation_bound=ta+tb,
                original_union_q=cuts.tolist(),fine_union_q=q.tolist(),
                original_common_mass=orig.tolist(),inherited_common_mass_by_new_union=x.tolist())


def run(root, accepted):
    raw=accepted.read_bytes()
    v=json.loads(raw)
    out=[]
    for mesh in (1,2):
        o=[]; n=[]
        for g in ('fine','third'):
            o.append(geometry(root/f'{g}-d2-endpoint-continuum-cold-{mesh}.jsonl.gz'))
            n.append(geometry(root/f'{g}-d2-placement-eight-endpoint-continuum-cold-{mesh}.jsonl.gz'))
        rows={}
        for name,i,j in CHANNELS:
            (g,G),(h,H)=o
            (e,A),(f,B)=n
            r=inherit(g,h,G[:,i,j],H[:,i,j],e,f,A[:,i,j],B[:,i,j])
            ref=v['compared_measures'][name]
            if (abs(r['parent_compatible_minimum_tv_bracket'][0]-ref['original']['minimum_TV_bracket'][0])>2e-9
                    or abs(r['alternative_projected_minimum_tv_bracket'][0]-ref['parent_conservative']['minimum_TV_bracket'][0])>2e-9
                    or not np.isclose(r['true_empirical_minimum_TV_perturbation_bound'],
                        sum(ref['child_acceptance']['normalized_panel_TV_correction']),atol=3e-13,rtol=0.)):
                raise ArithmeticError('Refinement input differs from the accepted account')
            rows[name]=r
        out.append(rows)
    return dict(accepted_receipt_sha256=hashlib.sha256(raw).hexdigest(),
                scope='Inherited best overlap under exact parent totals. No intra-child shape or physical transport assumption.',
                meshes=out)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True)
    p.add_argument('--accepted',type=Path,required=True)
    a=p.parse_args()
    print(json.dumps(run(a.root,a.accepted),indent=2,allow_nan=False))
