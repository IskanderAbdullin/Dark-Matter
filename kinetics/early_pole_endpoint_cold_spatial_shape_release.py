"""Bound the common-shape conclusion under independent per-chart releases."""
import argparse
import hashlib
import json
from math import fsum
from pathlib import Path

import numpy as np
from early_pole_endpoint_cold_spatial_common_shape_geometry import band


A=np.array([-.5,.5,-.5,.5])
B=np.array([-.5,-.5,.5,.5])


def adversary(e,p):
    w=np.diff(e)
    u=A@p
    v=B@p
    q=np.r_[0.,np.cumsum(u)]
    r=np.r_[0.,np.cumsum(v)]
    # Shared atom mixture attains the minimum population CDF L1 in the
    # limit of edge placement, before releasing mass by kernel row.
    z=np.zeros(len(u))
    for k,(l,h) in enumerate(zip(q[:-1],q[1:],strict=True)):
        look=[(abs(l),0.),(abs(h),1.)]
        if l*h<0.:
            f=-l/(h-l)
            look.append((0.,f))
        z[k]=min(look)[1]
    guard=2e-6
    def ratio(eps,reverse=False):
        F=np.array([1.,1.,0.,0.])[:,None]
        if reverse:
            F=1.-F
        mid=(1-eps)*z+eps*F
        # Use atoms at inward fractions guard and 1-guard. The remaining
        # tiny left/right lengths carry the exact cumulative edge values.
        v_mid=r[:-1]+B@(p*mid)
        u_mid=q[:-1]+A@(p*mid)
        def calc(cuts,interior):
            return fsum(w*(guard*(abs(cuts[:-1])+abs(cuts[1:]))+(1-2*guard)*abs(interior)))
        return calc(q,u_mid),calc(r,v_mid)
    return ratio


def run(rects,shapes):
    rows=[]
    ref=json.loads(shapes.read_text())
    e=np.asarray(ref['common_q_edges'])
    for p0,g in zip(rects,ref['pairs'],strict=True):
        raw=p0.read_bytes()
        x=json.loads(raw)
        if (x['levels']!=g['levels'] or list(e)!=x['common_q_edges']
                or hashlib.sha256(raw).hexdigest()!=g['source_sha256']):
            raise ArithmeticError('Release experiment does not bind the accepted spatial display')
        out={}
        for k,r in x['observables'].items():
            p=np.asarray(r['normalized_profile_rectangle']).reshape(4,-1)
            l=band(e,A@p)['flow_cdf_L1_q_lower_upper'][0]
            u=band(e,B@p)['flow_cdf_L1_q_lower_upper'][1]
            wire=g['observables'][k]['paired_mesh_cdf_L1_change_cap_q']
            # For the released CDFs, δP±δK are differences of two
            # corner perturbations (diagonal or off-diagonal). At a
            # given q, |δP|+|δK|=max(|δP+δK|,|δP-δK|).
            # For arbitrary Bernoulli CDF values x,y,z in [0,1],
            # |a(x-z)-b(y-z)|<=max(a,b), by the extrema of the cube.
            # Crucially both released corners use the SAME remaining
            # shape, so the sum of separate W1 shift bounds is too loose.
            defect=fsum(np.max(p,axis=0)*np.diff(e))
            lower=(l-u-2*wire)/defect
            make=adversary(e,p)
            c=[make(1.,i)[1]-make(1.,i)[0] for i in (False,True)]
            if min(lower,defect)<=0. or max(c)<=2*wire:
                raise ArithmeticError('No useful conservative bound or contrast reversal')
            reverse=(c[1]>c[0])
            a,b=0.,1.
            for _ in range(45):
                m=(a+b)/2
                v,w=make(m,reverse)
                if w-v>2*wire:
                    b=m
                else:
                    a=m
            b*=1.001
            v,w=make(b,reverse)
            if not lower<b<1. or w-v<=2*wire:
                raise ArithmeticError('No independent sufficient/reversal bracket')
            out[k]=dict(parent_l1_floor_q=l,kernel_l1_cap_q=u,release_shift_budget_q=defect,
                        uniform_release_sufficient_to_preserve_dominance=lower,
                        constructive_release_exceeding_dominance=b,
                        witness_population_kernel_amplitudes_q=[v,w],release_advects_reversed=reverse)
        rows.append(dict(levels=x['levels'],observables=out))
    return dict(pairs=rows,scope='The same fraction in each bin can take an independent shape for each numerical corner; common remaining shape, unchanged bin counts. Bounds and explicit counterexample, not a physical distribution or an exact threshold.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--rectangles',nargs=3,type=Path,required=True)
    p.add_argument('--shapes',type=Path,required=True)
    a=p.parse_args()
    print(json.dumps(run(a.rectangles,a.shapes),indent=2,allow_nan=False))
