"""Audit parent-momentum support for the archived cold modes without occupation weights."""
import argparse
import gzip
import json

import numpy as np
from scipy.optimize import brentq

from early_pole_archived_shared_cold_window_reference import Cold


def intervals(C,E,q,pole,samples):
    es=np.hypot(q,pole)
    R=[]
    P,D=C.M.r[1:]
    for i in range(P.start,P.stop):
        A=C.M.cells[i]
        for j in range(D.start,D.stop):
            B=C.M.cells[j]
            lo=max(A.lo,B.lo-q); hi=min(A.hi,B.hi+q)
            if hi<=lo:
                continue
            assert min(A.lo,B.lo)>q
            ep=lambda p:E[i,0]+E[i,1]*A.basis(p)
            ed=lambda k:E[j,0]+E[j,1]*B.basis(k)
            def f(p,s):
                return ep(p)-es-ed(np.clip(p+s*q,B.lo,B.hi))
            # Angular lower and upper bounds are separated at the clipping transitions.
            z=np.unique(np.r_[lo,hi,np.linspace(lo,hi,samples+1),
                 [x for x in (B.lo+q,B.hi-q) if lo<x<hi]])
            nodes=z.tolist()
            for side in [-1,1]:
                y=f(z,side)
                for k in np.flatnonzero(np.signbit(y[:-1])!=np.signbit(y[1:])):
                    nodes.append(brentq(lambda v:f(v,side),z[k],z[k+1],xtol=5e-15))
            nodes=np.unique(nodes)
            m=(nodes[:-1]+nodes[1:])/2
            ok=(f(m,-1)>=0.)&(f(m,1)<=0.)
            seg=[(float(x),float(y)) for x,y in zip(nodes[:-1][ok],nodes[1:][ok],strict=True)]
            merged=[]
            for l,u in seg:
                if merged and l-merged[-1][1]<3e-14:
                    merged[-1][1]=u
                else:
                    merged.append([l,u])
            for l,u in merged:
                R.append(dict(parent_cell=[A.lo,A.hi],daughter_cell=[B.lo,B.hi],support=[l,u],width=u-l))
    return R


def run(path,samples):
    with gzip.open(path,'rt') as f:
        raw=json.load(f)['accepted_64']
    out={}
    for name,r in raw.items():
        C=Cold(r['grid'],1,0,len(r['snapshots'][0]['cold_bins']),True,exposure=r['duration'],
               reference_aperture_grid=r['aperture_grid'],archive_projection=name)
        rows=[]
        for t in (r['snapshots'][0],r['snapshots'][-1]):
            E=C.M.weights(np.array(t['hard_moments']),True)[0]
            pole=C.scalar(np.array(t['cold_bins']))[-1]
            modes=[]
            for q in C.q:
                cells=intervals(C,E,q,pole,samples)
                assert cells and all(r['width']>0 for r in cells)
                union=[]
                for l,h in sorted(z['support'] for z in cells):
                    if union and l<=union[-1][1]:
                        union[-1][1]=max(union[-1][1],h)
                    else:
                        union.append([l,h])
                width=sum(h-l for l,h in union)
                sumwidth=sum(z['width'] for z in cells)
                assert width<=sumwidth+3e-15
                modes.append(dict(q=q,parent_support_union_width=width,pairwise_support_width_sum=sumwidth,
                                  support_union=union,intervals=cells))
            rows.append(dict(step=t['step'],pole=pole,modes=modes))
        out[name]=rows
    return dict(samples_per_candidate=samples,scope='Bare angular support of finite dressed cell maps, not a thermal or broadened shell.',
                projections=out)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('bundle')
    p.add_argument('--samples',type=int,default=32)
    a=p.parse_args()
    if a.samples<8:
        p.error('Use at least eight root-bracketing cells.')
    print(json.dumps(run(a.bundle,a.samples),indent=2,allow_nan=False))
