"""Bound cold profile CDF separation without assuming within-panel positions."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_endpoint_cold_normalized_profile_envelope import build,fit_state,geometry,shares


def cut(card,keep,mode):
    if not np.any(keep):
        return [0.,0.]
    if np.all(keep):
        return [1.,1.]
    base=np.asarray(card['base'])
    if mode=='fixed':
        f=fsum(base[keep])
        return [f,f]
    c=np.asarray(card['offset'])
    place=lambda z:np.column_stack((np.sum(z[:,keep],axis=1),np.sum(z[:,~keep],axis=1)))
    P=[place(z) for z in card[mode]]
    R=np.array([fsum(c[keep]),fsum(c[~keep])])
    exact=np.array([fsum(base[keep]),fsum(base[~keep])])*card['total']
    return shares(R,P,exact)[0].tolist()


def distances(e,f,A,B,mode):
    edge=np.union1d(e,f)
    out=[]
    for l,h in zip(edge[:-1],edge[1:],strict=True):
        q=(l+h)/2
        r=[]
        for x,d in ((e,A),(f,B)):
            # Every q inside this open interval belongs to the same topology.
            # Entire-left panels must contribute; an intersected panel's response
            # may occur at either side, so do not interpolate.
            lo=cut(d,x[1:]<=q,mode)[0]
            hi=cut(d,x[:-1]<q,mode)[1]
            if not 0<=lo<=hi<=1:
                raise ArithmeticError('Invalid CDF support interval')
            r.append([lo,hi])
        u,v=r
        a=max(0.,u[0]-v[1],v[0]-u[1])
        b=max(abs(u[1]-v[0]),abs(v[1]-u[0]))
        out.append(dict(q=[l,h],left_right_CDF_intervals=r,guaranteed_absolute_separation=a,
                        largest_allowed_absolute_separation=b))
    u=fsum((k['q'][1]-k['q'][0])*k['guaranteed_absolute_separation'] for k in out)
    v=fsum((k['q'][1]-k['q'][0])*k['largest_allowed_absolute_separation'] for k in out)
    return dict(Wasserstein1_q_lower_upper=[u,v],
                maximum_guaranteed_CDF_separation=max(z['guaranteed_absolute_separation'] for z in out),
                nonzero_guaranteed_intervals=[z for z in out if z['guaranteed_absolute_separation']>0])


def run(a,b,states,archive,l):
    if tuple(l) not in ((2,4),(4,8),(8,16)):
        raise ArithmeticError('Need two accepted adjacent scalar charts')
    P=[fit_state(p,archive,n,h) for p,n,h in zip(states,l,(a,b),strict=True)]
    out=[]
    for j in (0,1):
        u,v=map(geometry,(a[j],b[j]))
        e,f=u[0],v[0]
        if [e[0],e[-1]]!=[f[0],f[-1]]:
            raise ArithmeticError('CDF bands differ')
        x,y=[build(g,r,g[0],raw=True) for g,r in zip((u,v),P,strict=True)]
        out.append({m:distances(e,f,x,y,m) for m in ('fixed','count_only','count_and_mean_excess')})
    a,b=out
    h=lambda x:np.asarray(x['Wasserstein1_q_lower_upper'])
    if max(np.max(abs(h(a[k])-h(b[k]))) for k in b)>2e-12:
        raise ArithmeticError('CDF support bound did not close on paired meshes')
    u,v,w=map(lambda k:h(b[k]),('count_only','count_and_mean_excess','fixed'))
    if not (u[0]-1e-14<=v[0]<=w[0]+1e-14 and u[1]+1e-14>=v[1]>=w[1]-1e-14):
        raise ArithmeticError('Nested CDF families have inconsistent bounds')
    return dict(x=.05,grid='fine',charge_div=2,scalar_levels=l,**b,
                scope='Outward CDF envelopes over finite panel amounts and unconstrained within-panel position. Lower bound zero does not exhibit a common distribution.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('left','right','states'):
        p.add_argument('--'+k,type=Path,nargs=2,required=True)
    p.add_argument('--levels',type=int,nargs=2,required=True)
    p.add_argument('--archive',type=Path,required=True)
    a=p.parse_args()
    print(json.dumps(run(a.left,a.right,a.states,a.archive,a.levels),indent=2,allow_nan=False))
