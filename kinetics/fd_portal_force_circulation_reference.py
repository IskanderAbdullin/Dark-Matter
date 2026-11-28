"""Check whether the displaced two-bath flow is the gradient of a force potential."""
import json
from math import fsum

import numpy as np

from fd_portal_differential_reciprocity_reference import integrate, PositivePortal
from positive_fd_portal_quadrature_reference import rule


def rectangle(portal,T,rect,n,fraction):
    (a,b),(c,d)=rect
    x,w=rule(n)
    u=(a+b)/2+(b-a)*x/2
    v=(c+d)/2+(d-c)*x/2
    ua,va=(b-a)*w/2,(d-c)*w/2
    def flux(p,q):
        if fraction==1.:
            # Use the pre-existing independent rate evaluator for the full
            # FD boundary, rather than the differentiated kernel.
            return portal.moments(T,T/(1-q),p/2)[:2]/[1.,T]
        return integrate(portal,T,p,q,dark_blocking_fraction=fraction)[0][:2]
    edges=[fsum(s*flux(p,c)[0] for p,s in zip(u,ua,strict=True)),
           fsum(s*flux(b,q)[1] for q,s in zip(v,va,strict=True)),
           -fsum(s*flux(p,d)[0] for p,s in zip(u,ua,strict=True)),
           -fsum(s*flux(a,q)[1] for q,s in zip(v,va,strict=True))]
    skew=[]
    for p,s in zip(u,ua,strict=True):
        for q,t in zip(v,va,strict=True):
            r=integrate(portal,T,p,q,dark_blocking_fraction=fraction)
            skew.append(-s*t*r[2])
    dS=fsum(skew)
    flow=fsum(edges)
    scale=fsum(abs(x) for x in edges)
    if (scale<=0. or abs(dS-flow)>3e-8*scale or
        fraction==0. and (dS!=0. or abs(flow)>3e-8*scale)):
        raise ArithmeticError(('Boundary and independent curl did not close',flow,dS,scale))
    return dict(boundary_CCW_circulation=flow,integrated_CCW_curl=dS,
                oriented_boundary_segments=edges,absolute_boundary_scale=scale,
                stokes_relative_closure=abs(flow-dS)/scale)


def run():
    q=(PositivePortal((13,36,26)),PositivePortal((17,44,32)))
    data=[]
    for label,T,rect,sign in (
            ('same_sign_hot',.05,((.24,.72),(.06,.18)),-1),
            ('same_sign_cold',.05,((-.72,-.24),(-.18,-.06)),1),
            ('mixed_force',.30,((-1.64,-1.28),(.11,.17)),None)):
        trials=[]
        for p,n in zip(q,(4,6),strict=True):
            ds={str(v):rectangle(p,T,rect,n,v) for v in (0.,.5,1.)}
            reference=ds['1.0']
            full=reference['integrated_CCW_curl']
            if full==0. or sign is not None and np.sign(full)!=sign:
                raise ArithmeticError('Observed circulation contradicted the pointwise sign identity')
            mix=ds['0.5']['integrated_CCW_curl']
            if abs(mix/full-.5)>3e-14:
                raise ArithmeticError('Only the differentiated dark blocking can carry this curl')
            trials.append(ds)
        old,new=trials
        u,v=[k['1.0']['integrated_CCW_curl'] for k in trials]
        if abs(u/v-1.)>3e-5:
            raise ArithmeticError('Independent quadrature does not resolve the same FD circulation')
        data.append(dict(label=label,fixed_T_SM=T,coordinate_box=rect,
                         paired_fd_curl_drift=abs(u/v-1.),quadratures=trials))
    return dict(coordinates=['2 eta_D','1 - T_SM / T_D'],oriented_flow=['F','Q / T_SM'],
                rectangles=data,scope='Local charged-lepton timelike FD kernel. Blocking toggle keeps the same visible occupations, phase-space kernel and affinity; no alteration of numerical bath trajectory or omitted physical processes.')


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
