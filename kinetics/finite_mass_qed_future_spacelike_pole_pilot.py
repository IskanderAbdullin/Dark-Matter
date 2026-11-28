"""Verify continuation through zero charged transfer energy at a timelike pole."""
import json

import mpmath as mp
import numpy as np

from finite_mass_qed_timelike_pole_angular_reference import angle,angular_knots
from finite_mass_qed_timelike_pole_batched_absorption_reference import (
    load,pole,radial,spacelike_parent_minimum)


def high_precision():
    mp.mp.dps=65
    rows=[]
    for t in (0.,.001,-.001,.4,-.4):
        for delta in (1e-3,3e-8,3e-13):
            a,b,v,t0=1.4,1.26,abs(t)+delta,t
            s=(t0-v)*(t0+v)
            X=spacelike_parent_minimum(a,b,t0,v,s)
            M,m,T,V=map(mp.mpf,(a,b,t0,v))
            mag=mp.sqrt((V-T)*(V+T))
            p=-(M*M-mag*mag-m*m)/(2*mag)
            # Construct the on-shell parent in a frame in which T has
            # no energy, then transform back at the collinear threshold.
            Y=(V*mp.sqrt(M*M+p*p)+T*p)/mag
            err=float(abs(X/Y-1))
            if err>2e-13:
                raise ArithmeticError(('Spacelike threshold mismatch',t,delta,X,Y,err))
            rows.append(dict(t=t,s=s,minimum=X,boost_frame_defect=err))
    return rows


def discrepancy(a,b):
    return float(np.max(abs(a-b)/np.maximum(np.maximum(abs(a),abs(b)),1e-22)))


def audit(grid):
    d,ch,_=load(grid)
    m=d['checkpoints'][-1]['scalar_pole_over_T']
    k=.0015
    w,R=pole(ch,'T',k)
    q0=np.sqrt((w-m)*(w+m))
    mass=np.sqrt((w-k)*(w+k))
    assert mass>m
    # The past-null pocket is excluded, but the E_s crossing itself
    # must not introduce a discontinuity in the finite signed kernel.
    c=0.
    changes=[]
    for diff in (1e-7,1e-9,1e-11):
        qrange=(q0-diff,q0+diff)
        pp=[]
        for q in qrange:
            Es=np.hypot(m,q)
            a,on,trace=radial(ch,k,w,R,Es,c,'T',8,18,audit=True,qscalar=q,
                              allow_future_spacelike=True)
            if not on or trace>2e-8:
                raise ArithmeticError(('No cross-checked charged continuum',grid,q,trace))
            if diff>0 and q<q0:
                old=radial(ch,k,w,R,Es,c,'T',8,18,qscalar=q)[0]
                if not np.array_equal(old,a):
                    raise ArithmeticError('Opt-in changed the established past-directed result')
            pp.append(a)
        changes.append(dict(q_offsets=diff,coefficient_pair=np.array(pp).tolist(),
                            across_transfer_zero_defect=discrepancy(*pp)))
    if changes[-1]['across_transfer_zero_defect']>1e-5 or changes[0]['across_transfer_zero_defect']<1e-6:
        raise ArithmeticError(('Transfer-crossing behavior',grid,changes))
    tail=[]
    for q in (.004,.015,.03,.08):
        Es=np.hypot(m,q)
        previous=None
        ladder=[]
        for mesh in ((4,(5,12)),(8,(8,18)),(12,(12,24)),(18,(16,32))):
            B,S,N,D=angle(ch,k,w,R,Es,'T',*mesh,q,allow_future_spacelike=True)
            r=discrepancy(previous,B) if previous is not None else None
            ladder.append(dict(mesh=mesh,coeff=B.tolist(),defect=r))
            if r is not None and r<=2e-5:
                break
            previous=B
        if r is None or r>2e-5 or max(B)==0:
            raise ArithmeticError(('Unresolved future spacelike tail',grid,q,N,ladder))
        try:
            angular_knots(ch,k,w,Es,q)
            extra=False
        except ArithmeticError as exc:
            if 'Unexpected sum-of-momenta' not in str(exc):
                raise
            extra=True
        ordinal=np.searchsorted(d['scalar_momentum_edges'],q)-1
        occ=d['endpoint_scalar_occupations'][ordinal]
        denom=B[2]+B[0]/occ
        neutral=float(B[1]/denom) if denom>0 else None
        tail.append(dict(q=q,t=Es-w,raw_A_B_D=B.tolist(),mesh_component_relative=r,
                         conditional_neutral_n=neutral,scalar_histogram_occupation=occ,
                         extra_sum_of_momenta_boundary=extra,angular_panels=N,angular_gate=D,mesh_ladder=ladder))
    if not any(z['extra_sum_of_momenta_boundary'] for z in tail):
        raise ArithmeticError('The new angular branch was not exercised')
    return dict(grid=grid,k=k,w=w,photon_to_scalar_mass_ratio=mass/m,
                zero_charged_transfer_q=q0,crossing=changes,future_spacelike_tail=tail)


if __name__=='__main__':
    print(json.dumps(dict(independent_boost=high_precision(),data=[audit(g) for g in ('fine','third')]),
                     indent=2,allow_nan=False))
