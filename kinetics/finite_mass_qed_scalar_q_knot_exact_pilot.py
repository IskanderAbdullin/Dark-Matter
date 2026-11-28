"""Algebraic enumeration of fixed-photon scalar support knots."""
import json

import numpy as np

from finite_mass_qed_fixed_pole_scalar_q_band_reference import q_knots
from finite_mass_qed_timelike_pole_batched_absorption_reference import load,pole


def roots(a,b,c,m,span):
    # a*sqrt(m*m+q*q) - b*q = c. Squaring gives at most two
    # candidates; return only roots of the unsquared relation.
    if a<=0:
        raise ValueError('Expected positive energy for the declared endpoint branch')
    aa=(a-b)*(a+b)
    bb=-2*b*c
    cc=(a*m-c)*(a*m+c)
    if aa==0:
        rr=() if bb==0 else (-cc/bb,)
    else:
        dd=bb*bb-4*aa*cc
        if dd<0:
            return ()
        if dd==0:
            rr=(-bb/(2*aa),)
        else:
            v=-.5*(bb+np.copysign(np.sqrt(dd),bb))
            rr=(v/aa,cc/v)
    return tuple(q for q in rr if span[0]<q<span[1] and (c+b*q)/a>0 and
                 abs(a*np.hypot(m,q)-b*q-c)<2e-16*max(abs(c),abs(a*m),1e-16))


def exact(d,ch,k,w,s,lo,hi):
    M,m=[np.longdouble(x.m) for x in ch]
    k,w,s,lo,hi=map(np.longdouble,(k,w,s,lo,hi))
    if not 0<lo<hi<k:
        raise ValueError('This enumeration holds within the stated scalar q<k branch')
    out=[]
    counts={'charged_same_forward_direction':0,'charged_other_orientation':0,'past_null':0}
    for c in (-1,1):
        v=roots(np.longdouble(1),np.longdouble(c),w-k,s,(lo,hi))
        counts['past_null']+=len(v)
        out.extend(v)
        for j,x in enumerate(ch):
            own=m if j else M
            for E in map(np.longdouble,x.E):
                if E<own:
                    continue
                p=np.sqrt((E-own)*(E+own))
                a=E-w if j else E+w
                for sign in (-1,1):
                    z=sign*p-k if j else sign*p+k
                    if j:
                        target=(M*M-(a-z)*(a+z)-s*s)/2
                    else:
                        target=((a-z)*(a+z)+s*s-m*m)/2
                    for q in roots(a,z*c,target,s,(lo,hi)):
                        t=np.hypot(s,q)-w
                        father,daughter=(E+t,E) if j else (E,E-t)
                        if father<=M or daughter<=m:
                            continue
                        pp=np.sqrt((father-M)*(father+M))
                        pb=np.sqrt((daughter-m)*(daughter+m))
                        partner=(sign*pb+q*c-k) if j else (sign*pp-q*c+k)
                        p_other=pp if j else pb
                        if abs(abs(partner)-p_other)>2e-13:
                            continue
                        arm='charged_same_forward_direction' if sign==1 and partner>0 else 'charged_other_orientation'
                        counts[arm]+=1
                        out.append(q)
    f=np.asarray(d['scalar_momentum_edges'])
    sampled=f[(f>lo)&(f<hi)].tolist()
    return np.unique(np.array([lo,*out,*sampled,hi],dtype=np.float64)),counts


def run():
    # Both crossings can be inside one bracket whose endpoint residuals
    # have the same sign. The linear case is the past-null boundary.
    wide=tuple(map(np.longdouble,(.3,.8)))
    check=sorted(roots(*map(np.longdouble,(2.,1.,1.75,1.)),wide))
    linear=roots(*map(np.longdouble,(1.,1.,.5,1.)),wide)
    if len(check)!=2 or np.max(abs(np.array(check)-[5/12,3/4]))>1e-15 or linear!=(.75,):
        raise ArithmeticError(('Quadratic or degenerate roots',check,linear))
    rows=[]
    for g in ('fine','third'):
        d,ch,_=load(g)
        k=.0015
        w,R=pole(ch,'T',k)
        m=d['checkpoints'][-1]['scalar_pole_over_T']
        for span in ((1e-6,1.5e-5),(1.5e-5,4e-5),(4e-5,3e-4),
                     (1e-6,3e-4)):
            a,counts=exact(d,ch,k,w,m,*span)
            b=q_knots(d,ch,k,w,m,span,89)
            if a.shape!=b.shape or max(abs(a-b))>1e-15 or counts['charged_other_orientation']:
                raise ArithmeticError(('q crossing mismatch',g,span,len(a),len(b),a.tolist(),b.tolist()))
            rows.append(dict(grid=g,k=k,omega=w,range=span,knots=len(a),
                             max_distance=float(max(abs(a-b))),root_counts=counts))
    return rows


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
