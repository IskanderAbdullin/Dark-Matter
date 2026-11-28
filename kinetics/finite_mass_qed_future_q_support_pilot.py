"""Enumerate scalar q support of the archived timelike pole beyond the soft bands."""
import json

import numpy as np

from finite_mass_qed_scalar_q_knot_exact_pilot import exact as soft_exact, roots
from finite_mass_qed_timelike_pole_batched_absorption_reference import load,pole


def crossings(ch,k,w,ms,lo,hi):
    """Exact charged parent/daughter energy-edge order changes in the angular interior."""
    k,w,ms=np.longdouble(k),np.longdouble(w),np.longdouble(ms)
    E,B=[np.asarray(x.E,dtype=np.longdouble) for x in ch]
    M,m=[np.longdouble(x.m) for x in ch]
    E=E[E>M]
    B=B[B>m]
    t=E[:,None]-B[None,:]
    es=w+t
    q=np.sqrt(np.maximum(0.,(es-ms)*(es+ms)))
    enabled=(es>ms)&(lo<q)&(q<hi)
    p=np.sqrt((E-M)*(E+M))[:,None]
    b=np.sqrt((B-m)*(B+m))[None,:]
    out=[]
    for name,V in (('same',abs(p-b)),('opposite',p+b)):
        keep=enabled&(abs(q-k)<V)&(V<q+k)&(abs(t)<V)
        for i,j in zip(*np.where(keep),strict=True):
            c=((k-V[i,j])*(k+V[i,j])+q[i,j]**2)/(2*k*q[i,j])
            if not -1<c<1:
                raise ArithmeticError('Crossing located outside the physical angular interior')
            out.append(dict(q=float(q[i,j]),c=float(c),parent=float(E[i]),daughter=float(B[j]),
                            charged_collinear_arrangement=name))
    return out


def candidates(d,ch,k,w,s,lo,hi):
    """Return optional collinear and null nodes; no sign scan or rate evaluation."""
    M,m=map(np.longdouble,[x.m for x in ch])
    k,w,s,lo,hi=map(np.longdouble,(k,w,s,lo,hi))
    if not 0<=lo<hi or not M>m or (w-k)*(w+k)<=s*s:
        raise ValueError('Outside the past-timelike exclusion hypotheses')
    shell=np.asarray(d['scalar_momentum_edges'])
    out={'scalar_edges':shell[(lo<shell)&(shell<hi)].tolist(),
         'charged_collinear':[],'charged_opposite':[],'null':[],'zero_transfer':[],
         'interior_radial_edge_order':[r['q'] for r in crossings(ch,k,w,s,lo,hi)]}
    z=np.sqrt((w-s)*(w+s))
    if lo<z<hi:
        out['zero_transfer'].append(z)
    for c in (-1,1):
        for sig in (-1,1):
            for q in roots(np.longdouble(1),np.longdouble(sig*c),w-sig*k,s,(lo,hi)):
                out['null'].append(q)
        for isb,charm in enumerate(ch):
            own=m if isb else M
            for E0 in charm.E:
                E=np.longdouble(E0)
                if E<own:
                    continue
                p=np.sqrt((E-own)*(E+own))
                a=E-w if isb else E+w
                for sign in (-1,1):
                    z=sign*p-k if isb else sign*p+k
                    C=(M*M-(a-z)*(a+z)-s*s)/2 if isb else ((a-z)*(a+z)+s*s-m*m)/2
                    for q in roots(a,z*c,C,s,(lo,hi)):
                        t=np.hypot(s,q)-w
                        father,daughter=(E+t,E) if isb else (E,E-t)
                        if father<=M or daughter<=m:
                            continue
                        pp=np.sqrt((father-M)*(father+M))
                        pb=np.sqrt((daughter-m)*(daughter+m))
                        partner=(sign*pb+q*c-k) if isb else (sign*pp-q*c+k)
                        p_other=pp if isb else pb
                        # Squaring can create a root on an unphysical shell
                        # branch. Require the signed one-dimensional partner.
                        if abs(abs(partner)-p_other)>3e-13*max(1.,p_other):
                            continue
                        out['charged_collinear' if partner*sign>0 else 'charged_opposite'].append(q)
    nodes=np.unique(np.asarray([lo,hi,*[x for y in out.values() for x in y]],dtype=float))
    return nodes,out



def verify_order(ch,k,w,ms,entry):
    # Reconstruct the two angular edges from *separate* fixed charged
    # energies at nearby q. They must exchange order at the reported q.
    k,w,ms=map(np.longdouble,(k,w,ms))
    M,m=map(np.longdouble,[z.m for z in ch])
    P,B,q0,c0=map(np.longdouble,[entry[x] for x in ('parent','daughter','q','c')])
    def edges(q):
        t=np.hypot(ms,q)-w
        e=np.array([P,B+t])
        pe=np.sqrt((e-M)*(e+M))
        be=np.sqrt((e-t-m)*(e-t+m))
        v=abs(pe-be) if entry['charged_collinear_arrangement']=='same' else pe+be
        return ((k-v)*(k+v)+q*q)/(2*k*q)
    eps=1e-7*q0
    a,b,c=edges(q0-eps),edges(q0),edges(q0+eps)
    if (a[0]-a[1])*(c[0]-c[1])>=0 or np.max(abs(b-c0))>2e-12:
        raise ArithmeticError(('The two angular edges did not exchange',entry,a,b,c))
    entry['left_and_right_angle_differences']=[float(a[0]-a[1]),float(c[0]-c[1])]
    return entry

def run():
    data=[]
    for g in ('fine','third'):
        d,ch,_=load(g)
        k=.0015
        w,R=pole(ch,'T',k)
        ms=d['checkpoints'][-1]['scalar_pole_over_T']
        d0=ch[1].m
        pmax=np.hypot(ch[0].m,max(d['endpoint_parent_momentum_edges']))
        # If f_P is nonzero in the archive, total incoming energy is bounded;
        # the final charged daughter always carries at least its rest mass.
        qmax=np.sqrt((pmax+w-d0-ms)*(pmax+w-d0+ms))
        Qcap=max(qmax,d['scalar_momentum_edges'][-1])
        spans=((1e-6,.0003),(.0003,.0015),(.0015,.004),(.004,.015),
               (.015,.1),(.1,1.),(1.,Qcap))
        allnodes,allkinds=candidates(d,ch,k,w,ms,0.,Qcap)
        null=sorted(float(v) for v in allkinds['null'])
        if len(null)!=2 or not null[0]<null[1]<float(allkinds['zero_transfer'][0]):
            raise ArithmeticError(('Unexpected past-null pocket geometry',g,null))
        for q in null:
            t=np.hypot(ms,q)-w
            if not t<0 or min(abs(t+abs(q*c-k)) for c in (-1,1))>3e-16:
                raise ArithmeticError(('Uncertified extreme null',g,q,t))

        def excluded(q):
            t=np.hypot(ms,q)-w
            c0=((k+t)*(k-t)+q*q)/(2*k*q)
            return (1.-np.clip(c0,-1.,1.))/2.

        area=[float(excluded(q)) for q in (null[0]*.99,sum(null)/2,null[1]*1.01)]
        if not (area[0]==area[2]==0. and 0.<area[1]<1.):
            raise ArithmeticError(('Wrong past timelike aperture',g,area))
        pieces=[]
        for span in spans:
            a,count=candidates(d,ch,k,w,ms,*span)
            if span==spans[0]:
                b,s=soft_exact(d,ch,k,w,ms,*span)
                if not np.array_equal(a,b):
                    raise ArithmeticError(('Soft limit lost',g,len(a),len(b),np.max(abs(a-b))))
            c={i:len(v) for i,v in count.items()}
            if np.min(np.diff(a))<=0:
                raise ArithmeticError('Past-null and endpoint brackets fell out of order')
            pieces.append(dict(q_span=span,panels=len(a)-1,root_counts=c))
        data.append(dict(grid=g,k=k,omega=w,archive_scalar_support=d['scalar_momentum_edges'][-1],
                         forward_max_q_energy_bound=qmax,all_source_max_q_bound=Qcap,
                         kinematic_null_q_endpoints=null,
                         excluded_cosine_fraction_at_three_probes=area,
                         first_opposite_charged_orientation=float(min(allkinds['charged_opposite'])),
                         interior_radial_edge_order=[verify_order(ch,k,w,ms,v)
                                                     for v in crossings(ch,k,w,ms,0.,Qcap)],
                         exact_retained_panel_counts=pieces))
    return data


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
