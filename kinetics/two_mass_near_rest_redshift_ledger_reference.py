"""Manufactured comoving decay sequence with an explicitly split kinetic FRW ledger."""
import json

import numpy as np
from scipy.integrate import quad
from scipy.special import roots_legendre

from two_mass_kinetic_offset_deposition_reference import DT, deposit, event, shell


def pressure_shells(edges,m,order=24):
    x,w=(np.asarray(v,dtype=DT) for v in roots_legendre(order))
    lo,hi=edges[:-1,None],edges[1:,None]
    p=lo+(hi-lo)*(1+x)/2
    p2=p*p
    W=p2*w
    W/=W.sum(axis=1,keepdims=True)
    E=np.sqrt(m*m+p2)
    K=p2/(E+m)
    return np.sum(W*K,axis=1),np.sum(W*p2/(3*E),axis=1),np.sum(W*m*K/E,axis=1)


def test(name,e,m,b,comoving_parent):
    e,m,b=map(DT,(e,m,b))
    M=m+b+e
    # Choose daughter grids once in comoving momentum. Do not resize or
    # advect their coordinates after injecting particles.
    anchor=max(np.sqrt(event(e,m,b,0.,0.)[0][0]),DT(comoving_parent))
    q=np.geomspace(anchor/DT(3e4),anchor*20,90)
    p=np.asarray([.72*comoving_parent,1.21*comoving_parent],dtype=DT)
    edges=(p,q,q)
    masses=(M,m,b)
    counts=[np.array([1.],dtype=DT),np.zeros(len(q)-1,dtype=DT),np.zeros(len(q)-1,dtype=DT)]

    def state(x):
        a=np.exp(DT(x))
        z=[pressure_shells(v/a,k) for v,k in zip(edges,masses,strict=True)]
        K,P,T=[np.array([arr[j]@n for arr,n in zip(z,counts,strict=True)],dtype=DT) for j in range(3)]
        if np.max(abs(3*P-K-T))>1e-15*sum(K):
            raise ValueError("Specieswise kinetic pressure and two-mass trace-deficit closure failed.")
        return K,P,T

    begin=state(0.)[0].sum()
    normal=begin+e
    records=[]
    W=DT(0.)
    D=DT(0.)
    x0=0.
    for a,fraction,mu in ((1.2,.12,.2),(1.9,.27,-.4),(3.6,.2,.1),(6.,.23,.62),(8.,0.,.9)):
        x=np.log(a)
        qa=np.exp(DT(x))
        pre=state(x)
        recovered=DT(quad(lambda t:float(3*state(t)[1].sum()/normal),
                          x0,x,epsabs=3e-13,epsrel=3e-13)[0])*normal
        if abs(state(x0)[0].sum()-pre[0].sum()-recovered)>8e-13*normal:
            raise ValueError("Exact redshift trace ledger did not reproduce the independent pressure integral.")
        W+=recovered
        event_err=0.
        if fraction:
            d=DT(fraction)
            kp=shell(p/qa,M)[0][0]
            lab=np.sqrt(kp*(kp+2*M))
            _,K,_=event(e,m,b,lab,mu)
            projected=deposit(e,m,b,lab,mu,(q/qa,q/qa))
            for c,off,mass in zip(counts[1:],K,(m,b),strict=True):
                vals=shell(q/qa,mass)[0]
                i=int(np.searchsorted(vals,off,side='right')-1)
                w=(off-vals[i])/(vals[i+1]-vals[i])
                c[i:i+2]+=d*np.array([1-w,w])
            counts[0][0]-=d
            D+=d
            post=state(x)
            event_err=float(abs(post[0].sum()-pre[0].sum()-e*d)/normal)
            if event_err>1e-14 or projected['energy_relative_error']>2e-15:
                raise ValueError("Collision kinetic gain did not match lost parent rest excess.")
        r=state(x)
        if min(c.min() for c in counts)<0. or max(abs(D-c.sum()) for c in counts[1:])>1e-16:
            raise ValueError("Paired daughter number or positive distribution failed.")
        residual=float(abs(r[0].sum()-begin+W-e*D)/normal)
        if residual>1e-12:
            raise ValueError("Combined pressure, kinetic and rest-change ledger failed.")
        rest_rounded=float(sum(k*c.sum() for k,c in zip(masses,counts,strict=True)))
        naive=float(rest_rounded+float(sum(r[0]))-rest_rounded)
        if name=='near_rest' and (naive!=0. or r[0].sum()<=0.):
            raise ValueError("Near-rest fixture no longer exposes the float64 cancellation.")
        records.append(dict(a=a,decayed_parent_count=float(D),cumulative_pressure_work=float(W),
                            per_species_comoving_kinetic_excess=list(map(float,r[0])),
                            per_species_comoving_pressure=list(map(float,r[1])),
                            per_species_comoving_trace_deficit=list(map(float,r[2])),
                            comoving_rest_change=float(-e*D),naive_recovered_kinetic_excess=naive,
                            fractional_split_ledger_residual=residual,fractional_event_ledger_residual=event_err))
        x0=x
    return dict(name=name,inputs=[float(e),float(m),float(b),comoving_parent],
                kinetic_error_norm=float(normal),stored_parent_rest_mass_excess_is_separate=True,
                initial_parent_absolute_shell_edges_collapse=bool(np.hypot(float(M),float(p[0]))==np.hypot(float(M),float(p[1]))),
                ledger=records)


if __name__=='__main__':
    print(json.dumps(dict(scope="Prescribed injection and redshift; no solved rates, Hubble history or abundance.",
        cases=[test('ordinary',.093,.061,.122,.33),test('near_rest',2e-23,.061,.122,4e-12)]),
        indent=2,allow_nan=False))
