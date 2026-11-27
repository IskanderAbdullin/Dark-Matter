"""Manufactured varying-mass stress-work and near-rest kinetic redshift ledger."""
import json

import numpy as np
from scipy.integrate import quad

from two_mass_kinetic_offset_deposition_reference import DT, event, shell
from two_mass_near_rest_redshift_ledger_reference import pressure_shells


def test(title,e0,m0,b0,p0):
    e0,m0,b0,p0=map(DT,(e0,m0,b0,p0))
    def masses(x):
        x=DT(x)
        a,b,e=m0*np.exp(DT(.15)*x),b0*np.exp(DT(-.21)*x),e0*np.exp(DT(.23)*x)
        rate=np.array([DT(.15)*a-DT(.21)*b+DT(.23)*e,DT(.15)*a,-DT(.21)*b],dtype=DT)
        return np.array([a+b+e,a,b],dtype=DT),rate,e
    grid=np.geomspace(p0/DT(3e4),p0*25,92)
    P=np.array([.75*p0,1.22*p0],dtype=DT)
    Q=(P,grid,grid)
    occ=[np.ones(1,dtype=DT),np.zeros(len(grid)-1,dtype=DT),np.zeros(len(grid)-1,dtype=DT)]
    def state(x):
        M,D,e=masses(x)
        a=np.exp(DT(x))
        tmp=[pressure_shells(q/a,m) for q,m in zip(Q,M,strict=True)]
        K,p,tr=(np.array([v[i]@n for v,n in zip(tmp,occ,strict=True)],dtype=DT) for i in range(3))
        count=np.array([np.sum(n) for n in occ],dtype=DT)
        if np.max(abs(K+tr-3*p))>3e-16*K.sum():
            raise ValueError("Both species' pressure/trace relations did not close.")
        return dict(excess=e,kinetic=K,pressure=p,trace_deficit=tr,
                    offset_work=(D/M)@tr,rest_work=D@count,rest=M@count)
    ini=state(0.)
    norm=ini['kinetic'].sum()+e0
    norms=(norm,norm,ini['rest'])
    W=np.zeros(3,dtype=DT)
    events=DT(0.)
    x0=0.
    C=DT(0.)
    out=[]
    for a,f,mu in ((1.3,.15,.3),(2.1,.22,-.4),(4.2,.18,.1),(8.,0.,.6)):
        x=np.log(a)
        at=state(x)
        def columns(t):
            v=state(t)
            return [3*v['pressure'].sum(),v['offset_work'],v['rest_work']]
        W+=np.array([DT(quad(lambda t:float(columns(t)[i]/n),x0,x,epsabs=1e-12,epsrel=1e-12)[0])*n
                     for i,n in enumerate(norms)],dtype=DT)
        if f:
            f=DT(f)
            M,_,ee=masses(x)
            a=np.exp(DT(x))
            kp=shell(P/a,M[0])[0][0]
            mom=np.sqrt(kp*(kp+2*M[0]))
            _,values,_=event(ee,M[1],M[2],mom,mu)
            for n,m,v in zip(occ[1:],M[1:],values,strict=True):
                K=shell(grid/a,m)[0]
                j=int(np.searchsorted(K,v,side='right')-1)
                if not 0<=j<len(K)-1:
                    raise ValueError("Escaped daughter offset gamut.")
                z=(v-K[j])/(K[j+1]-K[j])
                n[j:j+2]+=f*np.array([1-z,z])
            occ[0][0]-=f
            C+=f
            events+=ee*f
            after=state(x)
            if abs(after['kinetic'].sum()-at['kinetic'].sum()-ee*f)>4e-15*norm:
                raise ValueError("Independently redressed event missed the finite-excess transfer.")
            at=after
        if min(x.min() for x in occ)<0 or max(abs(x.sum()-C) for x in occ[1:])>1e-16:
            raise ValueError("Population or daughter number check failed.")
        ek=(at['kinetic'].sum()-ini['kinetic'].sum()+W[0]+W[1]-events)/norm
        er=(at['rest']-ini['rest']-W[2]+events)/norms[2]
        if max(abs(ek),abs(er))>2e-13:
            raise ValueError("Split varying mass/rest/trace ledger did not close.")
        out.append(dict(a=float(np.exp(x)),
                        comoving_kinetic_excess=list(map(float,at['kinetic'])),
                        comoving_trace_deficit=list(map(float,at['trace_deficit'])),
                        decayed_parent_count=float(C),
                        accumulated_mass_offset_work=float(W[1]),
                        accumulated_rest_work_from_external_mass_protocol=float(W[2]),
                        pressure_work=float(W[0]),cumulative_decay_rest_excess=float(events),
                        normalized_kinetic_defect=float(abs(ek)),normalized_rest_defect=float(abs(er))))
        x0=x
    return dict(case=title,initial=[float(e0),float(m0),float(b0),float(p0)],rows=out)


if __name__=='__main__':
    print(json.dumps(dict(scope="Prescribed masses and events with externally supplied work, not physical evolution or rates.",
        cases=[test('ordinary',.093,.061,.122,.33),test('near-rest',2e-23,.061,.122,4e-12)]),
                     indent=2,allow_nan=False))
