"""Certify a continuous observation-mass range using accepted G6/G7 x=.02 states."""
import argparse
from functools import lru_cache
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_endpoint_gauss_warm_comparison import chart_for,gauss
from early_pole_endpoint_scalar_subshell_chain_receipt import read
from early_pole_endpoint_scalar_warm_evolution_audit import unpack
from early_pole_endpoint_warm_stress_conditional_hull import measure,values


def model(card,base,m0):
    gauss(card,base,m0,.02,expected_steps=2048)
    _,d,repeats=chart_for(card)
    pop=d.ds.w*np.asarray(card['nodal_occupation'][:d.ns])
    ptr=np.r_[0,np.cumsum(repeats)]
    out=[]
    for i,(a,b) in enumerate(zip(d.ds.ed[:-1],d.ds.ed[1:],strict=True)):
        s=slice(ptr[i],ptr[i+1])
        q=d.ds.p[s]
        w=pop[s]
        N=fsum(w)
        h=np.hypot(a,m0)
        E=(b-a)*(b+a)/(np.hypot(b,m0)+h)
        u=(q-a)*(q+a)/(np.hypot(q,m0)+h)/E
        avg=0. if N==0. else fsum(w*u)/N
        if not 0.<=avg<=1.:
            raise ArithmeticError('Invalid reference energy mean')
        out.append((a,b,N,avg,q,w,repeats[i]>1))
    return out


def tighter_error(a,b,mu,m0,n):
    """Bound d² oscillation in reference energy after exact algebraic cancellation."""
    z,m=mu*mu,m0*m0
    u,v=a*a,b*b
    h=(b-a)*(b+a)/(np.hypot(b,m0)+np.hypot(a,m0))
    # Twice differentiate q^(2k)/(3*(q²+mu²)^(k-1/2)) in
    # E0, using q²=E0²-m0². Bound each polynomial term separately.
    n1=(z+m)*v+abs(2*z*z-4*m*z)
    n2=(3*z+m)*v*v+abs(12*z*z-8*m*z)*v+8*m*z*z
    return h*h/(24*n*n)*np.array([n1/(u+z)**2.5,n2/(u+z)**3.5])


def run(native, states, prefix, maximum=.1):
    r=prefix['rows'][-1]
    c,h=read(native)
    if (r['x'],r['common_native_sha256'],c['steps'])!=(.02,h,2048) or not 0.<maximum<=1.:
        raise ArithmeticError('Expected the accepted .02 state and finite observation range')
    d,y=unpack(c,False,.02)
    m0=d.state(.02,y)[0]
    if r['mass']!=m0:
        raise ArithmeticError('Common mass changed')
    charts=[]
    for k,s in enumerate(states):
        z,h=read(s)
        if h!=prefix['audit_receipts'][k][-1]['uncompressed_state_sha256'] or z['warm_gauss_order']!=6+k:
            raise ArithmeticError('Formal observation does not bind the actual orders')
        charts.append(model(z,d,m0))

    cold=[[(a,b,n) for a,b,n,*_,warm in u if not warm] for u in charts]
    if [x[:2] for x in cold[0]] != [x[:2] for x in cold[1]]:
        raise ArithmeticError('Cold constant shapes must be on the same chart')
    nums=np.array([b[2]-a[2] for a,b in zip(*cold,strict=True)])[:,None]
    e=np.array([x[:2] for x in cold[0]])**2

    @lru_cache(None)
    def bands(mu):
        cards=[]
        for u in charts:
            warm=[a for a in u if a[-1]]
            cold=[a for a in u if not a[-1]]
            _,_,_,inner=measure(warm,mu,m0,96)
            er=np.sum([u[2]*tighter_error(u[0],u[1],mu,m0,96) for u in warm],axis=0)
            lo,hi=inner+np.stack((-er,er))
            left,right=np.sum([a[2]*values(np.array(a[:2])**2,mu) for a in cold],axis=0)
            guard=256*np.finfo(float).eps*(hi+right)
            cards.append((lo-guard,hi+guard,np.maximum(0.,lo+left-guard),hi+right+guard))
        return cards

    records=[]
    failed=[]
    def visit(a,b,level):
        u,v=bands(a),bands(b)
        # For either pressure kernel, mu increases the
        # denominator only, so P_G7(b) <= P_G7(mu) and P_G6(mu) <= P_G6(a).
        # Shared cold constant shapes let us combine count differences
        # before bounding their common q dependence.  Do not bound the two
        # potentially large cold backgrounds independently.
        cl=values(e[:,0],b)*nums
        cu=values(e[:,1],a)*nums
        red=np.sum(np.minimum(cl,cu),axis=0)
        blue=np.sum(np.maximum(cl,cu),axis=0)
        bottom=v[1][0]-u[0][1]+red
        upper=u[1][1]-v[0][0]+blue
        bottom=bottom/np.where(bottom>=0.,u[0][3],v[0][2])
        upper=upper/np.where(upper>=0.,v[0][2],u[0][3])
        ready=np.all(bottom>0.)
        if ready or level==8:
            dest=records if ready else failed
            dest.append(dict(mass=[a,b],relative_pressure_and_pseudopressure_lower=bottom.tolist(),
                             relative_pressure_and_pseudopressure_upper=upper.tolist()))
            return
        mid=(a+b)/2
        visit(a,mid,level+1)
        visit(mid,b,level+1)

    visit(0.,maximum,0)
    all_intervals=sorted([r['mass'] for r in records+failed])
    if (all_intervals[0][0] != 0. or all_intervals[-1][1] != maximum
            or any(a[1]!=b[0] for a,b in zip(all_intervals[:-1],all_intervals[1:],strict=True))):
        raise ArithmeticError('Continuous observation partition has a gap or overlap')
    return dict(x=.02,reference_mass=m0,proposed_observation_range=[0.,maximum],
                certified_intervals=records,uncertified_intervals=failed,
                worst_certified_lower=(np.min([p['relative_pressure_and_pseudopressure_lower']
                                               for p in records],axis=0).tolist() if records else None),
                scope='Same frozen actual states; entire mass intervals use pointwise kernel monotonicity, conditional warm count/mean and common cold constant shapes bounded after count subtraction. No physical mass or temporal mapping.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--native',type=Path,required=True)
    p.add_argument('--prefix',type=Path,required=True)
    p.add_argument('--states',type=Path,nargs=2,required=True)
    p.add_argument('--max-mass',type=float,default=.1)
    a=p.parse_args()
    print(json.dumps(run(a.native,a.states,read(a.prefix)[0],a.max_mass),indent=2,allow_nan=False))
