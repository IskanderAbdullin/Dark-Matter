"""Conditional pressure envelopes for exact G8, G9 and mixed300 checkpoints."""
import argparse
import hashlib
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_endpoint_gauss_warm_comparison import chart_for,gauss
from early_pole_endpoint_scalar_subshell_chain_receipt import read
from early_pole_endpoint_scalar_warm_evolution_audit import unpack
from early_pole_endpoint_warm_stress_conditional_hull import measure


def conditioned(base,card,m0,x,control):
    obs=gauss(card,base,m0,x,expected_steps=2048,budget_control=control)[0]
    tag,d,repeats=chart_for(card,budget_control=control)
    e=d.ds.ed
    y=np.asarray(card['nodal_occupation'])
    w=d.ds.w*y[:d.ns]
    K=np.r_[0,np.cumsum(repeats)]
    output=[]
    for i,(a,b) in enumerate(zip(e[:-1],e[1:],strict=True)):
        j=slice(K[i],K[i+1]); q=d.ds.p[j]
        N=fsum(w[j])
        ha=np.hypot(a,m0)
        E=(b-a)*(b+a)/(np.hypot(b,m0)+ha)
        u=(q-a)*(q+a)/(np.hypot(q,m0)+ha)/E
        mu=0. if not N else fsum(w[j]*u)/N
        if not 0.<=mu<=1.:
            raise ArithmeticError('Nodal conditional energy is outside its own cell')
        output.append((a,b,N,mu,q,w[j],repeats[i]>1))
    return tag,output,obs


def state(a,b):
    return ['strict_positive' if l>0 else 'strict_negative' if h<0 else 'both_constructively_feasible'
            if f<0<g else 'inconclusive' for (l,h),(f,g) in zip(a,b,strict=True)]


def run(native,control,accepted,states):
    n,nh=read(native)
    raw=accepted.read_bytes()
    lock=json.loads(raw)
    r=lock['rows'][-1]
    x=r['x']
    D=[read(y) for y in states]
    if (x not in (.04,.05) or lock['charts']!=['G8','G9','mixed300']
            or r['steps']!=2048 or [h for _,h in D]!=r['states_sha256']):
        raise ArithmeticError('Need exact verified three-chart terminal state')
    d,y=unpack(n,False,x)
    m0=d.state(x,y)[0]
    if m0!=r['common_native_mass']:
        raise ArithmeticError('Must measure at the same native archive mass')
    charts=[conditioned(d,blob,m0,x,c) for (blob,_),c in zip(D,[None,None,control],strict=True)]
    if [k for k,_,_ in charts]!=[8,9,'mixed300']:
        raise ArithmeticError('Chart selections do not match the archive')
    v=np.asarray([m for _,_,m in charts])
    if max(abs((v[1:]/v[0]-1.-r['ratios_to_G8']).ravel()))>2e-15:
        raise ArithmeticError('Published reference count and pressure do not replay')
    out=[]
    for m in (0.,m0,.01,.1,1.,10.):
        a,b=[[measure(c,m,m0,nq) for _,c,_ in charts] for nq in (64,256)]
        rows=[]
        for (lo,hi),w,e,inner in b:
            rows.append(dict(outer=[lo.tolist(),hi.tolist()],constructive_inner=inner.tolist(),
                             nodal_witness=w.tolist(),summed_interpolation_remainder=e.tolist()))
        for s,t in zip(a,b,strict=True):
            (lb,ub),_,_,_=t
            (L,U),_,_,_=s
            if np.any(lb>U) or np.any(ub<L):
                raise ArithmeticError('Disjoint independently meshed enclosures')
        if m==m0:
            for (_,_,obs),((l,h),_,_,_) in zip(charts,b,strict=True):
                # gauss reports the pressure numerator; measure includes its factor 1/3.
                pressure=obs[1]/3.
                if not l[0]-5e-14*pressure<=pressure<=h[0]+5e-14*pressure:
                    raise ArithmeticError('Archived representative outside its own conditional enclosure')
        pairs=[]
        for j,k in ((0,1),(0,2),(1,2)):
            (a0,a1),aw,_,(ai,aj)=b[j]
            (b0,b1),bw,_,(bi,bj)=b[k]
            R=np.array([b0/a1-1,b1/a0-1]).T
            F=np.array([bi/aj-1,bj/ai-1]).T
            if np.any(F[:,0]<R[:,0]-3e-15) or np.any(F[:,1]>R[:,1]+3e-15):
                raise ArithmeticError('Constructive ratio witnesses violate outer enclosure')
            pairs.append(dict(comparison=[lock['charts'][k],lock['charts'][j]],
                              ratio_outer=R.tolist(),ratio_constructively_attained=F.tolist(),
                              nodal_ratio=(bw/aw-1).tolist(),sign_at_fixed_retained_moments=state(R,F)))
        out.append(dict(probe_mass=m,charts=rows,ratios=pairs))
    return dict(x=x,native_measurement_mass=m0,charts=lock['charts'],fields=['pressure','pseudo_pressure'],
                rows=out,input_sha256=[nh,hashlib.sha256(control.read_bytes()).hexdigest(),
                                      hashlib.sha256(raw).hexdigest(),*[h for _,h in D]],
                scope='Positive in-cell measures with fixed archived warm count and mean excitation; cold convention held. Probe masses are not evolution or observational predictions. Ordinary floating-point enclosure.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('native','control','accepted'):
        p.add_argument('--'+k,type=Path,required=True)
    p.add_argument('--states',type=Path,nargs=3,required=True)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
