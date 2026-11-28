"""Sharp worst-case q displacement given two fixed native panel-mass ledgers."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_endpoint_cold_non_nested_cdf_contrast import audit
from early_pole_endpoint_cold_profile_interval_transport import fixed,transport


def most(e,f,p,q):
    p=np.asarray(p,dtype=float)
    q=np.asarray(q,dtype=float)
    e=np.asarray(e,dtype=float)
    f=np.asarray(f,dtype=float)
    if (len(e)!=len(p)+1 or len(f)!=len(q)+1 or min(p)<0. or min(q)<0.
        or min(np.diff(e))<=0. or min(np.diff(f))<=0. or [e[0],e[-1]]!=[f[0],f[-1]]
        or max(abs(fsum(p)-1.),abs(fsum(q)-1.))>3e-14):
        raise ArithmeticError('Require two valid native normalized histograms on one band')
    n=len(p)
    m=len(q)
    A=np.r_[0.,np.cumsum(p)]
    B=np.r_[0.,np.cumsum(q)]
    A[-1]=B[-1]=1.
    A=np.clip(A,0.,1.)
    B=np.clip(B,0.,1.)
    unit=e[-1]-e[0]
    atlas=np.r_[np.stack((e[:-1],e[1:]),axis=1),np.stack((f[:-1],f[1:]),axis=1)]
    forest=[[] for _ in range(n+m)]
    pairs=[]
    def add(i,j,w):
        cost=w*abs(atlas[i,:,None]-atlas[j,None,:])/unit
        forest[i].append((j,cost))
        forest[j].append((i,cost.T))
        pairs.append((i,j,w,cost))
    i=j=0
    while i<n and j<m:
        w=min(A[i+1],B[j+1])-max(A[i],B[j])
        if w>0.:
            add(i,n+j,w)
        a,b=A[i+1],B[j+1]
        i+=a<=b
        j+=b<=a
    if abs(fsum(w for _,_,w,_ in pairs)-1.)>3e-14:
        raise ArithmeticError('Fractional transport path lost mass')
    seen=set()
    total=0.
    choice=np.zeros(n+m,dtype=int)
    for k in range(n+m):
        if k in seen:
            continue
        seen.add(k)
        order=[k]
        children={}
        for i in order:
            arr=[]
            for j,mat in forest[i]:
                if children.get((i,'parent'))==j:
                    continue
                if j in seen:
                    raise ArithmeticError('Ordered quantile-overlap graph is not a forest')
                seen.add(j)
                order.append(j)
                arr.append((j,mat))
                children[(j,'parent')]=i
            children[i]=arr
        DP={}
        for i in order[::-1]:
            DP[i]=sum((np.max(v+DP[j][None,:],axis=1) for j,v in children[i]),start=np.zeros(2))
        choice[k]=int(np.argmax(DP[k]))
        total+=float(DP[k][choice[k]])
        for i in order:
            for j,v in children[i]:
                choice[j]=int(np.argmax(v[choice[i]]+DP[j]))
    actual=fsum(v[choice[i],choice[j]] for i,j,_,v in pairs)
    if abs(total-actual)>2e-14 or min(p.min(),q.min())<0.:
        raise ArithmeticError('Endpoint witness does not reproduce the discrete extremum')
    return dict(sharp_upper_q=total*unit,endpoint_witness_q=actual*unit,
                left_right_endpoint_choices=[choice[:n].tolist(),choice[n:].tolist()],
                number_of_positive_quantile_overlaps=len(pairs))


def run(fine,third,states,archive,fixture):
    a,b=[audit(p,s,archive,g) for p,s,g in zip((fine,third),states,('fine','third'),strict=True)]
    if [a[1],b[1]]!=fixture['states_sha256_uncompressed'] or [a[2],b[2]]!=fixture['gates']:
        raise ArithmeticError('The response ledgers are not the accepted same-band states')
    out=[]
    for (e,A),(f,B) in zip(a[0],b[0],strict=True):
        row={}
        for key,i,j in (('forward_number',0,0),('inverse_number',0,1),('forward_excess',2,0)):
            P,Q=[z/fsum(z) for z in (A[:,i,j],B[:,i,j])]
            lower=transport(fixed(dict(base=P)),fixed(dict(base=Q)),e,f,'fixed',independent_bracket=True)
            upper=most(e,f,P,Q)
            if lower['independent_cost_lower_upper_q'][1]>upper['sharp_upper_q']:
                raise ArithmeticError('Independently feasible minimal cost exceeds sharp maximum')
            row[key]=dict(minimum_cost_independently_enclosed=lower['independent_cost_lower_upper_q'],
                          **upper)
        out.append(row)
    A,B=out
    err=max(abs(A[k]['sharp_upper_q']-z['sharp_upper_q']) for k,z in B.items())
    if err>2e-12:
        raise ArithmeticError('Worst-case cost not independently resolved across q meshes')
    return dict(grids=['fine','third'],fixed_x=.05,charge_subcells=2,
                states_sha256_uncompressed=fixture['states_sha256_uncompressed'],
                observables=B,maximum_paired_upper_q_defect=err,
                scope='Best and worst interval placement at fixed native panel masses; endpoint assignment is a limiting witness, not an observed intrapanel distribution.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('fine','third','states'):
        p.add_argument('--'+k,type=Path,nargs=2,required=True)
    for k in ('archive','fixture'):
        p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args()
    print(json.dumps(run(a.fine,a.third,a.states,a.archive,json.loads(a.fixture.read_text())),indent=2,allow_nan=False))
