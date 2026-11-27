"""Compare finite cold-window accumulation and the optional occupation snapshots."""
import argparse
import gzip
import json
from pathlib import Path

import numpy as np
from scipy.special import expit

from early_pole_live_joint_hard_pair_packets_reference import gauss_each


def load(path):
    p=Path(path)
    b=p.read_bytes()
    return json.loads(gzip.decompress(b) if p.suffix=='.gz' else b)


def discrepancy(A,B,a,b):
    edges=[]
    for grids in (A['hard_edge_sets'],B['hard_edge_sets']):
        start=0
        blocks=[]
        for g in grids:
            n=len(g)-1
            blocks.append((start,np.asarray(g)))
            start+=n
        edges.append(blocks)
    lhs,rhs=[np.asarray(z['hard_logits']) for z in (a,b)]
    rows=[]
    for j,mass in enumerate((0.,.05,.003)):
        ai,ea=edges[0][j]; bi,eb=edges[1][j]
        lo=max(ea[0],eb[0]); hi=min(ea[-1],eb[-1])
        for label,l,u in ([('full',lo,hi)] if j==0 else [('full',lo,hi),('aperture',max(lo,1.39),min(hi,1.48))]):
            if l>=u: continue
            es=np.unique(np.r_[l,u,ea[(ea>l)&(ea<u)],eb[(eb>l)&(eb<u)]])
            q,w=gauss_each(es,16)
            def reconstruct(F,i,e):
                j0=np.searchsorted(e,q)-1
                elo,ehi=np.hypot(e[j0],mass),np.hypot(e[j0+1],mass)
                t=2*(np.hypot(q,mass)-.5*(elo+ehi))/(ehi-elo)
                x=F[i+j0,0]+F[i+j0,1]*t
                return (1/np.expm1(x)) if j==0 else expit(x)
            f,g=reconstruct(lhs,ai,ea),reconstruct(rhs,bi,eb)
            measure=w*q*q
            L1=np.sum(measure*abs(f-g))/(measure@(.5*(f+g)))
            rows.append(dict(species=('photon','parent','daughter')[j],region=label,
                             occupation_relative_L1=float(L1)))
    return rows


def compare(A,B):
    if A['duration']!=B['duration'] or A['steps']!=B['steps']:
        raise ValueError('Use matched time intervals and the same macro grid for source increments.')
    u=np.asarray(A['state_history']);v=np.asarray(B['state_history'])
    if not len(u)==len(v)==A['steps']+1:
        raise ValueError('Inconsistent time series.')
    dt=A['duration']/A['steps']
    du,dv=np.diff(u)/dt,np.diff(v)/dt
    idx=sorted({0,1,2,3,len(du)//4-1,len(du)//2-1,3*len(du)//4-1,len(du)-1} & set(range(len(du))))
    mem=[dict(after_step=i+1,exposure=(i+1)*dt,average_last_interval_rate=[float(du[i]),float(dv[i])],
              accumulation_difference=float(u[i+1]-u[0]-v[i+1]+v[0])) for i in idx]
    gap=abs(du-dv)/np.maximum(abs(du),1e-30)
    closure={}
    for eps in (.01,.001):
        bad=np.flatnonzero(gap>=eps)
        at=int(bad[-1]+1) if len(bad) else 0
        closure[str(eps)]=None if at==len(gap) else float(at*dt)
    out=dict(first_time_after_which_every_rate_gap_below=closure,exposure=A['duration'],steps=A['steps'],gains=[float(u[-1]-u[0]),float(v[-1]-v[0])],
             largest_sample_relative_rate_gap=float(np.max(abs(du-dv)/np.maximum(abs(du),1e-30))),
             minimum_sample_rate=[float(np.min(du)),float(np.min(dv))],history=mem)
    if A.get('snapshots') and B.get('snapshots'):
        pa={z['step']:z for z in A['snapshots']}
        pb={z['step']:z for z in B['snapshots']}
        out['occupation']=[dict(step=k,rows=discrepancy(A,B,pa[k],pb[k])) for k in sorted(pa.keys()&pb.keys())]
    return out


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('first')
    p.add_argument('second')
    a=p.parse_args()
    print(json.dumps(compare(load(a.first),load(a.second)),indent=2,allow_nan=False))
