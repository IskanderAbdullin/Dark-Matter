"""Search the exact shellwise affine blocking inequality at tagged F=1/2."""

import argparse
import itertools as it
import json
from pathlib import Path

import numpy as np

from three_species_decay_onsager_reference import m,mp
from warm_light_occupation_stage_audit_reference import digest,read


def from_archive(path):
    r=read(path)
    r=r[-1] if isinstance(r,list) else r
    S=r["snapshot"]
    assert digest(S)==r["snapshot_sha256"]
    p=np.asarray(S["lab_momentum_edges_GeV"])
    e=[np.hypot(p,z) for z in (mp,mp,m)]
    a,b=[np.asarray(S["entropy_variable_"+k]) for k in ("scalar","fermion")]
    t=[a,a,b]
    N=len(p)-1
    A=np.array(list(it.product(range(N),repeat=3)))
    A=A[A[:,0]<=A[:,1]]
    lo=np.stack([v[A[:,i]] for i,v in enumerate(e)],axis=1)
    hi=np.stack([v[A[:,i]+1] for i,v in enumerate(e)],axis=1)
    offsets=np.stack([x[0,A[:,i]] for i,x in enumerate(t)],axis=1)
    slopes=np.stack([x[1,A[:,i]] for i,x in enumerate(t)],axis=1)/(hi-lo)
    fmin,fmax=e[2][0],e[2][-1]
    sign=np.array([1.,1.,-1.])
    zbest=(np.inf,None)
    feasible=(np.inf,None)
    tolerance=1e-15

    def study(C):
        nonlocal zbest,feasible
        test=(C@sign>=fmin-tolerance)&(C@sign<=fmax+tolerance)
        test &= np.all((C>=lo-tolerance)&(C<=hi+tolerance),axis=1)
        if not np.any(test):return
        u=np.flatnonzero(test)
        Q=C[u]
        z=(offsets[u]+slopes[u]*(Q-lo[u]))@sign
        k=np.argmin(z)
        if z[k]<zbest[0]:
            zbest=(z[k],(u[k],Q[k]))
        P=Q@sign
        pe=np.sqrt(np.maximum(0.,(P-m)*(P+m)))
        qe=np.sqrt(np.maximum(0.,(Q[:,2]-m)*(Q[:,2]+m)))
        smin=2*m*m+2*(P*Q[:,2]-pe*qe)
        smax=2*m*m+2*(P*Q[:,2]+pe*qe)
        W=Q[:,0]+Q[:,1]
        optimal=np.clip(2*mp*W,smin,smax)
        dq2=np.maximum(0.,W*W-optimal)*(1-4*mp*mp/optimal)
        good=(abs(Q[:,0]-Q[:,1])**2<dq2)
        # Ceiling statements distinguish feasible strict interior from a
        # limiting shell-edge witness, where a piecewise profile has two sides.
        k2=np.argmin(np.where(good,z,np.inf))
        if good[k2] and z[k2]<feasible[0]:
            feasible=z[k2],(u[k2],Q[k2])

    for bits in it.product((0,1),repeat=3):
        study(np.where(np.array(bits),hi,lo))
    for j in range(3):
        fixed=[k for k in range(3) if k!=j]
        for bits in it.product((0,1),repeat=2):
            Q=np.zeros_like(lo)
            Q[:,fixed]=np.where(np.array(bits),hi[:,fixed],lo[:,fixed])
            for edge in (fmin,fmax):
                Q[:,j]=(edge-Q[:,fixed]@sign[fixed])/sign[j]
                study(Q)

    def report(x):
        gap,pair=x
        j,q=pair
        et=q@sign
        return dict(gap=float(gap),tagged_fermion_energy_GeV=float(et),shells=A[j].tolist(),
                    outgoing_scalar1_scalar2_partner_energy_GeV=q.tolist())

    return dict(scope="Conservative energy-only cube minimum and a kinematically feasible limiting event. Not an integrated cap.",
                T_GeV=S["T_visible_GeV"],snapshot_sha256=digest(S),outer_p_GeV=float(p[-1]),
                min_energy_only_entropy_gap=report(zbest),min_feasible_among_vertices_entropy_gap=report(feasible))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("paths",nargs="+",type=Path)
    a=p.parse_args()
    print(json.dumps([from_archive(f) for f in a.paths],indent=2,allow_nan=False))
