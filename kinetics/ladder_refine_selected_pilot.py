"""Retry incomplete scalar-momentum points with componentwise refinement."""
import argparse
import json
from pathlib import Path

import numpy as np

from finite_mass_qed_frozen_transverse_plasmon_gate_reference import Q
from finite_mass_qed_timelike_pole_angular_reference import angle
from finite_mass_qed_timelike_pole_batched_absorption_reference import load,pole


def resolved(a,b):
    # A/B are separate positive channels. Fixing a violation by
    # comparing only to max(|A|,|B|) can bury an inverse-channel error.
    den=np.maximum(np.maximum(abs(a),abs(b)),1e-25)
    err=abs(a-b)/den
    return np.max(err)<3e-6,err


def run(g,base,result):
    retained={}
    for f in base:
        for s in f.read_text().splitlines():
            x=json.loads(s)
            key=(x['pol'],x['k'],x['q'])
            if key in retained or x['grid']!=g:
                raise ArithmeticError(('Wrong or repeated input panel',f,key))
            retained[key]=x
    endpoint,ch,_=load(g)
    m=endpoint['checkpoints'][-1]['scalar_pole_over_T']
    bins=np.array(endpoint['scalar_momentum_edges'])
    occ=endpoint['endpoint_scalar_occupations']
    points=(0.,Q/2,Q,Q*2,4e-5,1e-4,3e-4,8e-4,.0012)
    modes=(('T',.0015),('T',.00256),('L',.0012),('L',.00256))
    targets=[(p,k,q) for p,k in modes for q in points if (p,k,q) not in retained]
    report=[]
    for p,k,q in targets:
        w,R=pole(ch,p,k)
        E=np.hypot(m,q)
        if not E<w:
            raise ArithmeticError(('Past-directed energy assumption failed',g,p,k,q))
        if q>=.0012:
            meshes=((16,(24,36)),(20,(32,48)),(26,(38,56)))
        else:
            meshes=((6,(6,12)),(10,(10,20)),(16,(16,28)))
        previous=None
        attempts=[]
        for mesh in meshes:
            X,cnt,n,topology=angle(ch,k,w,R,E,p,*mesh,q)
            if not np.all(np.isfinite(X)) or max(abs(X))==0 or np.min(X[:2])<0:
                raise ArithmeticError(('Unexpected signed or empty channel',g,p,k,q,X))
            if abs(X[0]-X[1]-X[2])>5e-13*max(abs(X)):
                raise ArithmeticError(('Affine ledger has lost separated forward/inverse rates',g,p,k,q,X))
            attempts.append(dict(mesh=mesh,raw=X.tolist(),panels=n,angle_error=topology))
            if previous is not None:
                ok,err=resolved(previous,X)
                attempts[-1].update(componentwise_defect=err.tolist())
                if ok:
                    f=occ[np.searchsorted(bins,q)-1] if q>0. else None
                    threshold=X[1]/(X[2]+X[0]/f) if f and (X[2]+X[0]/f)>0. else None
                    report.append(dict(grid=g,pol=p,k=k,q=q,E_scalar=E,omega=w,scalar_occup=f,
                                       raw_A_B_D=X.tolist(),isotropic_conditional_population=threshold,
                                       resolved_componentwise=True,attempts=attempts))
                    break
            previous=X
        else:
            report.append(dict(grid=g,pol=p,k=k,q=q,resolved_componentwise=False,attempts=attempts))
        with result.open('a',encoding='utf-8') as f:
            f.write(json.dumps(report[-1],allow_nan=False)+'\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('grid',choices=('fine','third'))
    parser.add_argument('result',type=Path)
    parser.add_argument('previous',nargs='+',type=Path)
    o=parser.parse_args()
    if o.result.exists():
        raise FileExistsError(o.result)
    run(o.grid,o.previous,o.result)
