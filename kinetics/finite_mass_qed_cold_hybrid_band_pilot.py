"""Integrate on a shared scalar mesh with separate charged/dispersion factors."""
import argparse
import json
import os
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre

from finite_mass_qed_cold_sensitivity_pilot import GRIDS,K,load_all,scalar_occup
from finite_mass_qed_future_q_support_pilot import candidates
from finite_mass_qed_timelike_pole_angular_reference import angle


LO,HI=1e-6,1.5e-5
MESH=((2,3,3,8),(4,6,6,12),(8,8,8,16))


def run(charged,geo,mesh,out):
    inputs=load_all()
    archive,ch,_,_,_=inputs[charged]
    _,_,ms,w,R=inputs[geo]
    marks,counts=candidates(archive,ch,K,w,ms,LO,HI)
    both=np.unique(np.r_[marks,*[v for g in GRIDS for v in inputs[g][0]['scalar_momentum_edges'] if LO<v<HI]])
    n,cos,E,az=MESH[mesh]
    x,y=roots_legendre(n)
    used=set()
    tags={'charged_grid':charged,'dispersion_grid':geo,'mesh':list(MESH[mesh]),
          'scalar_grids':list(GRIDS),'total':len(both)-1,
          'k':K,'w':w,'m_scalar':ms}
    if out.exists():
        for obj in map(json.loads,out.read_text(encoding='utf-8').splitlines()):
            j=obj['panel']
            if any(obj[k]!=v for k,v in tags.items()) or j in used or obj['bounds']!=both[j:j+2].tolist():
                raise ArithmeticError('Incompatible hybrid checkpoint')
            used.add(j)
    with out.open('a',encoding='utf-8',buffering=1) as f:
        for j,(l,h) in enumerate(zip(both[:-1],both[1:],strict=True)):
            if j in used:
                continue
            accum=np.zeros((2,3,3))
            for q,weight in zip((h+l)/2+(h-l)*x/2,(h-l)*y/2,strict=True):
                Es=np.hypot(ms,q)
                A,B,D=angle(ch,K,w,R,Es,'T',cos,(E,az),q)[0]
                s=np.array([scalar_occup(inputs[g][0],q) for g in GRIDS])
                stats=np.c_[(1+s)*A,s*B,s*D+A]
                scalar_q_weight=weight*q*q/Es
                accum+=scalar_q_weight*np.array([stats,Es*stats,(q*q/(Es+ms))*stats]).transpose(1,0,2)
            f.write(json.dumps(dict(tags,panel=j,bounds=[l,h],number_scalar_energy_excess=accum.tolist()),
                               allow_nan=False)+'\n')
            f.flush()
            os.fsync(f.fileno())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('charged',choices=GRIDS)
    p.add_argument('geometry',choices=GRIDS)
    p.add_argument('mesh',type=int,choices=(0,1,2))
    p.add_argument('out',type=Path)
    a=p.parse_args()
    run(a.charged,a.geometry,a.mesh,a.out)
