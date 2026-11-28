"""Resolve how positive charged-lepton elastic heat changes closed-bath FD modes."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from boost_axis_elastic_heat_reference import integrate as elastic, leptons
from fd_portal_differential_reciprocity_reference import PositivePortal, integrate
from fd_portal_finite_capacity_modes_reference import dark
from fd_portal_reservoir_exponential_filter_reference import symmetric_root
from sm_dark_entropy_exchange_reference import StandardModel


def sample(p,mesh,T,cv):
    y,A=dark(T,0.)
    H=np.diag([1.,1./T])@np.linalg.inv(A)+np.diag([0.,1/(T*T*cv)])
    P,_=symmetric_root(H)
    z,J,_,_=integrate(p,T,0.,0.)
    S=np.diag([2.,T])
    M=S@J@S
    n,cut=mesh
    q=T*sum(elastic(T,T,0.,m,n=n,cut=cut,linear=True) for m in leptons)
    if min(np.diag(M))<=0. or q<=0. or max(abs(z))!=0.:
        raise ArithmeticError('Finite-bath equilibrium dissipation cone not reproduced')
    lo,U=np.linalg.eigh(P@M@P)
    hi,V=np.linalg.eigh(P@(M+np.diag([0.,q]))@P)
    if min(lo)<=0. or np.min(hi/lo)<1.-3e-12:
        raise ArithmeticError('Positive elastic heat cannot reduce either closed-system mode rate')
    rotations=np.arccos(np.clip(abs(np.diag(U.T@V)),0.0,1.0))
    return dict(timelike_modes=lo.tolist(),combined_modes=hi.tolist(),
                combined_over_timelike=(hi/lo).tolist(),elastic_over_timelike_heat_mobility=q/M[1,1],
                mode_rotation_radians=rotations.tolist())


def run(path):
    eos=StandardModel(path)
    keys=('timelike_modes','combined_modes','combined_over_timelike','elastic_over_timelike_heat_mobility')
    p,q=[PositivePortal(o) for o in ((13,36,26),(17,44,32))]
    out=[]
    for T in (.025,.05,.075,.1,.15,.2,.3,.45,.55,.64135,.6882487043,.8,.95):
        cv=eos.state(T)[3]
        a,b=[sample(z,m,T,cv) for z,m in ((p,(18,29.)),(q,(25,36.)))]
        er=max(np.max(abs(np.asarray(a[k])/b[k]-1.)) for k in keys)
        if er>2e-5 or np.max(abs(np.asarray(a['mode_rotation_radians'])-b['mode_rotation_radians']))>2e-5:
            raise ArithmeticError('Elastic mode modification not resolved by independent rules')
        out.append(dict(T_SM=T,independent_mesh_relative_difference=er,**b))
    return dict(rows=out,table_input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                scope='Local conservative bath kinetic modes measured for the restricted vacuum charged-lepton portal. No expansion or WIMP density inference.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('eos',type=Path)
    print(json.dumps(run(p.parse_args().eos),indent=2,allow_nan=False))
