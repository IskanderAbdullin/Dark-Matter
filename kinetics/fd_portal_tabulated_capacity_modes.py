"""Compare local closed-bath FD mode rates using a tabulated heat capacity."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from fd_portal_differential_reciprocity_reference import PositivePortal, integrate
from fd_portal_reservoir_exponential_filter_reference import symmetric_root
from sm_dark_entropy_exchange_reference import StandardModel

# Independent finite capacity reference from the same diagnostic benchmark.
from fd_portal_finite_capacity_modes_reference import dark


def spectrum(p,T,cv):
    _,A=dark(T,0.)
    H=np.diag([1.,1./T])@np.linalg.inv(A)
    _,J,_,_=integrate(p,T,0.,0.)
    S=np.diag([2.,T])
    P,_=symmetric_root(S@J@S)
    result=np.array([np.linalg.eigvalsh(P@(H+np.diag([0.,v]))@P) for v in (
        0.,1./(T*T*cv),1./(4*(10.75*np.pi**2/30.)*T**5))])
    if min(result[0])<=0 or np.min(result[1:]/result[0])<1.-3e-12:
        raise ArithmeticError('Positive heat-capacity addenda must not slow either mode')
    return result


def run(path):
    eos=StandardModel(path)
    a,b=[PositivePortal(o) for o in ((13,36,26),(17,44,32))]
    rows=[]
    for T in (.025,.05,.075,.1,.15,.2,.3,.45,.55,.64135,.6882487043,.8,.95):
        rho,press,s,cv=eos.state(T)
        offsets=T*np.array([-2.,-1.,1.,2.])*8e-6
        independent=sum(k*eos.state(T+q)[0] for k,q in zip((1,-8,8,-1),offsets,strict=True))/(12*8e-6*T)
        if abs(independent/cv-1.)>1e-6:
            raise ArithmeticError('Tabulated heat capacity does not differentiate its own energy EOS')
        levels=[spectrum(p,T,cv) for p in (a,b)]
        err=np.max(abs(levels[0]/levels[1]-1.))
        if err>2e-6:
            raise ArithmeticError('Finite-capacity modes are not resolved by both phase-space meshes')
        z=levels[1]
        rows.append(dict(T_SM=T,effective_ideal_dof_from_heat_capacity=cv/(4*np.pi**2/30*T**3),
                         capacity_derivative_relative_defect=abs(independent/cv-1.),
                         rates_by_reservoir=z.tolist(),rates_relative_to_fixed=(z[1:]/z[0]).tolist(),
                         independent_mesh_relative_difference=err))
    return dict(rows=rows,cases=['infinite_capacity','tabulated_finite','ideal_10p75_finite'],
                table_input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                scope='Near-equilibrium timelike vacuum FD portal with an entropic tabulated reservoir; not a decoupling or abundance prediction.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('eos',type=Path)
    print(json.dumps(run(p.parse_args().eos),indent=2,allow_nan=False))
