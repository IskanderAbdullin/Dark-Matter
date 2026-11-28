"""Separate expansion and collision derivatives in the full comoving coordinates."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, root

from fd_portal_combined_entropy_split_expansion import CombinedExpansion
from fd_portal_differential_reciprocity_reference import PositivePortal
from fd_portal_live_fixed_energy_collision_spectrum import fixed_energy, jac
from fd_portal_tabulated_entropy_split_expansion import Stage, hubble
from sm_dark_entropy_exchange_reference import StandardModel


def full(stage, center):
    """Use center-normalized comoving variables, with a fixed diagonal basis at this x."""
    T0,e0,l0=center
    R0=stage.eos.state(T0)[0]
    d0=stage.bath([e0,l0])
    scale=np.array([d0['n'],d0['rho'],R0+d0['rho']])
    def rhs(v):
        N,R,U=scale*v
        if min(N,R,U-R)<=0.:
            raise ArithmeticError('Nonphysical duration-local perturbation')
        def dr(z):
            d=stage.bath(z)
            B,E,C=(d[k] for k in ('dn_deta','chemical_energy','cv_fixed_number'))
            t=np.exp(z[1])
            D=np.array([[B,B*E/t],[B*E,t*C+B*E*E/t]])/np.array([N,R])[:,None]
            return np.array([d['n']/N-1.,d['rho']/R-1.]),D
        z=root(lambda u:dr(u)[0],[e0,l0],jac=lambda u:dr(u)[1],tol=2e-12).x
        T=brentq(lambda q:stage.eos.state(q)[0]/(U-R)-1.,stage.lo,stage.hi,xtol=5e-15)
        H=hubble(U,0.)
        q,_=stage.collision(T,np.exp(z[1]),z[0])
        pv,pd=stage.eos.state(T)[1],stage.bath(z)['pressure']
        if max(abs(dr(z)[0]))>4e-13:
            raise ArithmeticError('Physical density inversion failed to resolve')
        return np.array([[0.,-3*pd,-3*(pv+pd)],[-q[0]/H,-q[1]/H,0.]])/scale
    return rhs


def diff(f,h):
    v=np.ones(3)
    return np.stack([(f(v-2*d)-8*f(v-d)+8*f(v+d)-f(v+2*d))/(12*h) for d in h*np.eye(3)],axis=-1)


def measure(stage, p, h):
    z=p[1:4]
    f=full(stage,z)
    d=diff(f,h)
    earlier=diff(f,2*h)
    fixed,_,_=fixed_energy(stage,z)
    defect=np.max(abs(jac(fixed,h)-d[1,:2,:2]))
    if defect>3e-8 or np.max(abs(d-earlier))>5e-7:
        raise ArithmeticError('Full comoving Jacobian does not recover the collision restriction')
    return d,defect,np.max(abs(d-earlier))


def run(eos_file, earlier, later):
    raw=[p.read_bytes() for p in (earlier,later)]
    A,B=map(json.loads,raw)
    h=hashlib.sha256(eos_file.read_bytes()).hexdigest()
    if h!=A['table_input_sha256'] or h!=B['table_input_sha256'] or hashlib.sha256(raw[0]).hexdigest()!=B['previous_segment_sha256']:
        raise ArithmeticError('Need exact joined histories')
    E=StandardModel(eos_file)
    P,Q=[PositivePortal(m) for m in ((13,36,26),(17,44,32))]
    results=[]
    for K,r,s in zip((Stage,CombinedExpansion),A['rows'],B['rows'],strict=True):
        if r['operator']!=s['operator']:
            raise ArithmeticError('History operator mismatch')
        for p in [*r['points'][-4::1],*s['points'][4::4]]:
            a,b=[measure(K(E,m,n),p,h0) for m,n,h0 in ((P,112,6e-4),(Q,168,3e-4))]
            delta=np.max(abs(a[0]-b[0]))
            if delta>9e-6:
                raise ArithmeticError('Independent meshes failed to resolve the full comoving Jacobian')
            d=b[0]
            J=d.sum(axis=0)
            e=np.linalg.eigvals(J)
            results.append(dict(operator=r['operator'],x=p[0],state=p[1:4],expansion_jacobian=d[0].tolist(),
                      collision_jacobian=d[1].tolist(),full_jacobian=J.tolist(),local_eigen_real=e.real.tolist(),
                      local_eigen_imaginary=e.imag.tolist(),absolute_defects=[b[1],b[2],delta]))
    return dict(rows=results,variables=['comoving_dark_number','comoving_dark_energy','comoving_total_energy'],
                defect_fields=['restricted_collision','finite_difference','independent_mesh'],
                input_sha256=[hashlib.sha256(v).hexdigest() for v in raw],table_input_sha256=h,
                scope='Local Jacobians in full comoving physical coordinates; no stability statement for the time-ordered expanding flow or a physical threshold.')


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    for k in ('eos_file','earlier','later'):
        ap.add_argument(k,type=Path)
    print(json.dumps(run(**vars(ap.parse_args())),indent=2,allow_nan=False))
