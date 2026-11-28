"""Propagate a local linearized FD mode while its temperature changes along a stored path."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline,PchipInterpolator

from boost_axis_elastic_heat_reference import integrate as scatter,leptons
from fd_number_heat_state_reference import evaluate
from fd_portal_differential_reciprocity_reference import PositivePortal,integrate
from fd_portal_tabulated_entropy_split_expansion import hubble
from sm_dark_entropy_exchange_reference import StandardModel


def matrix(eos,p,T,n,combined):
    d=evaluate(T,0.,nodes=n)
    rs,_,_,cs=eos.state(T)
    B,E,C=[d[k] for k in ('dn_deta','chemical_energy','cv_fixed_number')]
    A=np.array([[B,B*E/T],[B*E,T*C+B*E*E/T]])
    H=np.diag([1.,1./T])@np.linalg.inv(A)+np.diag([0.,1/(cs*T*T)])
    _,J,_,_=integrate(p,T,0.,0.)
    S=np.diag([2.,T])
    M=S@J@S
    if combined:
        k,cut=(18,29.) if n==112 else (25,36.)
        M[1,1]+=T*sum(scatter(T,T,0.,m,n=k,cut=cut,linear=True) for m in leptons)
    # All eigenvectors are reported in the same fixed force scaling.
    D=np.array([1.,.70])
    return D[:,None]*(H@M)/D[None,:]/hubble(rs,d['rho'])


def modes(A):
    d,V=np.linalg.eig(A)
    if not np.isrealobj(d) or min(d)<=0.:
        raise ArithmeticError('No resolved positive decay modes')
    j=np.argsort(d)
    return d[j],V[:,j]


def follow(path,eos,p,n,combined,N):
    q=np.asarray(path)
    if q[0,0]!=.25 or q[-1,0]!=.55 or np.any(np.diff(q[:,1])>=0.):
        raise ArithmeticError('Expected the exact sampled cooling trajectory')
    T=PchipInterpolator(q[:,0],q[:,1])
    gx=np.linspace(.25,.55,N)
    y=CubicSpline(gx,np.array([matrix(eos,p,float(T(x)),n,combined) for x in gx]),axis=0)
    _,W=modes(y(.25))
    keys=np.array([.25,.30,.40,.50,.55])
    def f(x,z):
        M=y(x)
        d,_=modes(M)
        return np.r_[-M@z[:2],d]
    v=solve_ivp(f,(keys[0],keys[-1]),np.r_[W[:,0],0.,0.],t_eval=keys,method='Radau',rtol=2e-10,atol=5e-13)
    if not v.success:
        raise ArithmeticError('No resolved ordered tangent diagnostic')
    out=[]
    for x,a in zip(keys,v.y.T,strict=True):
        _,V=modes(y(x))
        # Fix only the sign conventions, retaining each unit Euclidean norm.
        for j in (0,1):
            V[:,j]*=np.sign(V[:,j]@W[:,j])
        value=np.linalg.solve(V,a[:2])
        out.append(dict(log_a=x,ordered_slow_and_fast_components=value.tolist(),local_mode_integrals=a[2:].tolist(),
                        independent_scalar_slow_survival=np.exp(-a[2])))
    return out


def run(path,archive):
    b=archive.read_bytes(); d=json.loads(b)
    h=hashlib.sha256(path.read_bytes()).hexdigest()
    if d['table_input_sha256']!=h or [r['operator'] for r in d['rows']]!=['timelike','combined']:
        raise ArithmeticError('Need the original separated finite histories')
    eos=StandardModel(path)
    p,q=[PositivePortal(o) for o in ((13,36,26),(17,44,32))]
    out=[]
    for j,r in enumerate(d['rows']):
        a,z=[follow(r['points'],eos,o,n,bool(j),N) for o,n,N in ((p,112,33),(q,168,65))]
        rms=max(np.max(abs(np.asarray(x['ordered_slow_and_fast_components'])-y['ordered_slow_and_fast_components']))
                for x,y in zip(a,z,strict=True))
        if rms>3e-6:
            raise ArithmeticError('Ordered propagation did not resolve with independent meshes')
        out.append(dict(operator=r['operator'],max_mesh_and_interpolation_component_difference=rms,rows=z))
    return dict(rows=out,source_history_sha256=hashlib.sha256(b).hexdigest(),table_input_sha256=h,
                scope='Ordered local linear diagnostic on temperature interpolated from direct restricted histories, not a nonlinear perturbation equation or abundance result.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('eos','archive'):
        p.add_argument(k,type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.archive),indent=2,allow_nan=False))
