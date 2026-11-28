"""Time-order the actual expanding FD Jacobian in a single comoving basis."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline

from fd_portal_combined_entropy_split_expansion import CombinedExpansion
from fd_portal_differential_reciprocity_reference import PositivePortal
from fd_portal_direct_comoving_variational_jacobian import diff,full
from fd_portal_tabulated_entropy_split_expansion import Stage
from sm_dark_entropy_exchange_reference import StandardModel


def quantities(stage,p,x0,start):
    x,T,e,l=p[:4]
    d=stage.bath([e,l])
    return np.exp(3*(x-x0))*np.array([d['n'],d['rho'],d['rho']+stage.eos.state(T)[0]])/start


def solve(stage, points, tol, step):
    x0,T,e,l=points[0][:4]
    x1=points[-1][0]
    dark=stage.bath([e,l])
    start=np.array([dark['n'],dark['rho'],dark['rho']+stage.eos.state(T)[0]])
    z=full(stage,[T,e,l])
    def rates(x,v):
        a=np.exp(3*(x-x0))
        return a*z(v/a)
    def rhs(x,v):
        return rates(x,v).sum(axis=0)
    def nonlinear(initial):
        r=solve_ivp(rhs,(x0,x1),initial,method='DOP853',atol=tol*.03,rtol=tol,max_step=step,dense_output=True)
        if not r.success or np.min(r.y)<=0.:
            raise ArithmeticError('Unable to follow physical comoving trajectory')
        return r.sol
    control=nonlinear(np.ones(3))
    relative_end=np.array([quantities(stage,p,x0,start) for p in points])
    defect=np.max(abs(relative_end-control(np.array(points)[:,0]).T))
    if defect>5e-8:
        raise ArithmeticError('Unperturbed comoving trajectory does not recover exact saved solution')
    return rates,nonlinear,control,defect


def ordered(rates,control,x,eps):
    tables=[]
    for t in x:
        y=control(t)
        local=diff(lambda v:rates(t,v*y)/y,eps)
        tables.append(local*y[None,:,None]/y[None,None,:])
    tables=np.asarray(tables)
    A=CubicSpline(x,tables)
    def ode(t,z):
        M=A(t).sum(axis=0)
        return np.r_[(M@z[:9].reshape(3,3)).reshape(-1),np.trace(M)]
    u=solve_ivp(ode,(x[0],x[-1]),np.r_[np.eye(3).ravel(),0.],
                method='DOP853',rtol=2e-11,atol=2e-13,max_step=.008)
    if not u.success:
        raise ArithmeticError('Ordered variational solve did not finish')
    P=u.y[:9,-1].reshape(3,3)
    sign,ln=np.linalg.slogdet(P)
    check=abs(ln-u.y[-1,-1])
    if sign!=1 or check>4e-10:
        raise ArithmeticError('Violated Liouville volume closure')
    return P,check,float(u.y[-1,-1])


def run(eos_path,path,direct,operator="combined"):
    raw=[x.read_bytes() for x in (path,direct)]
    D,A=map(json.loads,raw)
    h=hashlib.sha256(eos_path.read_bytes()).hexdigest()
    if h!=D['table_input_sha256'] or h!=A['table_input_sha256'] or D['previous_segment_sha256']!=hashlib.sha256(raw[1]).hexdigest():
        raise ArithmeticError('Must use exact direct source histories and unchanged EOS')
    data=next(r['points'] for r in D['rows'] if r['operator']==operator)
    previous=next(r['points'][-1] for r in A['rows'] if r['operator']==operator)
    if data[0][0]!=.96 or data[-1][0]!=1.20 or max(abs(np.array(data[0])-previous))>2e-15:
        raise ArithmeticError('Expected selected exact archived branch')
    eos=StandardModel(eos_path)
    K={'combined':CombinedExpansion,'timelike':Stage}[operator]
    low,high=[K(eos,PositivePortal(q),n) for q,n in (((13,36,26),112),((17,44,32),168))]
    field0,_,control0,ref0=solve(low,data,7e-11,.021)
    field1,nonlinear,control1,ref1=solve(high,data,1e-11,.014)
    x0,x1=data[0][0],data[-1][0]
    coarse,_,_=ordered(field0,control0,np.linspace(x0,x1,13),6e-4)
    median,_,_=ordered(field1,control1,np.linspace(x0,x1,25),3e-4)
    fine,trace_error,trace=ordered(field1,control1,np.linspace(x0,x1,49),3e-4)
    n1=np.max(abs(coarse-median)); n2=np.max(abs(median-fine))
    if max(n1,n2)>5e-5:
        raise ArithmeticError(('Sampled tangent did not resolve',n1,n2))
    eps=3e-4
    v=np.ones(3)
    direct=np.column_stack([(nonlinear(v+dv)(x1)-nonlinear(v-dv)(x1))/(2*eps) for dv in eps*np.eye(3)])
    residual=np.max(abs(fine-direct))
    if residual>5e-6:
        raise ArithmeticError(('Ordered prediction disagrees with separate nonlinear perturbations',residual))
    initial_basis=np.linalg.svd(fine,compute_uv=False)
    final_basis=np.linalg.svd(fine/control1(x1)[:,None],compute_uv=False)
    return dict(interval=[x0,x1],operator=operator,fields=['comoving_number','comoving_dark_energy','comoving_total_energy'],
                ordered_response=fine.tolist(),finite_perturbation_response=direct.tolist(),
                mesh_defects=[n1,n2],separate_nonlinear_max_defect=residual,bg_defects=[ref0,ref1],
                liouville_abs_error=trace_error,integrated_trace=trace,
                final_baseline_in_start_units=control1(x1).tolist(),singular_values_start_units=initial_basis.tolist(),
                singular_values_endpoint_fractions=final_basis.tolist(),
                input_sha256=[hashlib.sha256(v).hexdigest() for v in raw],table_input_sha256=h,
                scope='Finite-interval tangent of direct expanding restricted FD dynamics. Gains depend on chosen dimensionless comoving norm; not a complete interval or physical stability criterion.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('eos_path','path','direct'):
        p.add_argument(key,type=Path)
    p.add_argument('--operator',choices=('timelike','combined'),default='combined')
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
