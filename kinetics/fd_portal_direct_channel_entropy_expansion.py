"""Reintegrate the combined direct path with separate entropy ledgers by channel."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from fd_portal_combined_tabulated_bath_entropy import elastic_heat
from fd_portal_differential_reciprocity_reference import PositivePortal
from fd_portal_tabulated_entropy_split_expansion import Stage,hubble
from sm_dark_entropy_exchange_reference import StandardModel


def path(eos,p,n,eta_controls,method,rtol):
    T0=.90
    S=Stage(eos,p,n)
    ini=np.array([T0,0.,np.log(T0)])
    S0=S.entropy(0.,ini)
    mesh=(18,29.) if n==112 else (25,36.)

    def flow(x,y):
        T=y[0]*T0
        e,l=y[1:3]
        Td=T0*np.exp(l)
        t=1./T-1./Td
        f,s=S.collision(T,Td,e)
        q=elastic_heat(T,Td,e,mesh)
        se=t*q
        if min(s,se)<0.:
            raise ArithmeticError('Negative channel entropy')
        f[1]+=q
        v,w,_,cv=eos.state(T)
        d=S.bath([e,np.log(Td)])
        nd,rd,pd,B,E,C=[d[k] for k in ('n','rho','pressure','dn_deta','chemical_energy','cv_fixed_number')]
        A=np.array([[B,B*E/Td],[B*E,Td*C+B*E*E/Td]])
        H=hubble(v,rd)
        return np.r_[(-3*(v+w)+f[1]/H)/(T0*cv),np.linalg.solve(A,-3*np.array([nd,rd+pd])-f/H),
                     np.exp(3*x)*np.array([s,se])/(S0*H)]

    a=np.asarray(eta_controls)
    x=a[:,0]
    solver=solve_ivp(flow,(0.,x[-1]),[1.,0.,0.,0.,0.],method=method,rtol=rtol,atol=rtol*.018,t_eval=x)
    if not solver.success:
        raise ArithmeticError('Channel-split entropy trajectory did not finish')
    h=solver.y.T
    h[:,0]*=T0
    h[:,2]+=np.log(T0)
    diff=np.max(abs(h[:,:3]-a[:,1:4]))
    cs=np.max(abs(np.sum(h[:,3:5],axis=1)-a[:,5]))
    if min(h[-1,3:])<=0. or diff>9e-9 or cs>4e-12:
        raise ArithmeticError('Independent channel integration did not reproduce the combined archive')
    return h,dict(state_archive_defect=diff,total_entropy_archive_defect=cs,evaluations=solver.nfev)


def run(eosfile,reference):
    a=reference.read_bytes()
    d=json.loads(a)
    digest=hashlib.sha256(eosfile.read_bytes()).hexdigest()
    if d['table_input_sha256']!=digest or d['initial_T']!=.90 or d['rows'][1]['operator']!='combined':
        raise ArithmeticError('Require the exact .90 direct two-channel reference')
    eos=StandardModel(eosfile)
    r0,r1=[np.asarray(z['controlled_full_path']) for z in d['rows']]
    p,q=[PositivePortal(o) for o in ((13,36,26),(17,44,32))]
    (a1,b1),(a2,b2)=[path(eos,k,n,r1,m,t) for k,n,m,t in ((p,112,'BDF',4e-9),(q,168,'Radau',4e-10))]
    error=np.max(abs(a1[:,3:]-a2[:,3:]))
    if error>2e-12:
        raise ArithmeticError('Independent entropy integrals disagree')
    return dict(reference_sha256=hashlib.sha256(a).hexdigest(),table_input_sha256=digest,
                fields=['log_a','combined_path_timelike_entropy','combined_path_elastic_entropy',
                        'timelike_path_timelike_entropy','combined_minus_timelike_total_entropy'],
                rows=np.column_stack((r1[:,0],a2[:,3:5],r0[:,5],a2[:,3]+a2[:,4]-r0[:,5])).tolist(),
                maximum_integral_mesh_defect=error,replays=[b1,b2],
                scope='Separate integrated histories; the extra positive elastic term does not imply larger total integrated entropy than a different solution.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('eos',type=Path)
    p.add_argument('tracking',type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.tracking),indent=2,allow_nan=False))
