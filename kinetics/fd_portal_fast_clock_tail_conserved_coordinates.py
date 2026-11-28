"""Independent comoving-coordinate check of the restricted FD fast-clock tail."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq,root

from fd_portal_combined_entropy_split_expansion import CombinedExpansion
from fd_portal_differential_reciprocity_reference import PositivePortal
from fd_portal_tabulated_entropy_split_expansion import Stage,hubble
from sm_dark_entropy_exchange_reference import StandardModel


def integrate(stage,old,xout=None,*,maximum_step=.002):
    Ts=.90
    root_case=np.array([Ts,0.,np.log(Ts)])
    E,S=stage.energy(0.,root_case),stage.entropy(0.,root_case)
    N=stage.bath(root_case[1:])['n']
    x0,Vi,e0,l0,_,entropy,work,number=old
    d0=stage.bath([e0,l0])
    D=d0['rho']*np.exp(3*x0)
    B=stage.energy(x0,old[1:4])
    u0=[np.exp(3*x0)*d0['n']/N,D/E,B/E,entropy,work]

    def state(x,y):
        a=np.exp(-3*x)
        n,r=(N*y[0]*a,E*y[1]*a)
        s=E*y[2]*a-r
        if min(n,r,s)<=0.:
            raise ArithmeticError('Unphysical extensive stage')
        T=brentq(lambda x:stage.eos.state(x)[0]/s-1.,stage.lo,stage.hi,xtol=5e-15)
        def model(z):
            v=stage.bath(z)
            B,C,d=(v[k] for k in ('dn_deta','cv_fixed_number','chemical_energy'))
            td=np.exp(z[1])
            J=np.array([[B,B*d/td],[B*d,td*C+B*d*d/td]])/np.array([n,r])[:,None]
            return (np.array([v['n'],v['rho']])-[n,r])/[n,r],J
        a=root(lambda z:model(z)[0],[e0,l0+.25*np.log(r/d0['rho'])],
               jac=lambda z:model(z)[1],tol=2e-12).x
        if np.max(abs(model(a)[0]))>4e-13:
            raise ArithmeticError('Independent extensive dark inverse failed')
        return np.r_[T,a]

    def rhs(x,y):
        T,e,q=state(x,y)
        Td=np.exp(q)
        es,ps,_,_=stage.eos.state(T)
        d=stage.bath([e,q])
        ed,pd=[d[k] for k in ('rho','pressure')]
        H=hubble(es,ed)
        Q,P=stage.collision(T,Td,e)
        a=np.exp(3*x)
        work=3*a*(ps+pd)/E
        return [-a*Q[0]/(N*H),-a*(3*pd+Q[1]/H)/E,-work,a*P/(S*H),work]

    x=np.array([.91,.93,.96]) if xout is None else np.asarray(xout)
    if x[0]!=x0 or np.any(np.diff(x)<=0.) or maximum_step<=0.:
        raise ValueError('Need strictly ordered outputs starting at the exact archived endpoint')
    t=solve_ivp(rhs,(x[0],x[-1]),u0,method='Radau',atol=1e-13,rtol=4e-12,max_step=maximum_step,t_eval=x)
    if not t.success:
        raise ArithmeticError('Independent extensive Radau solve failed')
    rows=[]
    for p,z in zip(x,t.y.T,strict=True):
        a=state(p,z)
        h=[abs(stage.energy(p,a)/E-1.+z[4]),abs(stage.entropy(p,a)/S-1.-z[3])]
        if max(h)>2e-7 or z[3]<entropy-1e-13:
            raise ArithmeticError(('Extensive clock drift ledger failure',p,h))
        rows.append(dict(log_a=p,state=a.tolist(),comoving_number_ratio=z[0],relative_entropy_added=z[3],
                         normalized_work=z[4],entropy_energy_absolute_defects=h[::-1]))
    if np.max(abs(np.asarray(rows[0]['state'])-old[1:4]))>2e-13:
        raise ArithmeticError('Independent intensive quantities did not reproduce the archived start')
    return rows


def run(path,lock):
    b=lock.read_bytes(); d=json.loads(b); h=hashlib.sha256(path.read_bytes()).hexdigest()
    if d['table_input_sha256']!=h or [x['operator'] for x in d['rows']]!=['timelike','combined']:
        raise ArithmeticError('Need exact original physical histories and EOS input')
    eos=StandardModel(path)
    p=PositivePortal((17,44,32))
    rows=[]
    for old,Q in zip(d['rows'],(Stage,CombinedExpansion),strict=True):
        v=np.asarray(old['points'][-1])
        if v[0]!=.91:
            raise ArithmeticError('Wrong archived restart')
        rows.append(dict(operator=old['operator'],rows=integrate(Q(eos,p,168),v)))
    return dict(rows=rows,previous_archive_sha256=hashlib.sha256(b).hexdigest(),table_input_sha256=h,
                scope='Independent extensive numerical continuation of restricted FD histories, not a full reaction model or cosmological threshold.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('eos',type=Path)
    p.add_argument('archive',type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.archive),indent=2,allow_nan=False))
