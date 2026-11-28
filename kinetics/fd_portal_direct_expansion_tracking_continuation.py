"""Continue the exact direct expanding FD archives toward the local slow clock."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from fd_portal_combined_entropy_split_expansion import CombinedExpansion
from fd_portal_differential_reciprocity_reference import PositivePortal
from fd_portal_tabulated_entropy_split_expansion import Stage,hubble
from sm_dark_entropy_exchange_reference import StandardModel


def follow(s,old,xout,method,tol):
    T=.90
    original=np.array([T,0.,np.log(T)])
    S,E=s.entropy(0.,original),s.energy(0.,original)
    N=s.bath(original[1:])['n']
    x0,ts,eta,ld,_,ep,w,q=old
    initial=np.array([ts/T,eta,ld-np.log(T),ep,w,q])

    def rhs(x,u):
        Ts,eta,Td=u[0]*T,u[1],T*np.exp(u[2])
        d=s.bath([eta,np.log(Td)])
        rs,ps,_,cv=s.eos.state(Ts)
        n,r,p,B,dE,C=[d[k] for k in ('n','rho','pressure','dn_deta','chemical_energy','cv_fixed_number')]
        H=hubble(rs,r)
        M=np.array([[B,B*dE/Td],[B*dE,Td*C+B*dE*dE/Td]])
        z,sigma=s.collision(Ts,Td,eta)
        du=np.linalg.solve(M,-3*np.array([n,r+p])-z/H)
        a=np.exp(3*x)
        return np.r_[(-3*(rs+ps)+z[1]/H)/(cv*T),du,a*sigma/(H*S),3*a*(ps+p)/E,-a*z[0]/(H*N)]

    t=solve_ivp(rhs,(x0,xout[-1]),initial,t_eval=xout,method=method,rtol=tol,atol=.022*tol)
    if not t.success:
        raise ArithmeticError('Continuation did not finish')
    a=t.y.T.copy()
    a[:,0]*=T
    a[:,2]+=np.log(T)
    resolved=[]
    for x,v in zip(xout,a,strict=True):
        n=s.bath(v[1:3])['n']
        L=[abs(s.entropy(x,v[:3])/S-1.-v[3]),abs(s.energy(x,v[:3])/E-1.+v[4]),abs(np.exp(3*x)*n/N-1.-v[5])]
        if max(L)>2e-7 or v[3]<ep-2e-12:
            raise ArithmeticError('Combined continuation conservation check failed')
        resolved.append(L)
    b=np.column_stack((xout,a[:,:3],1.-a[:,0]/np.exp(a[:,2]),a[:,3:]))
    if max(abs(b[0]-old))>2e-15:
        raise ArithmeticError('Does not preserve the exact initial archive endpoint')
    return b,np.max(resolved,axis=0)


def run(path,lock):
    raw=lock.read_bytes()
    d=json.loads(raw)
    h=hashlib.sha256(path.read_bytes()).hexdigest()
    if d['table_input_sha256']!=h or d['initial_T']!=.90 or [r['operator'] for r in d['rows']]!=['timelike','combined']:
        raise ArithmeticError('Need exact original input and phase-space variants')
    x=np.array([.25,.30,.35,.40,.45,.50,.55])
    eos=StandardModel(path)
    p,q=[PositivePortal(o) for o in ((13,36,26),(17,44,32))]
    z=[]
    for old,St in zip(d['rows'],(Stage,CombinedExpansion),strict=True):
        init=np.array(old['controlled_full_path'][-1])
        a,b=[follow(St(eos,Q,n),init,x,m,t) for Q,n,m,t in ((p,112,'BDF',8e-9),(q,168,'Radau',9e-10))]
        defect=np.max(abs(a[0]-b[0]),axis=0)
        if max(defect)>3e-7:
            raise ArithmeticError('Two independent continuation paths failed to agree')
        z.append(dict(operator=old['operator'],points=b[0].tolist(),ledger_defects=b[1].tolist(),independent_defects=defect.tolist()))
    return dict(archive_sha256=hashlib.sha256(raw).hexdigest(),table_input_sha256=h,fields=d['fields'],rows=z,
                combined_minus_timelike=(np.asarray(z[1]['points'])-z[0]['points']).tolist(),
                scope='Finite numerical continuation of the restricted vacuum charged-lepton tests, not complete early-universe evolution or freeze-out.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('eos',type=Path)
    p.add_argument('archive',type=Path)
    o=p.parse_args()
    print(json.dumps(run(o.eos,o.archive),indent=2,allow_nan=False))
