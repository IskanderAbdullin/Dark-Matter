"""Direct stiff tracking through a finite, tabulated homogeneous expansion interval."""
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


def integrate(stage,method,rtol,initial_T,xvalues):
    a=np.array([initial_T,0.,np.log(initial_T)])
    S0,E0=stage.entropy(0.,a),stage.energy(0.,a)
    N0=stage.bath(a[1:])['n']

    def law(x,z):
        T=initial_T*z[0]
        eta,l=z[1:3]
        Td=np.exp(l)*initial_T
        rd=stage.bath([eta,np.log(Td)])
        rv,pv,_,cv=stage.eos.state(T)
        n,r,p,B,d,C=(rd[q] for q in ('n','rho','pressure','dn_deta','chemical_energy','cv_fixed_number'))
        J=np.array([[B,B*d/Td],[B*d,Td*C+B*d*d/Td]])
        H=hubble(rv,r)
        v,ep=stage.collision(T,Td,eta)
        q=np.linalg.solve(J,-3*np.array([n,r+p])-v/H)
        push=np.exp(3*x)
        return np.r_[(-3*(rv+pv)+v[1]/H)/(cv*initial_T),q,
                     push*ep/(H*S0),3*(pv+p)*push/E0,-v[0]*push/(H*N0)]

    t=solve_ivp(law,(0.,xvalues[-1]),[1.,0.,0.,0.,0.,0.],t_eval=xvalues,
                method=method,rtol=rtol,atol=rtol*.02)
    if not t.success:
        raise ArithmeticError('Direct stiff solve failed')
    out=[]
    local=[]
    for x,y in zip(t.t,t.y.T,strict=True):
        a=np.array([initial_T*y[0],y[1],y[2]+np.log(initial_T)])
        n=stage.bath(a[1:])['n']
        b=[abs(stage.entropy(x,a)/S0-1.-y[3]),abs(stage.energy(x,a)/E0-1.+y[4]),
           abs(np.exp(3*x)*n/N0-1.-y[5])]
        if max(b)>2e-7 or y[3]<-2e-12:
            raise ArithmeticError('Direct collision and expansion balances did not close')
        z=[x,*a.tolist(),1.-a[0]/np.exp(a[2]),*y[3:6]]
        out.append(z)
        local.append(b)
    return np.array(out),np.max(local,axis=0),t.nfev


def run(path):
    eos=StandardModel(path)
    stages=[(Stage, 'timelike'),(CombinedExpansion,'combined')]
    x=np.array([0.,.025,.05,.10,.15,.20,.25])
    out=[]
    for kind,name in stages:
        p,q=[PositivePortal(o) for o in ((13,36,26),(17,44,32))]
        a,b=[integrate(kind(eos,z,n),method,rtol,.90,x)
             for z,n,method,rtol in ((p,112,'BDF',7e-9),(q,168,'Radau',8e-10))]
        error=np.max(abs(a[0]-b[0]),axis=0)
        if max(error)>3e-7:
            raise ArithmeticError('Independent stiffness algorithm and phase-space mesh do not resolve the trajectory')
        out.append(dict(operator=name,controlled_full_path=b[0].tolist(),max_independent_absolute_defect=error.tolist(),
                        max_entropy_energy_number_ledger=b[1].tolist(),evaluation_counts=[a[2],b[2]]))
    d=np.asarray(out[1]['controlled_full_path'])-out[0]['controlled_full_path']
    return dict(initial_T=.90,fields=['log_a','T_SM','eta_D','log_T_D','beta','comoving_entropy_normalized',
                                    'normalized_pressure_work','normalized_count_source'],rows=out,combined_minus_timelike=d.tolist(),
                table_input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                scope='Finite restricted homogeneous vacuum charged-lepton history on 0<=log(a/a_start)<=0.25. No complete interactions or final abundance.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('eos',type=Path)
    print(json.dumps(run(p.parse_args().eos),indent=2,allow_nan=False))
