"""Track closed-system equilibrium mode rates as the visible heat capacity varies."""
import json

import numpy as np

from fd_number_heat_state_reference import evaluate
from fd_portal_differential_reciprocity_reference import PositivePortal, integrate
from fd_portal_reservoir_exponential_filter_reference import symmetric_root


def dark(T,eta):
    d=evaluate(T,eta,nodes=168)
    B,E,C=[d[k] for k in ('dn_deta','chemical_energy','cv_fixed_number')]
    return np.array([d['n'],d['rho']]),np.array([[B,B*E/T],[B*E,T*C+B*E*E/T]])


def sample(p,T):
    y,A=dark(T,0.)
    # Force g=(eta,1/Ts-1/Td); non-infinite visible C adds a rank-one
    # increment to the entropy Hessian at fixed total energy.
    fixed=np.diag([1.,1./T])@np.linalg.inv(A)
    if not np.allclose(fixed,fixed.T,rtol=3e-14,atol=0.):
        raise ArithmeticError('Reciprocal dark entropy Hessian required')
    _,J,_,_=integrate(p,T,0.,0.)
    S=np.diag([2.,T])
    M=S@J@S
    P,_=symmetric_root(M)
    base=np.linalg.eigvalsh(P@fixed@P)
    ans=[]
    prev=np.ones(2)
    for g in (10750.,107.5,10.75,3.):
        E=g*np.pi**2/30*T**4+y[1]
        cv=4*g*np.pi**2/30*T**3
        H=fixed+np.diag([0.,1/(cv*T*T)])
        d=np.linalg.eigvalsh(P@H@P)
        if min(d)<=0 or np.any(d/base<prev-2e-11):
            raise ArithmeticError('Positive rank-one heat-capacity correction did not order the eigenvalues')
        grad=[]
        # Check df/d(eta,ln(Td)) by centered displacements along isoenergy.
        def grad_u(eta,l):
            Td=T*np.exp(l)
            w,_=dark(Td,eta)
            Ts=((E-w[1])/(g*np.pi**2/30))**.25
            return np.array([eta,1/Ts-1/Td])
        for eps in (2e-5,1e-5):
            D=np.column_stack([(grad_u(*v)-grad_u(*(-v)))/(2*eps) for v in eps*np.eye(2)])
            grad.append(np.max(abs((D-H@A)/np.maximum(np.max(abs(H@A),axis=1)[:,None],1e-11))))
        if grad[-1]>3e-6:
            raise ArithmeticError('Independent isoenergy force derivative does not match predicted Hessian')
        prev=d/base
        ans.append(dict(ideal_visible_dof=g,relative_rates=prev.tolist(),slow_and_fast_rates=d.tolist(),
                        finite_difference_defect=grad))
    return dict(T_SM=T,infinite_capacity_slow_and_fast_rates=base.tolist(),finite_capacity=ans)


def run():
    p,q=[PositivePortal(o) for o in ((13,36,26),(17,44,32))]
    rows=[]
    for T in (.05,.3,.64135,.6882487043):
        a,b=[sample(z,T) for z in (p,q)]
        check=max(np.max(abs(np.array(x['slow_and_fast_rates'])/y['slow_and_fast_rates']-1.))
                  for x,y in zip(a['finite_capacity'],b['finite_capacity'],strict=True))
        if check>1e-6:
            raise ArithmeticError('Physical collision modes not reproduced by independent quadrature')
        rows.append(dict(independent_max_relative_rate_difference=check,**b))
    return dict(rows=rows,scope='Near-equilibrium collision eigenvalues at fixed total energy with explicitly ideal visible heat capacities; not a Standard Model EOS or decoupling estimate.')


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
