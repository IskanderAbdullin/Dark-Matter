"""Check normalized local FD inverse sensitivity along weak and strong directions."""
import json
import numpy as np

from fd_portal_differential_reciprocity_reference import PositivePortal,integrate
from fd_portal_physical_local_inverse_reference import inverse


def run(p,T,u):
    y,D,_,_=integrate(p,T,*u)
    f=y[:2]
    s=np.sqrt(np.diag(D))
    M=D/s[:,None]/s
    U,V,W=np.linalg.svd(M)
    if min(V)<=0 or not np.all(np.isfinite(V)):
        raise ArithmeticError('Degenerate local physical tangent')
    res=[]
    for i in (0,1):
        # A force perturbation with bounded amplitude defines an exact
        # rate target. Use its right-singular direction in scaled forces.
        dt=2e-6/np.max(abs(W[i]/s))
        v=dt*(W[i]/s)
        go=[integrate(p,T,*x)[0][:2] for x in (u+v,u-v)]
        # Independently invert the two true targets and also evaluate the
        # finite antisymmetric rate difference along the prescribed direction.
        sols=[inverse(p,T,u,x,'analytic')[0] for x in go]
        exact_rate_direction=(go[0]-go[1])/(2*dt*s)
        k=np.linalg.norm(exact_rate_direction-U[:,i]*V[i])/V[i]
        re=max(np.max(abs(sols[0]-(u+v))),np.max(abs(sols[1]-(u-v))))
        if k>3e-5 or re>1e-8:
            raise ArithmeticError('Finite shifts or independently recovered physical forces disagree')
        res.append(dict(direction='strong' if i==0 else 'weak',expected_normalized_gain=V[i],
                        measured_normalized_gain=float(np.linalg.norm(exact_rate_direction)),
                        finite_derivative_direction_error=k,roundtrip_force_max_absolute_error=re))
    return dict(T_SM=T,forces=u.tolist(),scaled_tangent_condition_number=V[0]/V[1],
                exact_physical_roundtrips=res)


def test_pairs():
    p,q=[PositivePortal(x) for x in ((13,36,26),(17,44,32))]
    rows=[]
    for T,a,b in ((.05,.6,.16),(.3,-1.5,.16),(.64135,-.034187877771,.002),(.6882487043,1.8,-.16)):
        u=np.array([a,b]); x,y=[run(k,T,u) for k in (p,q)]
        err=abs(x['scaled_tangent_condition_number']/y['scaled_tangent_condition_number']-1.)
        if err>2e-6:
            raise ArithmeticError('Independent quadrature does not reproduce local conditioning')
        rows.append(dict(independent_two_order_defect=err,**y))
    return dict(rows=rows,scope='Diagonal-normalized local sensitivity of the true two-bath FD number/heat response. Numerical direction and finite oracle controls only; no assumed observational noise or cosmological error bar.')


if __name__=='__main__':
    print(json.dumps(test_pairs(),indent=2,allow_nan=False))
