"""Recover an exact convex test-force potential without assuming physical integrability."""
import json
import numpy as np

from fd_portal_differential_reciprocity_reference import integrate,PositivePortal


def potential(z,g,a,b):
    h=(a+b*z)/2
    y=8*np.sum(g*np.sinh(h/2)**2)
    t=2*g*np.sinh(h)
    v=g*np.cosh(h)
    F=np.array([np.sum(t),np.sum(z*t)])
    H=np.array([[np.sum(v),np.sum(z*v)],[np.sum(z*v),np.sum(z*z*v)]])
    return y,F,H


def verify(p,T,a,b):
    f,D,skew,_,z,g,E=integrate(p,T,a,b,return_frozen_activity=True)
    y,F,H=potential(z,g,a,b)
    n=max(np.max(abs(D)),np.max(abs(H)))
    if y<0. or min(np.linalg.eigvalsh(H))<=0. or np.min(g)<0.:
        raise ArithmeticError('Expected a convex even frozen activity potential')
    err=max(np.max(abs(F-f[:2])),np.max(abs(D-H-E)))/n
    if err>3e-14 or abs(E[0,1]-E[1,0]-skew)>3e-14*n:
        raise ArithmeticError('Frozen activity potential fails its actual state derivative split')
    t=8e-5
    def values(i,h):
        x,u=a+h*(i==0),b+h*(i==1)
        y,F,_=potential(z,g,x,u)
        return np.r_[y,F]
    W=np.array([(values(i,-2*t)-8*values(i,-t)+8*values(i,t)-values(i,2*t))/(12*t)
                for i in (0,1)]).T
    err2=np.max(abs(W-np.vstack((F,H))))/n
    if err2>3e-8:
        raise ArithmeticError('Frozen test-force derivatives failed their independent numerical probe')
    return dict(actual_flows=F.tolist(),convex_test_force_hessian=H.tolist(),physical_differential_jacobian=D.tolist(),
                separately_derived_activity_response=E.tolist(),potential_value=y,
                max_exact_split_error_scaled=float(err),max_frozen_finite_difference_error=float(err2),
                max_activity_response_relative_to_physical=float(np.max(abs(E))/n))


def run():
    p,q=(PositivePortal(o) for o in ((13,36,26),(17,44,32)))
    out=[]
    for T,a,b in ((.05,0.,0.),(.05,.6,.16),(.05,-.6,-.16),(.3,-1.5,.16),(.3,1.5,-.16),
                  (.64135,-.03418787777103929,.002),(.6882487043,0.,0.),(.6882487043,1.8,-.16)):
        u,v=[verify(x,T,a,b) for x in (p,q)]
        h,k=[np.asarray(x['convex_test_force_hessian']) for x in (u,v)]
        if np.max(abs(h/k-1.))>2e-6:
            raise ArithmeticError('Two orders did not reproduce the same convex activity Hessian')
        out.append(dict(T_SM=T,forces=[a,b],**v))
    return dict(rows=out,scope='The activity is held fixed only while varying a mathematical test force; the physical bath derivative changes it. No global potential of the physical force-dependent flow or cosmological stability claim.')


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
