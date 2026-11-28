"""Invert the finite local bath-dependent rate using the true FD tangent."""
import json
import numpy as np

from fd_portal_differential_reciprocity_reference import PositivePortal,integrate
from fd_portal_frozen_response_inverse_reference import solve


def inverse(p,T,initial,target,method):
    u=np.asarray(initial,dtype=float)
    def value(x):
        f,D,*_=integrate(p,T,*x)
        return f[:2],D
    F,D=value(u)
    z=np.sqrt(np.diag(D))
    norm=lambda v: np.linalg.norm(v/z)
    scale=max(norm(target),np.finfo(float).tiny)
    for it in range(32):
        F,D=value(u)
        err=norm(F-target)/scale
        if err<1.5e-12:
            return u,it,err
        if method=='finite_difference':
            a=8e-6
            C=np.empty((2,2))
            for k in range(2):
                v=np.eye(2)[k]*a
                C[:,k]=(p.moments(T,T/(1-u[1]-v[1]),(u[0]+v[0])/2)[:2]-
                        p.moments(T,T/(1-u[1]+v[1]),(u[0]-v[0])/2)[:2])/(2*a*np.array([1.,T]))
            D=C
        du=np.linalg.solve(D/z[:,None]/z,(F-target)/z)/z
        for n in range(24):
            x=u-du*2.**(-n)
            if x[1]<.6 and x[1]>-.6 and abs(x[0])<8. and norm(value(x)[0]-target)/scale<err:
                u=x
                break
        else:
            raise ArithmeticError(('local physical inverse not resolved',T,err))
    raise ArithmeticError('Not a verified local inverse')


def run():
    p,q=[PositivePortal(o) for o in ((13,36,26),(17,44,32))]
    rows=[]
    for T,a,b in ((.05,.6,.16),(.05,-.6,-.16),(.3,-1.5,.16),
                  (.64135,-.034187877771,.002),(.6882487043,1.8,-.16)):
        initial=np.array([a,b])
        want=initial+np.array([.03,-.0018])
        goal=integrate(q,T,*want)[0][:2]
        final=[]
        for y in (p,q):
            target=integrate(y,T,*want)[0][:2]
            x,k,e=inverse(y,T,initial,target,'analytic')
            if np.max(abs(x-want))>3e-8:
                raise ArithmeticError('Could not reproduce oracle held-model forces')
            final.append((x,k,e))
        x,n,v=inverse(q,T,initial,goal,'finite_difference')
        if np.max(abs(x-final[1][0]))>3e-8:
            raise ArithmeticError('Independent finite-difference solution is inconsistent')
        f,D,skew,_,z,g,extra=integrate(q,T,a,b,return_frozen_activity=True)
        w,i,h=solve(z,g,goal,initial)
        got=integrate(q,T,*w)[0][:2]
        rel=np.linalg.norm((got-goal)/np.sqrt(np.diag(D)))/np.linalg.norm(goal/np.sqrt(np.diag(D)))
        if rel<=1e-7:
            raise ArithmeticError('Toy control did not resolve a distinct physical and frozen response')
        rows.append(dict(T_SM=T,initial_forces=initial.tolist(),target_forces=want.tolist(),
                        physical_inverse_forces=final[-1][0].tolist(),
                        paired_analytic_iterations=[d[1] for d in final],
                        paired_analytic_residuals=[d[2] for d in final],
                        independent_finite_difference_iterations=n,independent_finite_difference_residual=v,
                        frozen_activity_surrogate_forces=w.tolist(),frozen_surrogate_actual_physical_residual=rel))
    return dict(rows=rows,scope='Local numerical inverse at fixed visible temperature and one timelike two-bath kernel. Starts from nearby specified states; no global inverse or cosmological evolution/stability claim.')


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
