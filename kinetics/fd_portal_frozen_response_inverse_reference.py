"""Conditioned inverse of the frozen-activity FD source; no evolving bath inversion."""
import json
import numpy as np
from fd_portal_differential_reciprocity_reference import integrate, PositivePortal
from fd_portal_frozen_activity_potential_reference import potential


def solve(z,g,target,initial):
    u=np.asarray(initial,dtype=float)
    h=potential(z,g,*u)[2]
    scale=np.sqrt(np.diag(h))
    if min(scale)<=0:
        raise ArithmeticError('Missing fixed activity support')
    K=max(np.linalg.norm(target/scale),np.finfo(float).tiny)
    for it in range(40):
        val,F,H=potential(z,g,*u)
        res=(F-target)/scale
        r=np.linalg.norm(res)/K
        if r<2e-13:
            return u,it,r
        du=np.linalg.solve(H/scale[:,None]/scale,res)/scale
        for j in range(32):
            v=u-du*2.**(-j)
            f=potential(z,g,*v)[1]
            if np.linalg.norm((f-target)/scale)/K < r:
                u=v
                break
        else:
            raise ArithmeticError(('Damped Newton did not improve actual residual',r))
    raise ArithmeticError('Frozen force inverse did not converge')


def case(p,T,a,b):
    f,D,skew,_,z,g,_=integrate(p,T,a,b,return_frozen_activity=True)
    u=np.array([a,b])
    val,F,H=potential(z,g,*u)
    target=u+np.array([.12,-.008])
    hot,goal,_=potential(z,g,*target)
    x,it,res=solve(z,g,goal,(0.,0.))
    if np.max(abs(x-target))>3e-9:
        raise ArithmeticError('Cannot recover a prescribed nearby force')
    y,_,_=solve(z,g,F,(0.,0.))
    if np.max(abs(y-u))>3e-9:
        raise ArithmeticError('Cannot recover the actual fixed-state force')
    # At fixed F, the heat response to beta is the Schur complement
    # of the frozen Hessian; the inverse probes the same held activity.
    S=H[1,1]-H[0,1]**2/H[0,0]
    if S<=0 or abs(u@F-f[2])>3e-14*max(np.max(abs(f)),np.max(abs(H))):
        raise ArithmeticError('Fixed entropy pairing or Schur complement failed')
    d=.00035*S
    s,_i,_r=solve(z,g,F+np.array([0.,d]),u)
    t,_i,_r=solve(z,g,F-np.array([0.,d]),u)
    measured=(s[1]-t[1])/(2*.00035)
    if abs(measured-1)>3e-6:
        raise ArithmeticError('Independent heat-only inverse disagrees with the frozen Schur response')
    # The convex conjugate evaluated at the exact inverse goal.
    dual=x@goal-potential(z,g,*x)[0]
    gap_a=hot+dual-target@goal
    trial_gap=val+dual-u@goal
    if abs(gap_a)>3e-13*max(abs(hot),abs(dual)) or trial_gap<=0.:
        raise ArithmeticError('Fenchel equality or strict off-target convex gap failed')
    return dict(visible_temperature=T,original_forces=u.tolist(),oracle_test_force=target.tolist(),
                recovered_test_force=x.tolist(),iterations=it,scaled_residual=res,
                frozen_thermal_schur_over_heat_diagonal=S/H[1,1],
                independent_heat_only_inverse_ratio=measured,off_target_fenchel_gap_over_target_activity=trial_gap / hot,
                original_force_recovery_max_error=float(np.max(abs(y-u))))


def run():
    p,q=[PositivePortal(v) for v in ((13,36,26),(17,44,32))]
    rows=[]
    for T,a,b in ((.05,.6,.16),(.05,-.6,-.16),(.3,-1.5,.16),(.64135,-.03418787777103929,.002),
                  (.6882487043,1.8,-.16)):
        x,y=[case(v,T,a,b) for v in (p,q)]
        err=abs(x['frozen_thermal_schur_over_heat_diagonal']/y['frozen_thermal_schur_over_heat_diagonal']-1.)
        if err>2e-6:
            raise ArithmeticError('Two positive quadratures did not reproduce the conditioned inverse')
        rows.append(dict(independent_order_defect=err,**y))
    return dict(rows=rows,scope='Inverse of the fixed-state convex test-force map. It is not an inversion of physical forces with changing activity, nor a bound on bath evolution.')


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
