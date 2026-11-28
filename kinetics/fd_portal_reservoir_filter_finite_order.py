"""Check second-order agreement of filtered FD steps at finite displacements."""
import json

import numpy as np
from scipy.integrate import solve_ivp

from fd_portal_differential_reciprocity_reference import PositivePortal, integrate
from fd_portal_implicit_reservoir_reference import dark
from fd_portal_discrete_gradient_reservoir_reference import trajectory
from fd_portal_reservoir_exponential_filter_reference import filtered


def run():
    rows=[]
    p,q=[PositivePortal(o) for o in ((13,36,26),(17,44,32))]
    for T,start in ((.05,[.6,.16]),(.3,[-1.5,.16]),(.6882487043,[1.8,-.16])):
        a=np.asarray(start)
        y0,_,_=dark(T,a,168)
        f=integrate(q,T,*a)[0]
        scale=np.array([2.,T])
        tau=.12/np.max(abs(f[:2]*scale/y0))

        def rhs(t,u):
            _,_,A=dark(T,u,168)
            return np.linalg.solve(A,-tau*scale*integrate(q,T,*u)[0][:2])

        r=solve_ivp(rhs,(0.,1.),a,method='DOP853',rtol=3e-12,atol=4e-14)
        if not r.success:
            raise ArithmeticError('Missing clean unfiltered reference on this interval')
        target=dark(T,r.y[:,-1],168)[0]

        def steps(chart,nodes,n):
            z=a.copy()
            defect=[]
            rate=[]
            for _ in range(n):
                z,h,j=filtered(chart,T,z,tau/n,nodes)
                defect.append(h)
                rate.append(j)
            return z,max(defect),[min(rate),max(rate)]

        errors=[]
        for n in (2,4,8,16):
            u,e,z=steps(q,168,n)
            unfiltered,v,_,_=trajectory(q,T,a,tau,n,168)
            err=np.max(abs(dark(T,u,168)[0]-target)/y0)
            raw=np.max(abs(v-target)/y0)
            if not 0.<z[0]<=z[1]<=1.+2e-12:
                raise ArithmeticError('Filtered rate escaped its positive contraction interval')
            errors.append(dict(steps=n,filtered_force=u.tolist(),unfiltered_force=unfiltered.tolist(),
                               filtered_moment_error=err,unfiltered_moment_error=raw,maximum_entropy_identity_error=e,
                               step_entropy_rate_ratio_range=z))
        err=np.asarray([x['filtered_moment_error'] for x in errors])
        ratio=err[:-1]/err[1:]
        w,_,_=steps(p,112,8)
        mesh=np.max(abs(w-errors[2]['filtered_force']))
        if (np.min(ratio)<3.2 or np.max(ratio)>5.0 or mesh>2e-7):
            raise ArithmeticError('Finite-displacement second-order or paired-mesh check failed')
        rows.append(dict(T_SM=T,initial_forces=start,local_time=tau,steps=errors,filtered_error_halving_ratios=ratio.tolist(),
                         independent_mesh_difference_in_eight_step_forces=mesh))
    return dict(rows=rows,scope='Fixed visible reservoir finite-displacement temporal refinement. The exact target is the unfiltered moment ODE, not a cosmological trajectory.')


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
