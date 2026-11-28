"""Separate heat and number nullclines of the local charged-lepton FD portal."""
import json
from functools import lru_cache

import numpy as np
from scipy.optimize import brentq

from boost_axis_elastic_heat_reference import integrate, leptons
from positive_fd_portal_quadrature_reference import PositivePortal
from finite_width_vector_pauli_rate_reference import integrator


def roots(T, ratio, fine):
    dark = ratio*T
    p = PositivePortal((13, 36, 26) if fine else (9, 24, 16), tail=65. if fine else 55.)
    eta_tol = 5e-13
    @lru_cache(maxsize=None)
    def timelike(eta):
        return p.flow(T,dark,eta)

    @lru_cache(maxsize=None)
    def elastic(eta):
        q = -sum(integrate(T,dark,eta,m,n=22 if fine else 16,cut=35. if fine else 29.)
                 for m in leptons)
        if ratio!=1. and q*(ratio-1.)<=0:
            raise ValueError("Elastic heat does not flow to the cooler bath.")
        return q

    def find(f):
        rt = brentq(f,-15.,15.,xtol=eta_tol)
        h = 3e-6
        lo, hi = f(rt-h), f(rt+h)
        der = (hi-lo)/(2*h)
        if not lo<0<hi or abs(f(rt))>2e-10*der:
            raise ValueError("Root or positive local slope failed.")
        return rt, float(der)

    en, dn = find(lambda x:timelike(x)[0])
    ee, de = find(lambda x:timelike(x)[1])
    ec, dc = find(lambda x:timelike(x)[1]+elastic(x))
    if ratio!=1.:
        if not (ratio-1.)*(en-ee)>0 or not (ratio-1.)*(ee-ec)>0:
            raise ValueError("Expected heat-number separation or elastic shift did not close.")

    def response(eta):
        nd, qd, sd = timelike(eta)
        q = elastic(eta)
        beta = (dark-T)/(dark*T)
        if abs(sd-(2*eta*nd+beta*qd))>2e-10*max(abs(2*eta*nd),abs(beta*qd),1e-290):
            raise ValueError("Event entropy did not match the moment ledger.")
        if sd<0 or beta*q<0:
            raise ValueError("Unexpected entropy production sign.")
        return dict(number=float(nd), timelike_heat=float(qd), elastic_heat=q,
                    timelike_entropy=float(sd), combined_entropy=float(sd+beta*q))

    return dict(number_eta=en,heat_eta=ee,combined_heat_eta=ec,number_local_slope=dn,
                heat_local_slope=de,combined_heat_local_slope=dc,
                responses={k:response(e) for k,e in
                           (('number_null',en),('timelike_heat_null',ee),('combined_heat_null',ec))})


def exact_common_derivatives(p,T):
    _,_,_,W,M1,M2,_ = p.moments(T,T,0.)
    kel=sum(integrate(T,T,0.,m,n=22,cut=35.,linear=True) for m in leptons)
    return np.array([-M1/(2*T*W),-M2/(2*T*M1),-(M2/T+kel)/(2*M1)])


def run():
    rows=[]
    adaptive=[]
    p=PositivePortal((13,36,26),tail=65.)
    keys=('number_eta','heat_eta','combined_heat_eta')
    for T in (.015,.05,.1,.3,.6882487043):
        for ratio in (.8,.95,1.05,1.2):
            low=roots(T,ratio,False)
            f=roots(T,ratio,True)
            diff=np.array([f[k]-low[k] for k in keys])
            if max(abs(diff))>2e-5:
                raise ValueError(f"Eta threshold precision gate failed: {T}, {ratio}, {diff}.")
            rows.append(dict(T_SM=T, ratio=ratio, solved=f, eta_order_refinement=diff.tolist()))
            if T in (.05,.6882487043) and ratio in (.8,1.2):
                err=[]
                for i,k in enumerate(keys[:2]):
                    eta=f[k]
                    ad=integrator(T*ratio,T,eta)[3,i]
                    response=f['number_local_slope' if i==0 else 'heat_local_slope']
                    err.append(float(ad/response))
                if max(map(abs,err))>5e-6:
                    raise ValueError(f"Independent adaptive threshold check failed: {T}, {ratio}, {err}.")
                adaptive.append(dict(T_SM=T,ratio=ratio,adaptive_residual_in_eta_units=err))
        h=3e-4
        a,b=(roots(T,1.+i*h,True) for i in (-1,1))
        fd=np.array([(b[k]-a[k])/(2*h) for k in keys])
        an=exact_common_derivatives(p,T)
        if np.max(abs(an-fd)/np.maximum(abs(an),1.))>3e-5:
            raise ValueError(f"Common-temperature moment limit disagrees: {T}, {an}, {fd}.")
        rows.append(dict(T_SM=T, ratio=1., exact_common_ratio_derivatives=an.tolist(),
                         finite_difference_common_ratio_derivatives=fd.tolist()))
    return dict(scope="Local symmetric pair FD bath, charged leptons and finite-width effective vacuum portal only.",
                roots_order=list(keys),rows=rows,adaptive_quadrature_cross_checks=adaptive)


if __name__ == '__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
