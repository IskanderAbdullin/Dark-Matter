"""Regular BE two-moment boundary on a near-rest lowest momentum shell."""
import json
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq


def capacity(Q,m,beta):
    high=Q*Q/(np.hypot(Q,m)+m)
    def fn(x,order):
        q=Q*x
        s=q*q/(np.hypot(q,m)+m)/high if x else 0.
        if x==0:
            return (6*m*high/(Q*Q*beta) if order==0 else 0.)
        with np.errstate(over="ignore"):
            return 3*x*x*(s if order else 1.)/np.expm1(beta*s)
    cuts=[1/np.sqrt(beta)] if beta>1 else []
    return np.array([quad(lambda z:fn(z,i),0.,1.,points=cuts,
                          epsabs=0.,epsrel=1e-10)[0] for i in (0,1)])


def independent(Q,m,beta):
    K=Q*Q/(np.hypot(m,Q)+m)
    # Changes variables from q to sqrt((E-m)/K), independently of the
    # normalized-q integrand in capacity().
    def fn(x,i):
        eps=K*x*x
        q=np.sqrt(eps*(eps+2*m))
        jac=6*q*(eps+m)*K*x/Q**3
        with np.errstate(over="ignore"):
            return jac*(x*x if i else 1)/np.expm1(beta*x*x)
    return np.array([quad(lambda x:fn(x,k),0.,1.,
                          points=[1/np.sqrt(beta)] if beta>1 else [],
                          epsrel=1e-10,epsabs=0.)[0] for k in (0,1)])


def edge(Q,m,reln):
    hi=1.
    while capacity(Q,m,hi)[0]>reln:
        hi*=2.
    lo=hi
    while capacity(Q,m,lo)[0]<reln:
        lo*=.5
    beta=brentq(lambda z:capacity(Q,m,np.exp(z))[0]/reln-1.,
                np.log(lo),np.log(hi),xtol=1e-12)
    beta=np.exp(beta)
    C=capacity(Q,m,beta)
    ref=independent(Q,m,beta)
    err=float(np.max(abs(C/ref-1)))
    assert abs(C[0]/reln-1.)<1e-9 and err<1e-9
    return dict(Q=Q,mass=m,mean_fill=reln,slope_boundary=beta,
                min_thermal_mean_kinetic_fraction=float(C[1]/C[0]),
                relative_independent_moment_error=err)


def run():
    rows=[edge(Q,m,x) for Q,m in [(1e-7,.1),(.002,.0004)] for x in (.001,.01,.1,1.,10.,100.,1.e5)]
    near=rows[6]
    assert abs(near["min_thermal_mean_kinetic_fraction"]-1/3)<5e-6
    return dict(rows=rows,near_rest_large_occupation_limit=1/3,
                status="Regular local thermal shape feasibility only, not BEC evidence.")

if __name__ == "__main__":
    print(json.dumps(run(),indent=2,allow_nan=False))
