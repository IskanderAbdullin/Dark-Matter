"""Closure-independent pressure bracket at fixed comoving shell N and kinetic K."""
import json
import numpy as np
from scipy.special import roots_legendre, expit


def k(q,m):
    return q*q/(np.hypot(q,m)+m)


def p3(kin,m):
    # epsilon - m^2/epsilon; stable for k/m << 1.
    return kin*(1+m/(m+kin))


def pressure_bracket(edges,n,Q,m):
    K=Q/n
    kl,kh=k(edges[:-1],m),k(edges[1:],m)
    if np.any(K<kl) or np.any(K>kh):
        raise ValueError("shell kinetic moment outside support")
    a=(K-kl)/(kh-kl)
    lower=n@((1-a)*p3(kl,m)+a*p3(kh,m))
    upper=n@p3(K,m)
    return lower,upper


def model(N,a,distribution):
    # q is comoving; m=m_physical/Tref increases like a, and f(q) is frozen.
    m=(.008 if distribution=="bose_nonthermal" else .018)*a
    q=np.geomspace(2e-5,24.,N+1)
    z,w=roots_legendre(128)
    l,h=q[:-1,None],q[1:,None]
    xx=l+(h-l)*(z+1)/2
    weights=(h-l)*w/2*xx*xx
    f=(.02/np.expm1(xx)+100*np.exp(-.5*(np.log(xx/.003)/.22)**2)
       if distribution=="bose_nonthermal" else
       expit(.3-xx)*(.7+.2*np.sin(np.log(xx+.006))))
    w=weights*f
    kk=k(xx,m)
    num=w.sum(1)
    Q=(w*kk).sum(1)
    P=float(np.sum(w*p3(kk,m)))
    low,up=pressure_bracket(q,num,Q,m)
    assert low<=P<=up
    m0=m/a
    old=(w*k(xx,m0)).sum(1)/num
    low0,high0=k(q[:-1],m0),k(q[1:],m0)
    mix=(old-low0)/(high0-low0)
    def propagate(v):
        momentum2=v*(v+2*m0)
        return momentum2/(np.sqrt(momentum2+m*m)+m)
    kin_lower=num@propagate(old)
    kin_upper=num@((1-mix)*propagate(low0)+mix*propagate(high0))
    actual=Q.sum()
    assert kin_lower<=actual*(1+2e-14) and actual<=kin_upper*(1+2e-14)

    z=2e-5
    rp=float(np.sum(w*(k(xx,m*np.exp(z))+m*np.exp(z))))
    rm=float(np.sum(w*(k(xx,m*np.exp(-z))+m*np.exp(-z))))
    kp=float(np.sum(w*k(xx,m*np.exp(z))))
    km=float(np.sum(w*k(xx,m*np.exp(-z))))
    trace=float(np.sum(w*(m*m/(m+kk))))
    R=float(Q.sum()+m*num.sum())
    Qprime=float(-np.sum(w*m*kk/(m+kk)))
    trace_error=abs((rp-rm)/(2*z)-trace)/R
    kinetic_error=abs((kp-km)/(2*z)-Qprime)/Q.sum()
    assert trace_error<2e-9 and kinetic_error<2e-9
    assert abs((R-P)-trace)/R<2e-15
    return dict(distribution=distribution,N=N,a=a,m=m,actual_three_pressure=P,
                lower=low,upper=up,relative_scope=(up-low)/P,
                relative_trace_fd_error=trace_error,relative_kinetic_fd_error=kinetic_error,
                collisionless_final_kinetic=float(actual),
                collisionless_kinetic_lower=float(kin_lower),
                collisionless_kinetic_upper=float(kin_upper),
                collisionless_relative_scope=float((kin_upper-kin_lower)/actual))


def run():
    rows=[model(n,a,f) for f in ("bose_nonthermal","bounded_fermion")
          for n in (8,16,32,64,128) for a in (1.,30.,150.,1500.,10000.)]
    for f in ("bose_nonthermal","bounded_fermion"):
        coarse=max(r["relative_scope"] for r in rows if r["N"]==8 and r["distribution"]==f)
        fine=max(r["relative_scope"] for r in rows if r["N"]==128 and r["distribution"]==f)
        assert fine<coarse*.03
    return dict(rows=rows,status="manufactured collisionless FRW pressure-work control")
if __name__=="__main__":
    print(json.dumps(run(),indent=2,allow_nan=False))
