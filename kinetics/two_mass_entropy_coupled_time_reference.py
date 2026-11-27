"""Manufactured coupled-event time convergence on a near-rest entropy chart."""
import json

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from two_mass_entropy_extent_reference import DT, batch, entropy, rates, update


def interval(n,g,bose,v):
    left,right=-np.inf,np.inf
    for u,b,d,s in zip(n,g,bose,v,strict=True):
        if s>0:
            left=max(left,-u/s)
            if not d: right=min(right,(b-u)/s)
        elif s<0:
            right=min(right,-u/s)
            if not d: left=max(left,(b-u)/s)
    return float(left)*(1-4e-13),float(right)*(1-4e-13)


def midpoint(n,g,bose,v,step):
    A,S=rates(n,g,bose,v)
    if A==0.:
        return n.copy(),0.,0
    I=interval(n,g,bose,v)
    eq=brentq(lambda t:rates(n+DT(t)*v,g,bose,v)[0],*I,xtol=6e-16)
    a,b=sorted((0.,eq))
    f=lambda t:t-step*rates(n+DT(t/2)*v,g,bose,v)[1]
    if f(a)*f(b)>0.:
        m,q=update(n,g,bose,v,step)
        return m,q,1
    q=brentq(f,a,b,xtol=6e-16)
    m=n+DT(q)*v
    rates(m,g,bose,v)  # Finite interior and correct collision sign.
    if entropy(m,g,bose)<entropy(n,g,bose)-2e-14 or abs(f(q))>2e-10:
        raise ValueError("Midpoint step did not close its local entropy or implicit residual.")
    return m,q,0


def one(title,e,m,b,P,T):
    K,g,bose,V=batch(e,m,b,P,(-.65,.21,.72))
    eta=np.where(bose,DT(-.4),DT(-.2))
    eta[0]=DT(-.4)+DT(-.2)-DT(e)/DT(T)
    n0=g/(np.exp(K/DT(T)-eta)+np.where(bose,-1,1))
    n0[0]=g[0]*DT(.7)
    Z0=entropy(n0,g,bose)
    s=V[0]
    h=DT(2e-7)
    numeric=(entropy(n0+h*s,g,bose)-entropy(n0-h*s,g,bose))/(2*h)
    A=rates(n0,g,bose,s)[0]
    if abs(numeric-A)>2e-8:
        raise ValueError("Independent directional entropy check failed.")

    def rhs(_t,n):
        return np.array(sum((v*rates(n,g,bose,v)[1] for v in V),np.zeros_like(n)),float)
    span=.8
    estimates=[]
    for tol in (1e-10,3e-13):
        u=solve_ivp(rhs,(0.,span),np.array(n0,float),method='DOP853',rtol=tol,atol=tol*.004,max_step=.008)
        if not u.success: raise RuntimeError(u.message)
        estimates.append(u.y[:,-1])
    ref=estimates[-1]
    if np.linalg.norm(ref-estimates[0],np.inf)>5e-12:
        raise ValueError("Adaptive unsplit reference was not stable under tolerance refinement.")

    def ledgers(n,extent):
        change=n-n0
        rest=DT(e)*extent
        fermi=(~bose)&(np.arange(len(n))>0)
        if (max(abs(change[bose].sum()-extent),abs(change[fermi].sum()-extent),
                abs(change[0]+extent))>5e-15
                or abs(change@K-rest)>5e-15*max(K)
                or entropy(n,g,bose)<Z0-1e-14):
            raise ValueError("Global conservation or entropy ledger failed.")

    stats=[]
    for order in (1,2):
        mesh=[]
        for N in (8,16,32,64):
            n=n0.copy()
            extent=DT(0)
            fallback=0
            dx=span/N
            route=((0,1.),(1,1.),(2,1.)) if order==1 else ((0,.5),(1,.5),(2,1.),(1,.5),(0,.5))
            for _ in range(N):
                old=entropy(n,g,bose)
                for j,w in route:
                    if order==1: n,q=update(n,g,bose,V[j],dx*w)
                    else:
                        n,q,fb=midpoint(n,g,bose,V[j],dx*w)
                        fallback+=fb
                    extent+=DT(q)
                ledgers(n,extent)
                if entropy(n,g,bose)<old-3e-14:
                    raise ValueError("Stepwise entropy decreased.")
            mesh.append(dict(steps=N,max_population_error_vs_unsplit=float(np.max(abs(n-ref))),
                             midpoint_fallbacks=fallback,entropy_gain=float(entropy(n,g,bose)-Z0)))
        v=[x['max_population_error_vs_unsplit'] for x in mesh]
        observed=np.log2(v[-2]/v[-1])
        if abs(observed-order)>.17 or any(x['midpoint_fallbacks'] for x in mesh):
            raise ValueError("Manufactured convergence order was not reached in its intended smooth regime.")
        stats.append(dict(order=order,observed_last_rate=float(observed),mesh=mesh))

    # Separate genuinely stiff controls: midpoint is not guaranteed to be
    # bracketed, so fall back to the entropy-bounded backward event step.
    n=n0.copy()
    for j in np.flatnonzero(~bose):
        n[j]=g[j]*DT(.999999 if j==0 else 1e-10)
    n[bose]=g[bose]*DT(1e-10)
    z=entropy(n,g,bose)
    n,s,fb=midpoint(n,g,bose,V[0],100.)
    if not fb or entropy(n,g,bose)<=z or s<=0:
        raise ValueError("Stiff midpoint fallback did not exercise its intended guard.")
    return dict(title=title,nominal_span=span,directional_entropy_error=float(abs(numeric-A)),
                mesh=stats,stiff_fallback=dict(used=bool(fb),signed_extent=s,entropy_gain=float(entropy(n,g,bose)-z)))


if __name__=='__main__':
    print(json.dumps(dict(scope="Manufactured time-order and stiff fallback control for the three offset-shell events.",
        cases=[one('ordinary',.093,.061,.122,.33,.07),one('near-rest',2e-23,.061,.122,4e-12,4e-23)]),
                     indent=2,allow_nan=False))
