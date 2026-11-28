"""Bounded discrete Bose/Pauli entropy-coordinate event updates."""
import json

import numpy as np
from scipy.optimize import brentq

from two_mass_kinetic_offset_deposition_reference import DT, event, shell


def exponents(n,g,bose):
    f=n/g
    if min(f)<=0. or np.max(f[~bose])>=1.:
        raise ValueError("Event must be evaluated inside the allowed occupations.")
    z=np.empty_like(f)
    z[bose]=np.log1p(1/f[bose])
    z[~bose]=np.log1p(-f[~bose])-np.log(f[~bose])
    return z


def entropy(n,g,bose):
    f=n/g
    S=np.empty_like(f)
    S[bose]=np.log1p(f[bose])+f[bose]*np.log1p(1/f[bose])
    u=f[~bose]
    S[~bose]=-u*np.log(u)-(1-u)*np.log1p(-u)
    return g@S


def batch(e,m,b,P,dirs):
    # Manufactured capacities for the participating charts, not a
    # density-of-states selection or an assertion of physical rates.
    p=[np.sqrt(event(e,m,b,P,mu)[0]) for mu in dirs]
    grids=[np.geomspace(min(z[i] for z in p)*.3,max(z[i] for z in p)*3,11)
           for i in range(2)]
    A=shell(grids[0],m)[0]
    B=shell(grids[1],b)[0]
    _,_,available=event(e,m,b,P,dirs[0])
    K=np.r_[available-DT(e),A,B]
    g=np.array([.9]+[.8+.03*i for i in range(len(K)-1)],dtype=DT)
    bose=np.array([False]+[True]*len(A)+[False]*len(B))
    vectors=[]
    for mu in dirs:
        _,offsets,_=event(e,m,b,P,mu)
        L=np.zeros(len(K),dtype=DT); L[0]=-1
        for a,c,k in zip((A,B),(1,1+len(A)),offsets,strict=True):
            i=int(np.searchsorted(a,k,side='right')-1)
            if not 0<=i<len(a)-1: raise ValueError("Outside the accepted energy stencil.")
            w=(k-a[i])/(a[i+1]-a[i])
            L[c+i:c+i+2]=[1-w,w]
        if abs(K@L-e)>3e-18*available or abs(L[1:1+len(A)].sum()-1)>3e-18 or abs(L[1+len(A):].sum()-1)>3e-18:
            raise ValueError("Event stoichiometry did not reproduce number or energy.")
        vectors.append(L)
    return K,g,bose,vectors


def rates(n,g,bose,v):
    a=exponents(n,g,bose)
    # The entropy conjugates are chosen at each receiving shell and
    # reconstruct the local effective Bose and Fermi occupations.
    x=-a[0]*v[0]
    y=(v*a)@bose
    z=(v*a)@(~bose)-v[0]*a[0]
    A=v@a
    log_fwd=-np.logaddexp(DT(0),x)-np.logaddexp(DT(0),-z)-np.log(-np.expm1(-y))
    r=np.sign(A)*np.exp(max(log_fwd,log_fwd-A))*(-np.expm1(-abs(A)))
    if r*A<0:
        raise ValueError("Discrete entropy production must be nonnegative.")
    return float(A),float(r)


def update(n,g,bose,v,step):
    A,S=rates(n,g,bose,v)
    if not S:
        return n.copy(),0.
    # Intersect all positive-population and Pauli capacity intervals.
    left,right=-np.inf,np.inf
    for x,c,u,s in zip(n,g,bose,v,strict=True):
        if s>0:
            left=max(left,-x/s)
            if not u: right=min(right,(c-x)/s)
        elif s<0:
            right=min(right,-x/s)
            if not u: left=max(left,(c-x)/s)
    if not left<0<right:
        raise ValueError("No feasible interval for the extent.")
    a,b=float(left)*(1-4e-13),float(right)*(1-4e-13)
    zeq=brentq(lambda t:rates(n+DT(t)*v,g,bose,v)[0],a,b,xtol=6e-16)
    if S*zeq<0:
        raise ValueError("Equilibrium lies opposite the signed event direction.")
    f=lambda t:t-step*rates(n+DT(t)*v,g,bose,v)[1]
    s,t=sorted([0.,zeq])
    if f(s)*f(t)>0.:
        # Equilibrium rounded to closer than the extent-resolution
        # threshold; never step outside the entropy-monotone segment.
        if abs(zeq)<2e-14:
            return n.copy(),0.
        raise ValueError("Extent was not bracketed by current and equilibrium states.")
    q=brentq(f,s,t,xtol=1e-15)
    n1=n+DT(q)*v
    sig=entropy(n1,g,bose)-entropy(n,g,bose)
    if sig<-2e-14 or abs(q-step*rates(n1,g,bose,v)[1])>2e-10:
        raise ValueError("Near-rest event update failed its entropy or implicit residual check.")
    return n1,float(q)


def run():
    groups=[]
    for title,e,m,b,P,T in (('ordinary',.093,.061,.122,.33,.07),
                            ('near-rest',2e-23,.061,.122,4e-12,4e-23)):
        K,g,bose,V=batch(e,m,b,P,(-.65,.21,.72))
        eta=np.where(bose,DT(-.4),DT(-.2))
        eta[0]=DT(-.4)+DT(-.2)-DT(e)/DT(T)
        a=K/DT(T)-eta
        occ=1/(np.exp(a)+np.where(bose,-1,1))
        eq=occ*g
        eaff=[rates(eq,g,bose,v)[0] for v in V]
        if max(abs(x) for x in eaff)>2e-17:
            raise ValueError(f"Thermal equilibrium affinity was not balanced: {eaff}.")
        experiments=[]
        for direction,scenario in (('forward','perturbed_parent'),('inverse','perturbed_parent'),
                                    ('forward','almost_empty_daughters'),('inverse','crowded_daughters')):
            n=eq.copy()
            filt=np.arange(len(K)); scalar=bose; fermion=(~bose)&(filt>0)
            n[0]=g[0]*DT(.7 if direction=='forward' else .002)
            if scenario=='almost_empty_daughters':
                n[0]=g[0]*DT(.999999)
                n[scalar|fermion]=DT(1e-10)*g[scalar|fermion]
            if scenario=='crowded_daughters':
                n[0]=g[0]*DT(1e-10)
                n[scalar]=DT(2.4)*g[scalar]
                n[fermion]=DT(.999999)*g[fermion]
            init=n.copy()
            hs=[]
            total=DT(0)
            for step in (.1,1.,8.,50.,10.):
                old=entropy(n,g,bose)
                for v in V:
                    n,d=update(n,g,bose,v,step)
                    total+=DT(d)
                k=(n-init)@K
                if (max(abs((n-init)[a].sum()-total) for a in (scalar,fermion))>3e-15 or
                      abs(n[0]-init[0]+total)>3e-15 or abs(k-DT(e)*total)>3e-15*max(K) or
                      entropy(n,g,bose)<old-2e-14):
                    raise ValueError("Global number, kinetic-rest, or entropy ledger failed.")
                hs.append(dict(dt=step,cumulative_event_extent=float(total),entropy_gain=float(entropy(n,g,bose)-entropy(init,g,bose)),
                               two_mass_excess_balance=float((k-DT(e)*total)/max(K)),
                               smallest_population_fraction=float(min(n/g)),max_fermi_population_fraction=float(max((n/g)[~bose]))))
            if np.sign(total)!=(1 if direction=='forward' else -1):
                raise ValueError("Expected both directions were not exercised.")
            experiments.append(dict(direction=direction,scenario=scenario,updates=hs))
        groups.append(dict(case=title,largest_thermal_equilibrium_affinity=max(abs(v) for v in eaff),experiments=experiments))
    return dict(scope="Manufactured entropy-coordinate stoichiometry, no physical source rate or actual cosmological dynamics.",
                cases=groups)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
