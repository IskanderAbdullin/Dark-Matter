"""Quasistatic isolated closure of chi1 <-> chi2 + phi conserving two counts.

The saturation marker means no positive-temperature normal Bose solution
with both conserved counts and the selected entropy. It is not the actual
coherent, condensate, or cosmological history of the interacting model.
"""
import json
from functools import lru_cache
import numpy as np
from scipy.special import roots_legendre,expit
from scipy.optimize import brentq
from three_species_decay_onsager_reference import M,m,mp

@lru_cache
def grid(n):
    z,w=roots_legendre(n)
    q=np.sqrt(65.)*(z+1)/2
    return q**2,np.sqrt(65.)*w*q**5/(2*np.pi**2)

def state(T,mi,g,eta,n=160,be=False):
    p,w=grid(n)
    e=np.hypot(p,mi/T)
    a=(p*p/(e+mi/T)-(eta-mi/T)) if be else e-eta
    if be:
        assert np.min(a)>0.
        f=1/np.expm1(a)
        entropy=a*f-np.log(-np.expm1(-a))
    else:
        f=expit(-a)
        entropy=a*f+np.logaddexp(0.,-a)
    return g*T**3*np.array([w@f,w@entropy,T*(w@(f*e))])

class Counts:
    def __init__(self,T0,eta2=0.,z0=None,nodes=160):
        self.n=nodes
        self.T0=T0
        self.z0=z0
        self.eta2=eta2
        z=-mp/T0 if z0 is None else z0
        f1,f2,b=self.states(T0,eta2,z)
        self.nf=f1[0]+f2[0]
        self.ns=f1[0]+b[0]
        self.S=f1[1]+f2[1]+b[1]

    def states(self,T,eta2,z):
        b=z+mp/T
        return state(T,M,4,eta2+b,self.n),state(T,m,4,eta2,self.n),state(T,mp,1,b,self.n,True)

    def atT(self,T,a):
        nf,ns=self.nf/a**3,self.ns/a**3
        def chi(z):
            b=z+mp/T
            s=lambda x: (state(T,M,4,x+b,self.n)[0]+state(T,m,4,x,self.n)[0])/nf-1
            ef=np.hypot((1.5*np.pi**2*nf)**(1/3),M)/T
            v=brentq(s,-100.,max(100.,ef+30),xtol=3e-13)
            return self.states(T,v,z),v
        (s1,s2,sb),eta2=chi(0.)
        deficit=ns-s1[0]-sb[0]
        if deficit>=0:
            # A bookkeeping complement at rest: it restores count but has
            # zero entropy. No dynamical inference uses this complement.
            return np.array([s1,s2,sb]),0.,eta2,deficit
        def target(z):
            zbd,_=chi(z)
            return (zbd[0][0]+zbd[2][0])/ns-1
        z=brentq(target,-100.,0.,xtol=4e-14)
        v,x=chi(z)
        return np.array(v),z,x,0.

    def solve(self,a):
        S=self.S/a**3
        def entropy(T):
            return self.atT(T,a)[0][:,1].sum()/S-1
        T=brentq(entropy,.09*self.T0/a,2.2*self.T0/a,xtol=2e-14)
        v,z,x,deficit=self.atT(T,a)
        n1,n2,nphi=v[:,0]
        r=max(abs((n1+n2)*a**3/self.nf-1),abs((n1+nphi+deficit)*a**3/self.ns-1),
              abs(v[:,1].sum()/S-1))
        assert r<4e-9,r
        return dict(a=a,T_dark=T,eta_2=x,fugacity_gap=z,
                    unsupported_normal_count_fraction=deficit/(deficit+nphi),
                    comoving_scalar_count=a**3*(nphi+deficit),normalized_ledger_error=r)

def test():
    out=[]
    for T,z in ((.6882487043,None),(.1,None),(.3,-.35),(.6882487043,-.35),(.6882487043,-2.)):
        S=Counts(T,z0=z)
        d=[S.solve(a) for a in (1.,1.03,1.15,1.5,2.,3.,5.,10.,20.,40.)]
        r=Counts(T,z0=z,nodes=232)
        cross=max(abs(r.solve(x["a"])["T_dark"]/x["T_dark"]-1) for x in d[::3])
        assert cross<2e-5,cross
        assert abs(d[0]["T_dark"]/T-1)<1e-8
        if z is not None: assert abs(d[0]["fugacity_gap"]-z)<5e-8
        hit=[o for o in d if o["unsupported_normal_count_fraction"]>1e-8]
        # Locate saturation in log a using the fixed-entropy saturated
        # normal count as a signed root. Report when a reference state
        # lies tangent to the threshold, not as a dynamical transition.
        def threshold(la):
            a=np.exp(la)
            nf=S.nf/a**3
            sm=S.S/a**3
            def trial(t):
                f=lambda e:sum(v[0] for v in S.states(t,e,0.)[:2])/nf-1
                e=brentq(f,-100.,max(100.,np.hypot((1.5*np.pi**2*nf)**(1/3),M)/t+30))
                return S.states(t,e,0.)
            t=brentq(lambda t:sum(o[1] for o in trial(t))/sm-1,.09*T/a,2.2*T/a)
            f1,f2,b=trial(t)
            return (f1[0]+b[0])*a**3/S.ns-1
        onset=None
        brackets=np.linspace(0.,np.log(40.),90)
        for a,b in zip(brackets[:-1],brackets[1:]):
            u,v=threshold(a),threshold(b)
            if u>0 and v<=0:
                onset=float(np.exp(brentq(threshold,a,b)))
                break
            if a==0 and abs(u)<1e-9 and v<0:
                onset=1.
                break
        if hit:assert onset is not None
        out.append(dict(initial_T=T,initial_gap=z,saturation_a=onset,values=d,
                        max_refined_T_relative_difference=cross))
    return dict(rows=out,scope="Quasistatic isolated decay-only adiabatic reference, no other collisions. A deficit of normal-state support is not a physical condensate prediction.")

if __name__=="__main__":
    print(json.dumps(test(),indent=2,allow_nan=False))
