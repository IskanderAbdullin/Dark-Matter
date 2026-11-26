"""Optimistic vacuum one-pass chi2 depletion; arbitrary momenta, conditional injection."""
import json,math
from scipy.integrate import solve_ivp,quad
import numpy as np
from scipy.interpolate import PchipInterpolator
from functools import lru_cache
from kinetics.phi_bose_portal_low_reheat_gate import gamma,LEPTONS,M1,M2,MPL,GSTAR,S0,RHOC,MPHI
YUK=.001
K=YUK**4/(8*math.pi*M1*M1)
Icoef=2*math.pi**2/45*math.sqrt(GSTAR)*MPL/1.66

def trace_bound(s,t,m=M2,M=M1,mphi=MPHI):
    return (M*M-t)*s/2-(m*m+t-mphi*mphi+2*m*M)**2/2

@lru_cache(None)
def make_table():
    ts=np.geomspace(.0005,.1,160)
    vals=[]
    for T in ts:
        ss=2*math.pi**2/45*GSTAR*T**3; h=1.66*math.sqrt(GSTAR)*T*T/MPL
        vals.append(sum(gamma(T,m,ml) for m in (M1,M2) for ml in LEPTONS)/(ss*h))
    return PchipInterpolator(np.log(ts),np.log(vals))

def solve(TR,enhance=1.):
    table=make_table()
    def injection(z):
        T=TR*math.exp(-z)
        return math.exp(float(table(math.log(T)))) if T>=.0005 else 0.
    Y_in=quad(injection,0,5,epsabs=1e-35,epsrel=1e-8)[0]
    scale=K*enhance*Icoef*TR*Y_in
    def rhs(z,v):return [injection(z)/Y_in-scale*math.exp(-z)*v[0]**2]
    sol=solve_ivp(rhs,(0,5),[0.],method='Radau',rtol=2e-8,atol=1e-11)
    assert sol.success
    Y_opt=Y_in/(1/sol.y[0,-1]+scale*math.exp(-5))
    Y_bound=Y_in/(1+scale)
    assert Y_opt>=Y_bound*(1-1e-8)
    survive=Y_opt/Y_in
    target=2*(.01/.99*MPHI/M2)/(1+.01/.99*MPHI/M2)
    F_needed=math.sqrt(max(0.,(1/target-1)/(K*Icoef*TR*Y_in)))-1
    return dict(TR_MeV=TR*1000,Y_in=Y_in,omega_no_depletion=2*M2*Y_in*S0/RHOC,scale=scale,survival_instant_bound=Y_bound/Y_in,survival_distributed=survive,omega_lower_distributed=2*M2*Y_opt*S0/RHOC,pointwise_bose_cap_must_exceed_for_99pct_phi=max(0,F_needed))

if __name__=='__main__':
    for s in [4*M2*M2,.002,.01,.1,1.,100.,1e4]:
        d=math.sqrt(max(0,(s-4*M2*M2)*(s-4*MPHI*MPHI)))/2
        for c in [-1.,-.3,0.,.7,1.]:
            t=M2*M2+MPHI*MPHI-s/2+d*c
            tr=trace_bound(s,t)
            assert t<=1e-12 and tr>=-1e-12*s*s and tr<=(M1*M1-t)*s/2+1e-12
    rows=[solve(TR) for TR in (.005,.008,.01,.015,.02,.025,.03,.05,.1)]
    print(json.dumps(dict(K=K,Icoef=Icoef,rows=rows),indent=2))
