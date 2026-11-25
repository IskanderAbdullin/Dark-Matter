"""First exterior excursion response for a prescribed constant injection at r_t.
Frozen spherical potential, no rescattering outside, not an evolved halo.
"""
import argparse,json,math
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import brentq
from lowered_antonov_hessian import Halo
from lowered_bose_exterior_source import density,exterior_min,period,R,V,GYR,TIME
from lowered_bose_collision_escape import rate

def run(H,n=24,nr=64,ne=12):
    e0=H.Qt/H.rt;emax=e0/2;x,w=leggauss(nr);x=(x+1)/2;w/=2;r=H.rt*x*x;ww=w*2*H.rt*x*r*r;ps=H.psi(r)
    rho=V**2/(4*math.pi*6.6743e-11*R**2);pref=.5*rho*.125*math.sqrt(2)*V*GYR/(H.Ic**2*H.Qt)
    cum=lambda e:pref*(ww@np.array([rate(H,p,e,n)[1] for p in ps]))
    S=lambda eps:pref*(ww@np.array([density(H,p,e0-eps,n) for p in ps]))
    z,wz=leggauss(ne);z=(z+1)/2;wz/=2
    low=lambda e:exterior_min(H.Qt,H.rt,e);high=lambda e:period(H.Qt,e)
    s0=S(0);total=cum(e0-emax)-cum(e0);C=TIME*2*math.pi*H.Qt/2**1.5;rows=[]
    for T in [1.,5.,13.8,100.,1000.]:
        root=lambda fun:emax if fun(emax)>=T else brentq(lambda e:fun(e)-T,emax*1e-12,emax)
        edges=sorted(set([0.,root(low),root(high),emax]));N=np.zeros(3)
        for a,b in zip(edges[:-1],edges[1:]):
            ep=a+(b-a)*z;s=np.array([S(e) for e in ep])
            N+=(b-a)*np.array([wz@s,wz@(s*np.minimum(T,[low(e) for e in ep])),wz@(s*np.minimum(T,[high(e) for e in ep]))])
        rows.append(dict(age_Gyr=T,entered_fraction=T*N[0],outside_lower=N[1],outside_upper=N[2],returned_lower=max(0.,T*N[0]-N[2]),returned_upper=max(0.,T*N[0]-N[1]),source_closure=abs(N[0]/total-1),asymptotic=3*s0*C**(2/3)*T**(1/3)))
    return dict(model=H.model,n=n,nr=nr,ne=ne,emax=emax,source_per_Gyr=total,source_at_zero=s0,rows=rows)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--model',default='linear');p.add_argument('--n',type=int,default=24);p.add_argument('--nr',type=int,default=64);p.add_argument('--ne',type=int,default=12);a=p.parse_args()
    print(json.dumps(run(Halo(model=a.model),a.n,a.nr,a.ne)))
