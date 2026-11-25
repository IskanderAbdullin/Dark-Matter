"""Conservative finite-age nonreturning exterior population from interior births."""
import argparse,json,math
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import brentq
from lowered_antonov_hessian import Halo
from lowered_bose_exterior_source import density,exterior_min,R,V,GYR,TIME
from lowered_bose_collision_escape import rate

def run(H,n=24,nr=80,ne=12):
    e0=H.Qt/H.rt;emax=e0/2;x,w=leggauss(nr);x=(x+1)/2;w/=2;r=H.rt*x*x;ww=w*2*H.rt*x*r*r;ps=H.psi(r)
    rho=V**2/(4*math.pi*6.6743e-11*R**2);pref=.5*rho*.125*math.sqrt(2)*V*GYR/(H.Ic**2*H.Qt)
    S=lambda eps:pref*(ww@np.array([density(H,p,e0-eps,n) for p in ps]))
    cum=lambda e:pref*(ww@np.array([rate(H,p,e,n)[1] for p in ps]))
    z,wz=leggauss(ne);z=(z+1)/2;wz/=2;rows=[]
    for T in [5.,10.,13.8]:
        eps_ext=emax if exterior_min(H.Qt,H.rt,emax)>=T else brentq(lambda e:exterior_min(H.Qt,H.rt,e)-T,emax*1e-8,emax)
        eps_delay=max(0.,(e0-(2*H.rt*TIME/T)**2)/2)
        lim=min(eps_ext,eps_delay)
        ep=lim*z;s=np.array([S(e) for e in ep]);delay=2*H.rt*TIME/np.sqrt(e0-2*ep)
        lower=lim*(wz@(s*(T-delay)));born=T*(cum(e0-eps_ext)-cum(e0))
        rows.append(dict(age_Gyr=T,eps_exterior=eps_ext,eps_delay=eps_delay,eps_used=lim,lower=lower,weak_born_upper=born,lower_over_born=lower/born))
    return dict(model=H.model,n=n,nr=nr,ne=ne,rows=rows)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--model',default='linear');p.add_argument('--n',type=int,default=24);p.add_argument('--nr',type=int,default=80);p.add_argument('--ne',type=int,default=12);a=p.parse_args()
    print(json.dumps(run(Halo(model=a.model),a.n,a.nr,a.ne)))
