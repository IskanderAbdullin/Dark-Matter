"""Initial weak-binding source and exact frozen Kepler exterior dwell-time bounds."""
import argparse,json,math
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import brentq
from lowered_antonov_hessian import Halo
from lowered_bose_collision_escape import rate

R=.6*3.085677581491367e19; V=44100.; GYR=1e9*365.25*86400
TIME=R/V/GYR

def density(H,psi,excess,n):
    """Minus derivative of cumulative stimulated birth integral wrt threshold."""
    if psi<=excess or psi<=0:return 0.
    x,w=leggauss(n);x=(x+1)/2;w=w/2
    a=excess/psi+(1-excess/psi)*x[:,None,None]
    b0=1+excess/psi-a;b=b0+(1-b0)*x[None,:,None]
    q1=psi*a;q2=psi*b;S=q1+q2;P=psi+excess;A=S/2
    cm=np.sqrt(np.maximum(0,P*(S-P)/(q1*q2)));c=cm*(2*x[None,None,:]-1)
    B=np.sqrt(np.maximum(0,A*A-q1*q2*c*c));lo=psi+P-S
    qfac=(psi**3/4)*(1-excess/psi)*(1-b0)*np.sqrt(a*b)*H.F(psi-q1)*H.F(psi-q2)
    rel=np.sqrt(np.maximum(0,S-2*np.sqrt(q1*q2)*c))
    return float(np.sum(w[:,None,None]*w[None,:,None]*w[None,None,:]*qfac*cm*rel*(1+H.F(lo))/B))

def period(Q,eps):return TIME*2*math.pi*Q/(2*eps)**1.5

def exterior_min(Q,rt,eps):
    # Bound only for a>=rt, where x=rt/a<=1.
    x=2*eps*rt/Q
    if not 0<x<=1:raise ValueError('weak-binding domain')
    M=math.asin(math.sqrt(x))-math.sqrt(x*(1-x))
    return period(Q,eps)*(1-M/math.pi)

def run(H,n=24,nr=96):
    e0=H.Qt/H.rt;u,w=leggauss(nr);u=(u+1)/2;w/=2;r=H.rt*u*u;ww=w*2*H.rt*u*r*r;ps=H.psi(r)
    rho=V**2/(4*math.pi*6.6743e-11*R**2);pref=.5*rho*.125*math.sqrt(2)*V*GYR/(H.Ic**2*H.Qt)
    cum=lambda e:pref*(ww@np.array([rate(H,p,e,n)[1] for p in ps]))
    dens=lambda e:pref*(ww@np.array([density(H,p,e,n) for p in ps]))
    esc=cum(e0);slope=dens(e0);rows=[]
    for age in [1.,5.,10.,13.8]:
        emax=e0/2
        eps=emax if exterior_min(H.Qt,H.rt,emax)>=age else brentq(lambda e:exterior_min(H.Qt,H.rt,e)-age,emax*1e-8,emax)
        weak=cum(e0-eps)-esc
        rows.append(dict(age_Gyr=age,eps=eps,eps_over_e0=eps/e0,min_exterior_Gyr=exterior_min(H.Qt,H.rt,eps),birth_per_Gyr=weak,birth_over_escape=weak/esc,linear_estimate=slope*eps,apocenter_lower_kpc=.6*H.Qt/(2*eps)))
    h=e0*1e-4;fd=(cum(e0-h)-cum(e0+h))/(2*h)
    return dict(model=H.model,n=n,nr=nr,rt=H.rt,Qt=H.Qt,e0=e0,escape=esc,slope=slope,fd=fd,fd_rel=abs(fd/slope-1),rows=rows)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--model',default='linear');p.add_argument('--n',type=int,default=24);p.add_argument('--nr',type=int,default=96);a=p.parse_args()
    print(json.dumps(run(Halo(model=a.model),a.n,a.nr)))
