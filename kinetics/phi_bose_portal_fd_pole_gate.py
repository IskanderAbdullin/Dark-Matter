"""FD electron inverse source at the vacuum vector pole, empty dark sector.
Narrow-width approximation only. No time-dependent quantum evolution.
"""
import math,json
import numpy as np
from scipy.integrate import quad
from scipy.special import roots_legendre,expit
from functools import lru_cache
from kinetics.phi_bose_portal_initial_gain_gate import M,M2,A,C,GX,M1,MCHI2,scalar_kinematics

nodes=lru_cache(None)(roots_legendre)
def angular_factor(E,Ep,s,m,T,n=32):
    x,w=nodes(n)
    P=np.sqrt(np.maximum(0.,(E+Ep)**2-s))
    beta=math.sqrt(1.-4*m*m/s)
    psi2=np.clip(np.divide((E-Ep)**2,P**2*beta**2,out=np.zeros_like(P),where=P>0),0.,1.)
    X=x.reshape((1,-1)); ps=psi2.reshape((-1,1)); PP=P.reshape((-1,1)); Es=(E+Ep).reshape((-1,1))
    c2=ps*X**2+(1-ps)*(1-X**2)/2
    aa=1+4*m*m/s; bb=1-4*m*m/s
    weight=expit((Es+PP*X)/(2*T))*expit((Es-PP*X)/(2*T))
    return np.sum((aa+bb*c2)*weight*w,axis=1)/(2*(aa+bb/3))

def source_nwa(E,T,m,fd=True,n=32):
    if E<=m:return 0.
    p=math.sqrt((E-m)*(E+m));b=(M2-2*m*m)/2;B=math.sqrt((b-m*m)*(b+m*m))
    lo=max(m,(b*b+p*p*m*m)/(E*b+p*B));hi=(E*b+p*B)/(m*m)
    if hi<=lo:return 0.
    cap=min((hi-lo)/T,50.)
    if fd:
        x,w=nodes(n)
        u=(1+x)*cap/2;Ep=lo+T*u
        integ=cap/2*np.dot(w,np.exp(-u)*angular_factor(E,Ep,M2,m,T,n))
    else:integ=-math.expm1(-cap)
    return GX*C*T*math.exp(-(E+lo)/T)*math.pi/A*M2*(M2+2*m*m)*integ/(192*math.pi**3*E*p)

def ratio(T,alpha,fd=True,n=32):
    e,lo,hi=scalar_kinematics(T*alpha)
    cap=min((hi-lo)/T,60.)
    j1=quad(lambda u:source_nwa(lo+u*T,T,M1,fd,n)*T,0,cap,epsabs=0.,epsrel=1e-7,limit=120)[0]
    j2=quad(lambda u:source_nwa(lo+u*T-e,T,MCHI2,fd,n)*T,0,cap,epsabs=0.,epsrel=1e-7,limit=120)[0]
    return dict(T=T,alpha=alpha,fd=fd,n=n,J1=j1,J2=j2,ratio=j1/j2)

def rate(T,m,fd=True,n=32):
    return 2/(2*math.pi**2)*quad(lambda u:source_nwa(m+T*u,T,m,fd,n)*(m+T*u)*math.sqrt(T*u*(2*m+T*u))*T,0,80,epsabs=0.,epsrel=2e-7,limit=120)[0]

if __name__=='__main__':
    rows=[]
    for T,alpha in [(.4,.012),(.6882487043,.001),(.6882487043,.012),(1,.001),(1,.012),(2,.001),(5,.0001),(10,.0001)]:
        for fd,n in [(False,32),(True,32),(True,64)]:
            rows.append(ratio(T,alpha,fd,n))
    rates=[dict(T=T,m=m,mb=rate(T,m,False),fd=rate(T,m,True,64)) for T in [1.,5.] for m in [M1,MCHI2]]
    print(json.dumps(dict(rows=rows,rates=rates),indent=2))
