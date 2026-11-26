"""Finite-width FD empty-sector electron source by swapped s and partner integrals."""
import math,json
import numpy as np
from scipy.integrate import quad
from scipy.special import expit
from kinetics.phi_bose_portal_initial_gain_gate import M,M2,A,C,GX,M1,MCHI2,MPHI,scalar_kinematics
from kinetics.phi_bose_portal_fd_pole_gate import nodes

def source_bw(E,T,m,fd=True,n=16):
    if E<=m:return 0.
    p=math.sqrt((E-m)*(E+m))
    epcap=m+60*T;ppcap=math.sqrt((epcap-m)*(epcap+m))
    smax=2*m*m+2*(E*epcap+p*ppcap)
    smin=4*m*m
    x,w=nodes(n)
    marks=[smin,smax]+[smin+E*T*k for k in [.01,.1,1,4,16,64]]+[M2+A*k for k in [-1024,-256,-64,-16,-4,-1,0,1,4,16,64,256,1024]]
    pts=np.unique(np.clip(marks,smin,smax));lo=pts[:-1];hi=pts[1:]
    s=((lo+hi)[:,None]/2+(hi-lo)[:,None]/2*x).reshape(-1)
    ws=((hi-lo)[:,None]/2*w).reshape(-1)
    b=(s-2*m*m)/2;B=np.sqrt(np.maximum(0.,(b-m*m)*(b+m*m)))
    emin=np.maximum(m,(b*b+p*p*m*m)/(E*b+p*B));emax=np.minimum(epcap,(E*b+p*B)/(m*m))
    cap=np.clip((emax-emin)/T,0.,50.)
    if fd:
        up=np.minimum(cap,8.);u1=up[:,None]*(1+x)/2
        rem=cap-up;u2=up[:,None]+rem[:,None]*(1+x)/2
        u=np.concatenate([u1,u2],axis=1)
        wu=np.concatenate([up[:,None]*w/2,rem[:,None]*w/2],axis=1)
        Ep=emin[:,None]+T*u
        P=np.sqrt(np.maximum(0.,(E+Ep)**2-s[:,None]))
        beta2=1-4*m*m/s
        den=P**2*beta2[:,None]
        psi2=np.clip(np.divide((E-Ep)**2,den,out=np.zeros_like(P),where=den>0),0.,1.)
        X=x[None,None,:];ps=psi2[:,:,None]
        c2=ps*X**2+(1-ps)*(1-X**2)/2
        aa=1+4*m*m/s;bb=1-4*m*m/s
        weight=expit((E+Ep[:,:,None]+P[:,:,None]*X)/(2*T))*expit((E+Ep[:,:,None]-P[:,:,None]*X)/(2*T))
        F=np.sum((aa[:,None,None]+bb[:,None,None]*c2)*weight*w,axis=2)/(2*(aa+bb/3))[:,None]
        integ=np.sum(np.exp(-u)*wu*F,axis=1)
    else:integ=-np.expm1(-cap)
    ss=s*(s+2*m*m)/((s-M2)**2+A*A)
    result=np.dot(ws,ss*np.exp(-emin/T)*integ)
    return GX*C*T*math.exp(-E/T)*result/(192*math.pi**3*E*p)

def ratio(T,alpha,fd=True,n=16):
    e,lo,hi=scalar_kinematics(T*alpha)
    cap=min((hi-lo)/T,60.)
    j1=quad(lambda u:source_bw(lo+u*T,T,M1,fd,n)*T,0,cap,epsabs=0.,epsrel=3e-6,limit=90)[0]
    j2=quad(lambda u:source_bw(lo+u*T-e,T,MCHI2,fd,n)*T,0,cap,epsabs=0.,epsrel=3e-6,limit=90)[0]
    return dict(T=T,alpha=alpha,fd=fd,n=n,J1=j1,J2=j2,ratio=j1/j2)

def gain(T,alpha,n=16):
    k=T*alpha;e,lo,hi=scalar_kinematics(k);cap=min((hi-lo)/T,60.)
    r=ratio(T,alpha,True,n);j1,j2=r['J1'],r['J2']
    est=(M1*M1+MPHI*MPHI-MCHI2*MCHI2)/(2*M1);ks=math.sqrt((est-MPHI)*(est+MPHI))
    rest=1e-6*((M1+MCHI2)**2-MPHI**2)*ks/(8*math.pi*M1*M1)
    c1=quad(lambda u:source_bw(lo+u*T,T,M1,True,n)*(lo+u*T)/(rest*M1)*T,0,cap,epsabs=0.,epsrel=3e-6)[0]
    pref=4*rest*M1*M1/(2*ks*k*e);delta=max(0.,(j1-j2)/j1);tcut=8*(c1/j1)*delta/(1-4*delta) if delta<.25 else float('inf')
    r.update(tau=c1/j1,tcut=tcut,upper=min(pref*.5*c1*c1*max(0.,j1-j2)/(j1*j2),pref*.5*max(0.,j1-j2)*tcut*tcut))
    return r

if __name__=='__main__':
    comparison=[ratio(T,alpha,fd,n) for T,alpha in [(.4,.012),(.6882487043,.012),(5,.0001)] for fd,n in [(False,16),(True,12),(True,16),(True,24)]]
    gains=[gain(T,a) for T,a in [(.4,.012),(.6882487043,.0001),(.6882487043,.001),(.6882487043,.012),(1,.001),(2,.001),(5,.0001),(10,.0001)]]
    print(json.dumps(dict(comparison=comparison,gains=gains),indent=2))
