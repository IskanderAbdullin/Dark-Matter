"""Exact angular reduction for vacuum spin-correlated vector 2->2 MB/FD gate."""
import json
from math import exp
import numpy as np
from scipy.special import expit, roots_legendre
from timelike_portal_onsager_coverage_reference import mchi,me,V

X,W=roots_legendre(36)
W=W/2
P2=(3*X**2-1)/2

def moments(E,q0,T,eta,m):
    assert E>2*m and q0>=E and T>0
    v2=1-4*m*m/E**2
    q=np.sqrt((q0-E)*(q0+E))
    eps=q*np.sqrt(v2)*X
    w=expit((q0+eps)/(2*T)-eta)*expit((q0-eps)/(2*T)-eta)
    return w@W,w@(P2*W),v2

def angular(E,q0,Td,Ts,eta,m=me):
    a,b,vi2=moments(E,q0,Td,eta,mchi)
    c,d,vf2=moments(E,q0,Ts,0.,m)
    # Tr[(I-vi^2 ni ni^T)(I-vf^2 nf nf^T)]. The azimuth average
    # of P2(ni.nf) has only the m=0 term, by axial symmetry.
    coef=(2/3)*vi2*vf2/(3-vi2-vf2+vi2*vf2/3)
    z=a*c+coef*b*d
    ld=2*eta-q0/Td
    ls=-q0/Ts
    ma=max(ld,ls)
    delta=2*eta+q0*((Td-Ts)/(Td*Ts))
    net=np.sign(delta)*exp(ma)*(-np.expm1(-abs(delta)))*z
    return net,z,net*delta,coef*b*d

def direct(E,q0,Td,Ts,eta,m=me):
    # Independently integrate both lab polar angles and relative azimuth,
    # explicitly subtracting f_i f_ibar B_f B_fbar and inverse.
    def wi(T,eta_,m_):
        v2=1-4*m_**2/E**2
        q=np.sqrt((q0-E)*(q0+E))
        eps=q*np.sqrt(v2)*X
        rd=np.stack((q0+eps,q0-eps)).astype(np.longdouble)/(2*np.longdouble(T))-np.longdouble(eta_)
        return np.prod(1/(1+np.exp(rd)),axis=0),np.prod(1/(1+np.exp(-rd)),axis=0),v2
    fi,bi,v=wi(Td,eta,mchi)
    ff,bf,u=wi(Ts,0.,m)
    si=np.sqrt(1-X*X)
    phi=np.arange(24)*2*np.pi/24
    dot=X[:,None,None]*X[None,:,None]+si[:,None,None]*si[None,:,None]*np.cos(phi)
    norm=3-v-u+v*u/3
    k=np.mean(3-v-u+v*u*dot*dot,axis=2)/norm
    weights=(W[:,None]*W[None,:])*k
    return np.sum(weights*(fi[:,None]*bf-bi[:,None]*ff)),np.sum(weights*(bi[:,None]*bf))

def test():
    rows=[]
    for E in (.045,.3,1.,V,11.):
        for m in (me,.1056583755,1.77686):
            if E<=2*m: continue
            for q0_ratio in (1.,1.6,4.):
                q0=E*q0_ratio
                for Td,Ts,eta in ((.7,.63,.03),(.68,.71,-.05),(.69,.69,0.),(.69,.69+1e-6,0.)):
                    n,w,p,corr=angular(E,q0,Td,Ts,eta,m)
                    nr,wr=direct(E,q0,Td,Ts,eta,m)
                    e=max(abs(w/wr-1),abs(n-nr)/max(abs(n),1e-200)) if n!=0 else abs(w/wr-1)
                    assert e<4e-9 and p>=0., (E,m,q0_ratio,Td,Ts,eta,e,n,nr)
                    if Td==Ts and eta==0.: assert n==0.
                    rows.append((e,corr/w,w))
    return dict(events=len(rows),max_relative_angular_check=float(max(r[0] for r in rows)),
                max_abs_spin_correlation_fraction=max(abs(r[1]) for r in rows),
                minimum_total_blocking_factor=min(r[2] for r in rows),
                scope="Local vacuum matrix element with exact thermal FD Pauli factors for each bath; no medium propagators, rates, hadronic final states, expansion or relic inference.")
if __name__=="__main__":
    print(json.dumps(test(),indent=2,allow_nan=False))
