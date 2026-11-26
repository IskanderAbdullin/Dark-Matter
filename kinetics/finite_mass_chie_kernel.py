"""Massive-lepton generalization of the exact invariant chi+lepton energy kernel.
Tree-level Dirac vector contact EFT, NOT full thermal SM or cosmological rate.
Units: masses/energies GeV, G GeV^-2, kernel GeV.
"""
import math
import numpy as np
from numpy.polynomial.legendre import leggauss

def kernel(M,m,G,E,w,Ep,wp,nq=24,gchi=4.,glep=4.):
    if not (M>0 and m>=0 and G>=0 and gchi>0 and glep>0 and isinstance(nq,int) and nq>=2):
        raise ValueError("invalid physical parameters")
    Ep=np.asarray(Ep,dtype=float);wp=np.asarray(wp,dtype=float)
    if (Ep.shape!=wp.shape or not np.all(np.isfinite(Ep)) or
        not np.all(np.isfinite(wp)) or E<=M or w<=m):
        raise ValueError("invalid energies")
    if np.any(abs(E+w-Ep-wp)>1e-10*max(1,E+w)):
        raise ValueError("energy mismatch")
    p=math.sqrt(E*E-M*M);k=math.sqrt(w*w-m*m)
    pp=np.sqrt(np.maximum(0,Ep*Ep-M*M));kp=np.sqrt(np.maximum(0,wp*wp-m*m))
    lo=np.maximum(abs(p-pp),abs(k-kp));hi=np.minimum(p+pp,k+kp)
    valid=(hi>lo)&(Ep>M)&(wp>m)
    a=np.where(valid,lo,0.);b=np.where(valid,hi,0.)
    x,weights=leggauss(nq)
    q=a[...,None]+(b-a)[...,None]*(x+1)/2
    qs=np.where(q>0,q,1.)
    cp=np.clip((p*p+qs*qs-pp[...,None]**2)/(2*p*qs),-1,1)
    ce=np.clip((kp[...,None]**2-k*k-qs*qs)/(2*k*qs),-1,1)
    t=(E-Ep[...,None])**2-qs**2
    D0=2*(E*w-p*k*cp*ce)
    D2=D0**2+2*p*p*k*k*(1-cp**2)*(1-ce**2)
    amp=2*G**2*(2*D2+2*D0*t+t*t+2*(M*M+m*m)*t)
    if np.any(amp[valid]<-1e-12):
        raise ArithmeticError("negative amplitude")
    return valid*gchi*glep/(256*math.pi**5)*np.sum(np.maximum(0,amp)*weights,axis=-1)*(b-a)/2
