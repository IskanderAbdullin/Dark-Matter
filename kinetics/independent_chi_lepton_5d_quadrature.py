"""Independent 5D relativistic finite-vector χℓ Pauli weak quadrature.

A separate implementation of the existing project's on-shell CM event,
intended to catch discrepancies in numerical rate integration; no physical
strong/time solver or cosmological Γ_el or 86-eV DM conclusion.
"""
from __future__ import annotations
import math
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import expit


def occupations(shape):
    if shape == "smooth":
        return (lambda E: expit(-np.asarray(E)/.9+.15*np.exp(-np.asarray(E)/.9)),
                lambda E: expit(-np.asarray(E)+.05*np.exp(-np.asarray(E))))
    if shape == "oscillatory":
        return (lambda E: expit(-np.asarray(E)+.3*np.sin(4*np.asarray(E))),
                lambda E: expit(-np.asarray(E)+.3*np.cos(4*np.asarray(E))))
    if shape == "fd_cold":
        return (lambda E: expit(-np.asarray(E)/.9),
                lambda E: expit(-np.asarray(E)))
    if shape == "fd_eq_chem":
        return (lambda E: expit(.2-np.asarray(E)),
                lambda E: expit(-.1-np.asarray(E)))
    raise ValueError(shape)


def legendre(a,b,n):
    x,w=leggauss(n)
    return a+(b-a)*(x+1)/2,w*(b-a)/2


def evaluate(shape="smooth", *, M=.1, m=0., g=.001, mv=.1,
             nr=4, nmu=4, nc=12, nphi=12, cutoff=12.):
    """Incoming p,k,mu + outgoing CM theta,phi on exactly on-shell 2->2."""
    if not (M>0 and m>=0 and g>0 and mv>0 and nr>=3 and
            nmu>=3 and nc>=3 and nphi>=8 and cutoff>0):
        raise ValueError("nonphysical parameters")
    f,h=occupations(shape)
    pmax=math.sqrt((M+cutoff)**2-M*M)
    kmax=math.sqrt((m+cutoff)**2-m*m)
    px,wp=legendre(0,pmax,nr)
    kx,wk=legendre(0,kmax,nr)
    mux,wmu=legendre(-1,1,nmu)
    cx,wc=leggauss(nc)
    ph=2*math.pi*np.arange(nphi)/nphi
    c=cx[:,None]
    sin=np.sqrt(np.maximum(0,1-c*c))
    cphi=np.cos(ph)[None,:]
    heat=entropy=0.
    neg_event_count=0
    max_energy_defect=0.
    for p,ap in zip(px,wp):
      E=math.hypot(p,M)
      for k,ak in zip(kx,wk):
        omega=math.hypot(k,m)
        for mu,au in zip(mux,wmu):
          tot=E+omega
          P=math.sqrt(p*p+k*k+2*p*k*mu)
          s=tot*tot-P*P
          if P<=0 or s<=(M+m)**2: continue
          ps2=(s-(M+m)**2)*(s-(M-m)**2)/(4*s)
          if ps2<=0: continue
          ps=math.sqrt(ps2)
          root=math.sqrt(s)
          beta=P/tot
          gamma=tot/root
          projection=p*(p+k*mu)/P
          ca=max(-1.,min(1.,gamma*(projection-beta*E)/ps))
          Estar=(s+M*M-m*m)/(2*root)
          aligned=ca*c+math.sqrt(max(0,1-ca*ca))*sin*cphi
          Ep=gamma*(Estar+beta*ps*aligned)
          op=tot-Ep
          t=-2*ps2*(1-c)
          D=s-M*M-m*m
          U=-D-t
          amp=2*g*g*(D*D+U*U+2*(M*M+m*m)*t)*(mv*mv/(mv*mv-t))**2
          fe=float(f(E))
          ho=float(h(omega))
          fp=f(Ep)
          hp=h(op)
          forward=fe*ho*(1-fp)*(1-hp)
          reverse=fp*hp*(1-fe)*(1-ho)
          J=forward-reverse
          if np.any(forward<=0) or np.any(reverse<=0):
              raise ArithmeticError("underflow in test distribution")
          entropy_micro=J*np.log(forward/reverse)
          neg_event_count+=int(np.count_nonzero(entropy_micro<-1e-26))
          de=Ep-E
          max_energy_defect=max(max_energy_defect,float(np.max(abs(de+op-omega))))
          vrel=math.sqrt(max(0,(E*omega-p*k*mu)**2-M*M*m*m))/(E*omega)
          coeff=ap*p*p*ak*k*k*au*vrel/s
          heat+=coeff*np.sum(wc[:,None]*amp*de*J)/nphi
          entropy+=coeff*np.sum(wc[:,None]*amp*entropy_micro)/nphi
    pref=16/(512*math.pi**5)
    return {"shape":shape,"M":M,"m":m,"nr":nr,"nmu":nmu,
            "nc":nc,"nphi":nphi,"cutoff":cutoff,
            "Qchi":float(heat*pref),"Qlepton":float(-heat*pref),
            "entropy":float(entropy*pref),
            "negative_entropy_events":neg_event_count,
            "max_energy_defect":max_energy_defect,
            "synthetic_weak_only":True}


if __name__=="__main__":
    import json
    for shape in ("fd_cold","fd_eq_chem","smooth","oscillatory"):
        print(json.dumps(evaluate(shape),ensure_ascii=False))
