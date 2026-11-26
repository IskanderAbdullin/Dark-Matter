"""Independent adaptive angular weak heat for one fixed chi+lepton event.

Not the full 5D integrated collision / gamma_el. At fixed incoming CM
kinematics, angle-integrated full-Pauli energy transfer is evaluated by:
(1) adaptive scipy QUAD in cos(theta*), (2) Gauss in cosine,
(3) Gauss in rational u=|t|/(mV**2+|t|), with exact Jacobian.
The outgoing Fermi factors are recomputed at EVERY angle node.
"""
from __future__ import annotations
import math
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad
from scipy.special import expit


def integrand_factory(*, M=.1, m=0., p=.9, k=1.2, mu=-.35,
                      Te=1., Tchi=.9, G=.001, mV=.03, nphi=32):
    if not (M>0 and m>=0 and p>0 and k>0 and -1<mu<1 and
            Te>0 and Tchi>0 and G>0 and mV>0 and nphi>=8):
        raise ValueError("physical incoming state, mediator and temperatures")
    E=math.hypot(p,M)
    omega=math.hypot(k,m)
    Et=E+omega
    P=math.sqrt(p*p+k*k+2*p*k*mu)
    s=Et*Et-P*P
    root=math.sqrt(s)
    pstar2=(s-(M+m)**2)*(s-(M-m)**2)/(4*s)
    if pstar2<=0:
        raise ValueError("open scattering phase space required")
    pstar=math.sqrt(pstar2)
    beta=P/Et
    gamma=Et/root
    parallel=p*(p+k*mu)/P
    ca=max(-1.,min(1.,gamma*(parallel-beta*E)/pstar))
    Estar=(s+M*M-m*m)/(2*root)
    D=s-M*M-m*m
    phi=np.linspace(0,2*math.pi,nphi,endpoint=False)
    cphi=np.cos(phi)
    f=expit(-E/Tchi)
    g=expit(-omega/Te)

    def angular(c):
        c=float(c)
        if abs(c)>1+1e-12:
            raise ValueError("outside CM angle")
        c=float(np.clip(c,-1.,1.))
        t=-2*pstar2*(1-c)
        U=-D-t
        amplitude=2*G*G*(D*D+U*U+2*(M*M+m*m)*t)*(
            mV*mV/(mV*mV-t))**2
        cp=ca*c+math.sqrt(max(0.,1-ca*ca))*math.sqrt(max(0.,1-c*c))*cphi
        Ep=gamma*(Estar+beta*pstar*cp)
        wp=Et-Ep
        fp=expit(-Ep/Tchi)
        gp=expit(-wp/Te)
        F=f*g*(1-fp)*(1-gp)
        R=fp*gp*(1-f)*(1-g)
        return float(amplitude*np.mean((Ep-E)*(F-R)))
    return angular,pstar2


def angular_integrals(mV=.03, *,
                      orders=(16,24,32,48,64,80,96,128,160), **kw):
    fun,ps2=integrand_factory(mV=mV,**kw)
    # scipy QUAD is a third angular algorithm independent of Legendre.
    reference,abserr=quad(fun,-1.,1.,epsabs=1e-24,
                          epsrel=2e-12,limit=400)
    umax=4*ps2/(mV*mV+4*ps2)
    rows=[]
    for n in orders:
        z,w=leggauss(n)
        direct=float(np.dot(w,[fun(c) for c in z]))
        u=umax*(z+1)/2
        c=1-mV*mV*u/(2*ps2*(1-u))
        jac=mV*mV/(2*ps2*(1-u)**2)
        mapped=float(np.dot(w,[fun(ci)*ji for ci,ji in zip(c,jac)])*umax/2)
        rows.append({"order":int(n),
                     "cosine_relative_error":direct/reference-1,
                     "transfer_relative_error":mapped/reference-1})
    return {"mV":mV,"pstar_squared":ps2,"quad_reference":reference,
            "adaptive_quad_error_estimate":abserr,"scan":rows,
            "not_full_five_dimensional_integrated_heat":True}


def propagator_angular_identity(ps2,mV,n=32):
    """Independent exact angular integral of squared propagator alone."""
    if ps2<=0 or mV<=0:
        raise ValueError("positive CM momentum and mediator mass")
    umax=4*ps2/(mV*mV+4*ps2)
    z,w=leggauss(n)
    # integrand propagator(1-u)^2 times dc/du -> constant mV²/(2p*²).
    mapped=float(np.sum(w)*umax/2*mV*mV/(2*ps2))
    analytic=2*mV*mV/(mV*mV+4*ps2)
    return mapped,analytic


if __name__=="__main__":
    import json
    print(json.dumps({str(mV):angular_integrals(mV)
                      for mV in (.03,.1,.3)},indent=2))
