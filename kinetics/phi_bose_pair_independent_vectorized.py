"""Independent vectorized invariant-trace Yukawa χ2χbar2→φφ scattering.

This implementation uses center-of-mass Minkowski contractions and a
logarithmic tanh angular transform, independently of the production
scalar dot-product loops. Output is a positive vacuum spin-averaged σ.
An MB thermal forward reaction rate is a diagnostic, NOT equilibrium
quantum collision kinetics, full SM, FRW, or a 86eV relic prediction.
"""
import json
import math
from functools import lru_cache
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad
from scipy.special import kv

MHEAVY=.1
MLIGHT=.02
MSCALAR=86e-9
Y=.001
MPLANCK=1.2209e19

@lru_cache(maxsize=16)
def nodes(n):
    z,w=leggauss(n)
    v=14*z
    c=np.tanh(v)
    return c,14*w*(1-c*c)

def sigma_independent(s,angle_order=128,y=Y):
    """Physical χ2+antiχ2→2 identical real φ vacuum tree σ in GeV^-2."""
    if s<=4*MLIGHT**2:return 0.
    e=math.sqrt(s)/2
    a=math.sqrt(e*e-MLIGHT**2)
    p=math.sqrt(e*e-MSCALAR**2)
    c,w=nodes(angle_order)
    qz=[a-p*c,a+p*c]
    qsq=[-a*a-p*p+2*a*p*c,-a*a-p*p-2*a*p*c]
    dab=e*e+a*a
    accum=np.zeros_like(c)
    for i,zi in enumerate(qz):
        for j,zj in enumerate(qz):
            daq=a*zi
            dar=a*zj
            dqb=-daq
            dbr=-dar
            qr=qsq[i] if i==j else (p*p-a*a)
            four=daq*dbr-dab*qr+dar*dqb
            two=(MLIGHT*MHEAVY*(daq+dar-dqb-dbr)
                 +MHEAVY*MHEAVY*dab-MLIGHT*MLIGHT*qr)
            tr=four+two-MLIGHT*MLIGHT*MHEAVY*MHEAVY
            accum+=tr/((qsq[i]-MHEAVY*MHEAVY)*(qsq[j]-MHEAVY*MHEAVY))
    return float(y**4*np.dot(w,accum)*p/(64*math.pi*s*a))

def reaction_density_MB(T,angle_order=128,y=Y):
    """Physical forward MB rate and total incoming-pair energy moment.

    γ=4 T/(32π^4) int_{4m2²} ds sigma(s)(s-4m2²) sqrt(s) K1.
    Q_event = -∂γ/∂(1/T) at fixed Yukawa masses and y.
    """
    if T<=0:raise ValueError("T positive")
    xmin=2*MLIGHT/T
    def pref(u):
        s=(T*u)**2
        return 2*T**5*u*u*(u*u-xmin*xmin)*sigma_independent(s,angle_order,y)
    gamma=quad(lambda u:pref(u)*kv(1,u),
               xmin,40.,epsabs=1e-33*(y/Y)**4,epsrel=3e-8)[0]
    energy=quad(lambda u:pref(u)*T*(2*kv(1,u)+u*kv(0,u)),
                xmin,40.,epsabs=1e-32*(y/Y)**4,epsrel=3e-8)[0]
    coeff=4*T/(32*math.pi**4)
    return dict(T_GeV=T,gamma_forward_GeV4=coeff*gamma,
                pair_energy_forward_GeV5=coeff*energy,
                mean_total_energy_per_event_GeV=energy/gamma,
                scattering_order=angle_order)

def ratio_to_H(T,gstar,angle_order=128):
    x=reaction_density_MB(T,angle_order)
    density=2*MLIGHT**2*T*kv(2,MLIGHT/T)/(2*math.pi**2)
    H=1.66*math.sqrt(gstar)*T*T/MPLANCK
    return x["gamma_forward_GeV4"]/density/H

def audit():
    rows=[]
    for order in (32,48,64,96,128):
        rows.append(dict(angle_order=order,
                         gamma_T0p1=reaction_density_MB(.1,order)["gamma_forward_GeV4"],
                         gamma_T1=reaction_density_MB(1.,order)["gamma_forward_GeV4"],
                         gamma_T100=reaction_density_MB(100.,order)["gamma_forward_GeV4"]))
    return dict(schema_version=1,
                scope="independent physically normalized vacuum Yukawa pair trace/cross section and forward MB thermal event energy, not full FD BE kinetic source",
                angular_refinement=rows,
                final_128_Hubble_ratios={
                    "T0p1_g17p25":ratio_to_H(.1,17.25),
                    "T1_g61p75":ratio_to_H(1.,61.75),
                    "T100_g106p75":ratio_to_H(100.,106.75)},
                class_ready=False)

if __name__=="__main__":
    print(json.dumps(audit(),indent=2,allow_nan=False))
