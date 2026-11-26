"""Closed ννbar ↔ Zprime equilibrium for conserved Nnu+2NV and energy.

One Zprime ↔ ν_f+νbar_f, two active ν flavors (4 one-helicity
degrees), three massive-vector polarisations. Both ν and νbar carry
the same *effective particle-number* chemical potential mu, NOT a
lepton asymmetry. μV=2μ from this restricted chemical reaction.
All other weak, e±γ, dark channels, expansion and actual time
integration are OMITTED. Initial FD ν mu=0 and EMPTY vector.
"""
from __future__ import annotations
import math
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import expit,zeta
from scipy.optimize import root


def quadrature(n=256):
    if n<16:raise ValueError
    z,w=leggauss(n)
    return 30*(z+1),30*w


def thermodynamics(T,y,mV,*,n=256):
    """Finite-mass vector + massless ν equilibrium, y=muν/T, muV=2muν."""
    if not(T>0 and mV>0 and math.isfinite(y) and 2*y*T<mV):
        raise ValueError
    q,w=quadrature(n)
    fv=expit(y-q)
    e=np.hypot(q,mV/T)
    bv=1/np.expm1(e-2*y)
    nnu=4*T**3/(2*math.pi**2)*float(np.dot(w,q*q*fv))
    enu=4*T**4/(2*math.pi**2)*float(np.dot(w,q**3*fv))
    nv=3*T**3/(2*math.pi**2)*float(np.dot(w,q*q*bv))
    ev=3*T**4/(2*math.pi**2)*float(np.dot(w,q*q*e*bv))
    return {"nnu":nnu,"enu":enu,"nV":nv,"eV":ev,
            "N_conserved":nnu+2*nv,"E_conserved":enu+ev,
            "S_nu":(4*enu/3-y*T*nnu)/T,
            "S_V":3*T**3/(2*math.pi**2)*float(np.dot(
                w,q*q*((1+bv)*np.log1p(bv)-bv*np.log(bv))))}


def equilibrium_from_empty_vector(mV_over_T0,*,n=256):
    if not(mV_over_T0>0):raise ValueError
    m=mV_over_T0
    initial=thermodynamics(1.,0.,m,n=n)
    n0=initial["nnu"]
    e0=initial["enu"]
    S0=initial["S_nu"]
    def equations(x):
        T=math.exp(x[0]);y=x[1]
        if 2*y*T>=m or not (-20<x[0]<20):
            return [1e5,1e5]
        v=thermodynamics(T,y,m,n=n)
        return [(v["N_conserved"]-n0)/n0,(v["E_conserved"]-e0)/e0]
    sol=root(equations,[math.log(1.246),-1.136],tol=1e-11)
    if not sol.success and max(abs(x) for x in equations(sol.x))>1e-9:
        raise RuntimeError(sol.message)
    T=math.exp(sol.x[0]);y=sol.x[1]
    v=thermodynamics(T,y,m,n=n)
    return {"mV_over_T0":m,"T_final_over_T0":T,
            "mu_nu_over_T_final":y,"mu_V_over_T_final":2*y,
            "n_V_final_over_thermal_mu0_initial":v["nV"]/initial["nV"],
            "nnu_final_over_initial":v["nnu"]/n0,
            "nV_over_initial_nnu":v["nV"]/n0,
            "relative_N_residual":v["N_conserved"]/n0-1,
            "relative_E_residual":v["E_conserved"]/e0-1,
            "S_total_final_over_initial":(v["S_nu"]+v["S_V"])/S0,
            "no_electron_weak_or_friedmann":True,
            "only_closed_Nnu_plus_2NV_charge_assumed":True}


def independent_massless_series_equilibrium():
    """Independent massless-vector polylog fugacity series, no radial quadrature.

    At final negative y sums are absolutely exponentially convergent.
    Exact initial muν=0 Fermi moments from ζ(3),π.
    """
    n0=3*zeta(3,1)/math.pi**2
    e0=7*math.pi**2/60
    k=np.arange(1,129,dtype=float)
    signs=(-1.)**(k+1)
    def vals(T,y):
        f3=float(np.sum(signs*np.exp(k*y)/k**3))
        f4=float(np.sum(signs*np.exp(k*y)/k**4))
        b3=float(np.sum(np.exp(2*k*y)/k**3))
        b4=float(np.sum(np.exp(2*k*y)/k**4))
        nnu=4*T**3/math.pi**2*f3
        enu=12*T**4/math.pi**2*f4
        nv=3*T**3/math.pi**2*b3
        ev=9*T**4/math.pi**2*b4
        return nnu,enu,nv,ev
    def fun(v):
        T=math.exp(v[0]);y=v[1]
        nnu,enu,nv,ev=vals(T,y)
        return [(nnu+2*nv)/n0-1,(enu+ev)/e0-1]
    r=root(fun,[math.log(1.246),-1.136],tol=1e-12)
    if not r.success and max(abs(v) for v in fun(r.x))>1e-10:
        raise RuntimeError
    T=math.exp(r.x[0]);y=r.x[1]
    nnu,enu,nv,ev=vals(T,y)
    return {"T_final_over_T0":T,"mu_nu_over_T_final":y,
            "nV_over_initial_nnu":nv/n0,
            "relative_number_error":(nnu+2*nv)/n0-1,
            "relative_energy_error":(enu+ev)/e0-1,
            "independent_massless_fugacity_series":True}


if __name__=="__main__":
    import json
    print(json.dumps({str(m):equilibrium_from_empty_vector(m)
                      for m in (.03,.1,.3)},indent=2))
