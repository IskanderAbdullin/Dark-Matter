"""Discrete microreversible Z' <-> nu + antinu quantum event graph.

NUMERICAL PILOT (NOT full continuum strong Boltzmann):
- Equal energy lattice E_V,j=j*h, k_nu,i=i*h, k_antinu,l=l*h, j=i+l.
- Exact massless two-body angular support |2k-E|<=sqrt(E²-mV²).
- Independently evolving f_nu, f_antinu, bV WITHOUT FD thermal reset.
- Eventwise scalar backward-Euler root strictly inside Pauli/Bose bounds.
- Exact discrete L=Nnu-Nantinu, C=Nnu+Nantinu+2NV, energy and H theorem.
- Positive approximate physical angular weights; finite-grid aliasing remains.
Spin degeneracies gV=3, gnu=2 and g_antinu=2 (two active flavors).
The E lattice has endpoint quadrature bias; it is NOT a physical
gauge/FRW/chi1-chi2-phi/86eV DM solver or calibrated Gamma_el.
"""
from __future__ import annotations

import math
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import brentq
from scipy.special import expit, xlogy

GV=3
G_NU=2

def graph(m=0.4,h=0.2,emax=4.,gamma0=1.):
    if not(all(math.isfinite(z) and z>0 for z in (m,h,emax,gamma0))
           and emax>m and abs(m/h-round(m/h))<1e-10
           and abs(emax/h-round(emax/h))<1e-10):
        raise ValueError("positive lattice with integer m/h, Emax/h")
    jv=np.arange(int(round(m/h))+1,int(round(emax/h))+1,dtype=int)
    E=jv*h
    p=np.sqrt(E*E-m*m)
    k=np.arange(1,int(round(emax/h))+1,dtype=float)*h
    wV=GV/(2*math.pi**2)*E*p*h
    wnu=G_NU/(2*math.pi**2)*k*k*h
    edges=[]
    for v,(j,pp) in enumerate(zip(jv,p)):
        for i in range(1,int(j)):
            if abs((2*i-j)*h)<=pp+1e-13:
                # A vector at rest has vanishing shell measure here.
                # dcos(theta)/2=dk/p, positive angular event weight.
                rate=wV[v]*gamma0*m/E[v]*h/pp
                edges.append((v,i-1,int(j-i-1),rate))
    return {"E":E,"p":p,"k":k,"wV":wV,"wnu":wnu,
            "edges":edges,"m":m,"h":h,"emax":emax,"gamma0":gamma0}

def synthetic_spectra(g):
    """Prescribed nonthermal initial occupations; no physical cosmology."""
    k=g["k"];E=g["E"]
    fn=expit(-k+0.15*np.exp(-k))
    fa=expit(-k-0.10*np.exp(-0.5*k))
    b=.75/np.expm1(E)
    return b,fn,fa

def ledger(g,b,fn,fa):
    E,k,wV,wnu=(g[a] for a in ("E","k","wV","wnu"))
    NV=float(np.dot(wV,b))
    Nn=float(np.dot(wnu,fn))
    Na=float(np.dot(wnu,fa))
    energy=float(np.dot(wV,E*b)+np.dot(wnu,k*(fn+fa)))
    Sb=float(np.dot(wV,(1+b)*np.log1p(b)-xlogy(b,b)))
    Sf=float(np.dot(wnu,-xlogy(fn,fn)-xlogy(1-fn,1-fn)
                        -xlogy(fa,fa)-xlogy(1-fa,1-fa)))
    return {"reduced_charge":Nn+Na+2*NV,"lepton_asymmetry":Nn-Na,
            "energy":energy,"entropy":Sb+Sf,
            "NV":NV,"Nnu":Nn,"Nantinu":Na}

def event_moments(g,b,fn,fa):
    """Integrated discrete reference rates, energy and entropy production."""
    E=g["E"]
    dNV=dEV=dS=0.
    for v,i,l,W in g["edges"]:
        F=b[v]*(1-fn[i])*(1-fa[l])
        R=(1+b[v])*fn[i]*fa[l]
        J=F-R
        q=W*J
        dNV-=q
        dEV-=E[v]*q
        if F>0 and R>0:
            dS+=J*W*math.log(F/R)
        elif F!=R:
            dS=math.inf
    return {"dNV":dNV,"dEV":dEV,"entropy_production":dS,
            "dNnu":-dNV,"dNantinu":-dNV,
            "dEnu_plus_antinu":-dEV}

def sweep(g,b,fn,fa,dt):
    """Implicit elementary reaction extent; NO clipping/thermal reset.

    For every microreversible quartet, F(x) strictly decreases and R(x)
    increases as reaction extent x increases. Thus x-dt*W*(F-R)=0 has
    a unique root between physical Pauli/Bose bounds; brentq solves it.
    Every event changes one V and one nu plus one antinu in exact
    number/energy ratios, independently of the finite step size.
    """
    if not(math.isfinite(dt) and dt>=0):
        raise ValueError("finite nonnegative dt")
    b=np.array(b,float,copy=True)
    fn=np.array(fn,float,copy=True)
    fa=np.array(fa,float,copy=True)
    wV=g["wV"]; wnu=g["wnu"]
    if (b.shape!=wV.shape or fn.shape!=wnu.shape or
        fa.shape!=wnu.shape or not np.all(np.isfinite(b)) or
        not np.all(np.isfinite(fn)) or not np.all(np.isfinite(fa)) or
        np.any(b<0) or np.any(fn<0) or np.any(fa<0) or
        np.any(fn>1) or np.any(fa>1)):
        raise ValueError("occupations must obey exact finite Pauli/Bose limits")
    if dt==0:return b,fn,fa
    for v,i,l,W in g["edges"]:
        lo=max(-fn[i]*wnu[i],-fa[l]*wnu[l])
        hi=min(b[v]*wV[v],(1-fn[i])*wnu[i],(1-fa[l])*wnu[l])
        def residual(x):
            bv=b[v]-x/wV[v]
            f=fn[i]+x/wnu[i]
            a=fa[l]+x/wnu[l]
            F=bv*(1-f)*(1-a)
            R=(1+bv)*f*a
            return x-dt*W*(F-R)
        x=brentq(residual,lo,hi,xtol=5e-324,rtol=1e-14)
        b[v]-=x/wV[v]
        fn[i]+=x/wnu[i]
        fa[l]+=x/wnu[l]
    return b,fn,fa

def evolve(g,b,fn,fa,total_t=4.,dt=.25):
    if not(total_t>=0 and dt>0 and math.isfinite(total_t)
           and math.isfinite(dt)):raise ValueError
    steps=int(round(total_t/dt))
    if abs(steps*dt-total_t)>1e-10:
        raise ValueError("choose dt dividing total_t")
    begin=ledger(g,b,fn,fa)
    min_entropy_change=math.inf
    last=begin["entropy"]
    for _ in range(steps):
        b,fn,fa=sweep(g,b,fn,fa,dt)
        ss=ledger(g,b,fn,fa)["entropy"]
        min_entropy_change=min(min_entropy_change,ss-last)
        last=ss
    after=ledger(g,b,fn,fa)
    return b,fn,fa,{"steps":steps,"dt":dt,
       "final_b_lowest_energy":float(b[0]),
       "max_abs_reduced_charge_error":abs(after["reduced_charge"]-begin["reduced_charge"]),
       "max_abs_lepton_asymmetry_error":abs(after["lepton_asymmetry"]-begin["lepton_asymmetry"]),
       "max_abs_energy_error":abs(after["energy"]-begin["energy"]),
       "entropy_gain":after["entropy"]-begin["entropy"],
       "min_step_entropy_gain":min_entropy_change,
       "min_b":float(np.min(b)),"min_fn":float(np.min(fn)),
       "min_fantinu":float(np.min(fa)),
       "max_fn":float(np.max(fn)),"max_fantinu":float(np.max(fa))}

def independent_continuum_weak(m,emax,*,nr=144,nangle=144):
    """Separate continuum vector p,CM angle weak N/E/entropy integration.

    Same synthetic initial shapes. Only used as a finite-grid error
    control, NOT a full physical Z' plasma or continuum time evolution.
    """
    z,w=leggauss(nr)
    pmax=math.sqrt(emax*emax-m*m)
    p=pmax*(z+1)/2
    wp=w*pmax/2
    c,wc=leggauss(nangle)
    E=np.hypot(p,m)
    k1=(E[:,None]+p[:,None]*c[None,:])/2
    k2=E[:,None]-k1
    f=expit(-k1+.15*np.exp(-k1))
    a=expit(-k2-.1*np.exp(-.5*k2))
    b=.75/np.expm1(E[:,None])
    F=b*(1-f)*(1-a)
    R=(1+b)*f*a
    wV=GV/(2*math.pi**2)*wp*p*p*m/E
    q=(F-R)@wc/2
    return {"dNV":-float(np.dot(wV,q)),
            "dEV":-float(np.dot(wV,E*q)),
            "entropy_production":float(np.dot(wV,((F-R)*np.log(F/R))@wc/2))}

if __name__=="__main__":
    import json
    rows=[]
    ref=independent_continuum_weak(.4,4)
    for h in (.4,.2,.1,.05):
        g=graph(.4,h,4)
        init=synthetic_spectra(g)
        x=event_moments(g,*init)
        rows.append({"h":h,"edges":len(g["edges"]),
          "rate_relative_to_continuum":x["dNV"]/ref["dNV"]-1,
          "energy_relative_to_continuum":x["dEV"]/ref["dEV"]-1,
          "entropy_relative_to_continuum":x["entropy_production"]/ref["entropy_production"]-1})
    g=graph(.4,.2,4)
    init=synthetic_spectra(g)
    trials=[evolve(g,*init,total_t=4,dt=dt)[3] for dt in (1,.5,.25,.125)]
    print(json.dumps({"continuum_weak":ref,"grid_scan":rows,
      "time_scan":trials,"class_ready":False},indent=2))
