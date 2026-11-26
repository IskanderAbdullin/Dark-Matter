"""DIRECT microreversible on-shell χ+e 2→2 Pauli collision on an energy lattice.

This is a DISTINCT discretisation from paired_energy_quartets.py:
No post-hoc Wraw/Wreverse averaging or outgoing energy interpolation.
Integrate the exact Lorentz-invariant four-body collision phase space
analytically over five angles to a ONE-DIMENSIONAL momentum transfer q
quadrature at fixed Eχ,k,Eχ',k'=Eχ+k-Eχ'. The resulting bare W is positive,
exactly symmetric under incoming/outgoing exchange, and contains the
actual vector-contact amplitude. Cell-center energy conservation is exact.

This is energy-COLLOCATION quadrature for the continuum operator and not
an exact finite-bin phase-space representation. It must be tested under
energy-spacing, q, cutoff, time and angular-independent full 5D refinement.
Massless e±, below-EW vector contact, fixed volume. Not physical gamma_el,
FRW, real full-SM plasma, Bose 1↔2 or a relic of 86 eV dark matter.
"""
from __future__ import annotations
import math
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import expit
from scipy.optimize import brentq


def invariant_kernel(M,G,E,k,Ep,kp,gchi=4.,ge=4.,nq=12):
    """dR/dE dk dE' for exact 2->2 phase space, symmetric time reversal.

    This is gχ ge/(256π⁵) ∫_{qlo}^{qhi}dq <|M(s,t)|²>_azimuth.
    Electron is massless. E,k scalar; Ep,kp same-shaped arrays.
    The common energy delta is already used to eliminate dk'.
    """
    z=(M,G,E,k,gchi,ge)
    if (not all(math.isfinite(float(t)) for t in z) or
         M<=0 or G<0 or E<=M or k<=0 or gchi<=0 or ge<=0 or
         type(nq) is not int or nq<2):
        raise ValueError("physical scattering input and nq>=2 required")
    Ep=np.asarray(Ep,float);kp=np.asarray(kp,float)
    if Ep.shape!=kp.shape or np.any(~np.isfinite(Ep)) or np.any(~np.isfinite(kp)):
        raise ValueError("matching finite outgoing energies required")
    if np.any(np.abs(E+k-Ep-kp)>1e-10*max(1.,E+k)):
        raise ValueError("incoming/outgoing total energy must match")
    p=math.sqrt(E*E-M*M)
    pp=np.sqrt(np.maximum(0.,Ep*Ep-M*M))
    low=np.maximum(abs(p-pp),abs(k-kp))
    high=np.minimum(p+pp,k+kp)
    ok=(high>low)&(Ep>M)&(kp>0)
    lo=np.where(ok,low,0.);hi=np.where(ok,high,0.)
    x,w=leggauss(nq)
    q=lo[...,None]+(hi-lo)[...,None]*(x+1)/2
    qs=np.where(q>0,q,1.)
    cp=np.clip((p*p+qs*qs-pp[...,None]**2)/(2*p*qs),-1.,1.)
    ce=np.clip((kp[...,None]**2-k*k-qs*qs)/(2*k*qs),-1.,1.)
    # χ momentum transfer q=p−p′=k′−k; in COM t=(E−E′)²−q².
    t=(E-Ep[...,None])**2-qs**2
    D0=2*(E*k-p*k*cp*ce)
    # Average over the relative χ/e cone azimuth ψ exactly:
    # <D²>=D0²+2 p²k²(1−cp²)(1−ce²), D=s−M².
    D2=D0*D0+2*p*p*k*k*(1-cp*cp)*(1-ce*ce)
    m2=2*G*G*(2*D2+2*D0*t+t*t+2*M*M*t)
    if np.any(m2[ok]<-1e-10):
        raise ArithmeticError("negative tree-level physical amplitude")
    integral=np.sum(np.maximum(0.,m2)*w,axis=-1)*(hi-lo)/2
    return gchi*ge/(256*math.pi**5)*integral*ok


def microreversible_lattice(*,M=.1,T=1.,G=.001,dE=.5,
                             nchi=44,ne=44,gchi=4.,ge=4.,nq=12):
    if (type(nchi) is not int or type(ne) is not int or nchi<3 or ne<3
        or nchi*ne>45000 or not all(math.isfinite(float(t)) for t in (M,T,G,dE,gchi,ge))
        or M<=0 or T<=0 or G<0 or dE<=0 or gchi<=0 or ge<=0):
        raise ValueError("finite inputs, moderate grids and dE>0 required")
    E=M+(np.arange(nchi)+.5)*dE
    k=(np.arange(ne)+.5)*dE
    p_edges=np.sqrt(np.maximum(0,(M+dE*np.arange(nchi+1))**2-M*M))
    Vchi=gchi/(6*math.pi**2)*np.diff(p_edges**3)
    Ve=ge*dE**3/(6*math.pi**2)*np.diff(np.arange(ne+1,dtype=float)**3)
    i_all=[];a_all=[];j_all=[];b_all=[];w_all=[]
    for i,e in enumerate(E):
        for a,ka in enumerate(k):
            j=np.arange(i+1,min(nchi,i+a+1),dtype=int)
            if not len(j):continue
            b=a+i-j
            w=invariant_kernel(M,G,float(e),float(ka),E[j],k[b],
                               gchi=gchi,ge=ge,nq=nq)*dE**3
            nonzero=w>0
            i_all.extend([i]*int(nonzero.sum()))
            a_all.extend([a]*int(nonzero.sum()))
            j_all.extend(j[nonzero].tolist())
            b_all.extend(b[nonzero].tolist())
            w_all.extend(w[nonzero].tolist())
    idx=np.array(i_all,dtype=int);a=np.array(a_all,dtype=int)
    j=np.array(j_all,dtype=int);b=np.array(b_all,dtype=int)
    W=np.array(w_all,dtype=float)
    if not np.allclose(E[idx]+k[a],E[j]+k[b],atol=1e-12,rtol=1e-13):
        raise ArithmeticError("discrete four-particle energy mismatch")
    return {"Echi":E,"Ee":k,"Vchi":Vchi,"Ve":Ve,
            "i":idx,"a":a,"j":j,"b":b,"W":W,
            "params":{"M":M,"T":T,"G":G,"dE":dE,"nchi":nchi,"ne":ne,
                      "gchi":gchi,"ge":ge,"nq":nq},
            "rates_are_direct_microreversible_quadrature":True,
            "not_exact_continuum_collision":True}


def sources(g,fchi=None,fe=None):
    E=g["Echi"];k=g["Ee"]
    if fchi is None:fchi=expit(-E/g["params"]["T"])
    if fe is None:fe=expit(-k/g["params"]["T"])
    f=np.asarray(fchi,float);v=np.asarray(fe,float)
    if (f.shape!=E.shape or v.shape!=k.shape
        or np.any(~np.isfinite(f)) or np.any(~np.isfinite(v))
        or np.any((f<0)|(f>1)) or np.any((v<0)|(v>1))):
        raise ValueError("valid physical Fermi occupations required")
    i=g["i"];a=g["a"];j=g["j"];b=g["b"];W=g["W"]
    F=f[i]*v[a]*(1-f[j])*(1-v[b])
    R=f[j]*v[b]*(1-f[i])*(1-v[a])
    flow=W*(F-R)
    dchi=np.zeros_like(f);de=np.zeros_like(v)
    np.add.at(dchi,i,-flow);np.add.at(dchi,j,flow)
    np.add.at(de,a,-flow);np.add.at(de,b,flow)
    entropy=float(np.dot(flow,np.log(F/R))) if np.all(F>0)&np.all(R>0) else None
    Qchi=float(np.dot(E,dchi));Qe=float(np.dot(k,de))
    return {"dchi":dchi,"de":de,"Qchi_GeV5":Qchi,"Qe_GeV5":Qe,
            "Nchi_source_GeV4":float(np.sum(dchi)),
            "Ne_source_GeV4":float(np.sum(de)),
            "energy_residual_GeV5":Qchi+Qe,
            "entropy_GeV4":entropy,
            "no_posthoc_rate_symmetrisation":True}


def _entropy(f):
    """Safe elementwise 0<=f<=1 Fermi entropy, continuous endpoints."""
    f=np.asarray(f,float)
    result=np.zeros_like(f)
    mid=(f>0)&(f<1)
    v=f[mid]
    result[mid]=-v*np.log(v)-(1-v)*np.log1p(-v)
    return result


def invariants(g,f,v):
    E=g["Echi"];k=g["Ee"];Vc=g["Vchi"];Ve=g["Ve"]
    return (float(np.dot(Vc,f)),float(np.dot(Ve,v)),
            float(np.dot(Vc*f,E)+np.dot(Ve*v,k)),
            float(np.dot(Vc,_entropy(f))+np.dot(Ve,_entropy(v))))


def evolve(g,fchi,fe,dt,substeps=1,reverse=False):
    """Backward-Euler Lie splitting, every quartet physical bare rate W.

    No posthoc sym: microscopic rate is directly computed q integral.
    Discrete charge/energy/entropy/Pauli exact but time-splitting error
    must be measured, and this is ONLY an energy-collocation model.
    """
    f=np.asarray(fchi,float).copy();v=np.asarray(fe,float).copy()
    if (f.shape!=g["Echi"].shape or v.shape!=g["Ee"].shape or
        np.any(~np.isfinite(f)) or np.any(~np.isfinite(v)) or
        np.any((f<0)|(f>1)) or np.any((v<0)|(v>1))
        or type(substeps) is not int or substeps<1
        or not math.isfinite(dt) or dt<0):
        raise ValueError("valid Fermi occupations, dt and substeps required")
    baseline=invariants(g,f,v)
    idx=range(len(g["W"])-1,-1,-1) if reverse else range(len(g["W"]))
    for _ in range(substeps):
        for ix in idx:
            i,a,j,b=(int(g[name][ix]) for name in ("i","a","j","b"))
            V1,V2,V3,V4=(g["Vchi"][i],g["Ve"][a],g["Vchi"][j],g["Ve"][b])
            fi,ga,fj,gb=f[i],v[a],f[j],v[b]
            W=float(g["W"][ix])
            lo=-min(V3*fj,V4*gb,V1*(1-fi),V2*(1-ga))
            hi=min(V1*fi,V2*ga,V3*(1-fj),V4*(1-gb))
            def reaction(z):
                x,y,u,h=fi-z/V1,ga-z/V2,fj+z/V3,gb+z/V4
                return z-dt/substeps*W*(x*y*(1-u)*(1-h)-u*h*(1-x)*(1-y))
            if hi>lo and W>0 and dt>0:
                transfer=float(brentq(reaction,lo,hi,xtol=5e-15,rtol=5e-15))
                f[i]=fi-transfer/V1;v[a]=ga-transfer/V2
                f[j]=fj+transfer/V3;v[b]=gb+transfer/V4
    final=invariants(g,f,v)
    return {"fchi":f,"fe":v,"Nchi_residual":final[0]-baseline[0],
            "Ne_residual":final[1]-baseline[1],
            "energy_residual_GeV4":final[2]-baseline[2],
            "entropy_gain_GeV3":final[3]-baseline[3],
            "Qchi_integrated_GeV4":float(np.dot(g["Vchi"]*f,g["Echi"])-
                                         np.dot(g["Vchi"]*fchi,g["Echi"])),
            "edge_count":len(g["W"]),
            "not_exact_continuum_collision":True}
