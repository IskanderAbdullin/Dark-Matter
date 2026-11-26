"""Full relativistic Pauli chi-e source for arbitrary isotropic fchi(p).

Momentum-resolved Cchi(p) and its weak number, energy and entropy moments
are integrated over thermal massless electrons and full outgoing CM angles
with the physical tree vector-contact |M|^2, v_Moller and final-state Pauli.
Arbitrary supplied fchi(p) is a callable; a PCHIP tabulation is available.
Charge and energy conservation are verified against the independent
symmetric weak-form rather than imposed by numerical projection.

A fixed electron FD bath is used as a reference, with Qe=-Qchi ledger;
it is NOT a coupled evolving SM sector, BGK rate, solved full Boltzmann
evolution, UV-completed SM model or 86-eV DM prediction.
"""
from __future__ import annotations
import math
import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.special import expit
from relativistic_weak_form_2to2 import cm_event,_gauss_span,_input_checks


def distribution_from_grid(q,f):
    """Monotone-piecewise cubic shape-preserving interpolator, 0 above qmax.

    The explicit finite support can create artificial entropy/heat even if
    the supplied samples were exact FD values; always scan qmax and nq.
    """
    q=np.asarray(q,float);f=np.asarray(f,float)
    if (q.ndim!=1 or f.shape!=q.shape or len(q)<3 or
            np.any(~np.isfinite(q)) or np.any(~np.isfinite(f)) or
            np.any(np.diff(q)<=0) or q[0]!=0 or
            np.any((f<0)|(f>1))):
        raise ValueError("q must start at zero, grow strictly, 0<=f<=1")
    interpolator=PchipInterpolator(q,f,extrapolate=False)
    def occupation(momentum):
        momentum=np.asarray(momentum,dtype=float)
        vals=interpolator(np.clip(momentum,0,q[-1]))
        return np.where(momentum>q[-1],0.,vals)
    return occupation


def collision_from_nonthermal_chi(chi_distribution,*,M=.1,Te=1.,G=.001,
                                  mue=0.,gchi=4.,ge=4.,
                                  n_p=12,n_k=12,n_mu=10,n_cm=12,n_phi=12,
                                  cutoff=24.,p_max=None):
    """Nonlinear 2->2 strong Cchi(p), weak Q, Qe and entropy at one instant.

    p_max sets physical incoming chi quadrature boundary separately from
    electron k_max=cutoff*Te; mode Cchi(p) uses full reverse Pauli term.
    For strictly zero occupations or truncated tails the entropy diagnostic
    regularises log(F/R) at 1e-290; do not promote it to an exact H theorem.
    """
    _input_checks(M,Te,Te,G,gchi,ge,n_p,n_k,n_mu,n_cm,n_phi,cutoff)
    if not math.isfinite(mue):
        raise ValueError("mue must be finite")
    if p_max is None:
        p_max=cutoff*Te
    if not math.isfinite(p_max) or p_max<=0:
        raise ValueError("valid p_max required")
    pgrid,pw=_gauss_span(n_p,0.,p_max)
    kgrid,kw=_gauss_span(n_k,0.,cutoff*Te)
    mugrid,mw=_gauss_span(n_mu,-1.,1.)
    ugrid,uw=_gauss_span(n_cm,-1.,1.)
    phigrid=2*np.pi*(np.arange(n_phi)+.5)/n_phi
    u=ugrid[:,None];phi=phigrid[None,:]
    one_w=uw[:,None]/n_phi
    strong_number=strong_energy=0.
    heat_weak=entropy_weak=0.
    strong_modes=[]
    max_negative_entropy=0.
    pref=gchi*ge/(8*math.pi**3)
    for p,pwgt in zip(pgrid,pw):
        E=math.hypot(M,p)
        fi=float(np.asarray(chi_distribution(p)))
        if not math.isfinite(fi) or not 0<=fi<=1:
            raise ValueError("invalid incoming chi occupation")
        pinc=0.
        for k,kwgt in zip(kgrid,kw):
            gi=float(expit((mue-k)/Te))
            for mu,muw in zip(mugrid,mw):
                Eout,keout,t,_=cm_event(M,p,k,mu,u,phi)
                pout=np.sqrt(np.maximum(0.,Eout*Eout-M*M))
                f1=np.asarray(chi_distribution(pout))
                if (f1.shape!=pout.shape or np.any(~np.isfinite(f1)) or
                        np.any((f1<0)|(f1>1))):
                    raise ValueError("invalid final-state chi occupation")
                g1=expit((mue-keout)/Te)
                F=fi*gi*(1-f1)*(1-g1)
                R=f1*g1*(1-fi)*(1-gi)
                net=F-R
                log_affinity=np.log(np.maximum(F,1e-290)/np.maximum(R,1e-290))
                local_entropy=net*log_affinity
                max_negative_entropy=min(max_negative_entropy,
                                         float(np.min(local_entropy)))
                s=M*M+2*k*(E-p*mu)
                D=s-M*M
                us=2*M*M-s-t
                amp=2*G*G*(D*D+(us-M*M)**2+2*M*M*t)
                dsdom=amp/(64*math.pi**2*s)
                measure=kwgt*k*k*muw*(1-p*mu/E)
                heat_weak+=pwgt*p*p*measure*float(
                    np.sum(dsdom*net*(Eout-E)*one_w))
                entropy_weak+=pwgt*p*p*measure*float(
                    np.sum(dsdom*local_entropy*one_w))
                pinc+=measure*float(np.sum(dsdom*(-net)*one_w))
        C=ge/(2*math.pi)*pinc
        strong_modes.append((float(p),float(C)))
        strong_number+=gchi/(2*math.pi**2)*pwgt*p*p*C
        strong_energy+=gchi/(2*math.pi**2)*pwgt*p*p*E*C
    Q=pref*heat_weak
    return {
      "Q_chi_GeV5":float(Q),"Q_e_GeV5":float(-Q),
      "total_entropy_production_GeV4":float(pref*entropy_weak),
      "negative_event_entropy_floor":float(max_negative_entropy),
      "number_strong_GeV4":float(strong_number),
      "energy_strong_GeV5":float(strong_energy),
      "strong_vs_weak_heat_relative":float((strong_energy-Q)/Q) if Q else float(strong_energy),
      "strong_modes":strong_modes,
      "quadrature_points":n_p*n_k*n_mu*n_cm*n_phi,
      "scope":"generic isotropic chi occupation at one instant; finite electron FD bath; no time evolution, SM feedback or BEST; entropy log regularised at 1e-290",
      "no_gamma_el_or_relic":True
    }
