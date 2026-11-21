"""Physical 5D chi-e scattering event -> energy-conserving fermionic quartets.

The chi and electron energy bins share spacing dE:
 Echi_i=M+(i+1/2)dE, Ee_a=(a+1/2)dE.
Every deposited event i,a -> j,b has b=i+a-j, so exact BIN CENTER
total energy conservation.  The actual continuous CM scattering uses
the exact vector-contact |M|^2 and full Pauli factors; its projected
positive *bare* Wraw is NOT automatically reversible.  Wpair=(Wraw +
Wraw_reversed)/2 is an explicitly DECLARED modification at finite
resolution, not the exact physical nonlinear continuum collision.
Fermion chi/e separately conserve particle number. A one-quartet
implicit Pauli/energy/entropy-safe kinetic step is provided, but no
global self-consistent cosmology, Bose relic or validated gamma_el.
"""
from __future__ import annotations
import math
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import expit
from microphysical_equilibrium_graph import cm_event


def span(n,a,b):
    x,w=leggauss(n)
    return (a+b)/2+(b-a)*x/2,(b-a)*w/2


def quartets(*,M=.1,T=1.,G=.001,gchi=4.,ge=4.,dE=.5,nchi=30,ne=30,
             nmu=8,ncm=8,nphi=8):
    if (not all(math.isfinite(float(z)) for z in (M,T,G,gchi,ge,dE))
        or M<=0 or T<=0 or G<0 or gchi<=0 or ge<=0 or dE<=0
        or any(type(z) is not int or z<3 for z in (nchi,ne,nmu,ncm,nphi))
        or nchi*ne>3000):
        raise ValueError("physical finite parameters and modest grid required")
    Echi=M+(np.arange(nchi)+.5)*dE
    Ee=(np.arange(ne)+.5)*dE
    cedges=M+dE*np.arange(nchi+1)
    pedges=np.sqrt(np.maximum(0.,cedges**2-M*M))
    Vchi=gchi/(6*math.pi**2)*np.diff(pedges**3)
    eedges=dE*np.arange(ne+1)
    Ve=ge/(6*math.pi**2)*np.diff(eedges**3)
    mus,mw=span(nmu,-1,1)
    us,uw=span(ncm,-1,1)
    phis=2*np.pi*(np.arange(nphi)+.5)/nphi
    U=us[:,None];Phi=phis[None,:]
    aw=uw[:,None]*np.ones((1,nphi))/nphi
    W=np.zeros((nchi,ne,nchi,ne))
    lost=total=lost_eq=total_eq=0.
    Dcontinuous=0.
    for i,E in enumerate(Echi):
        p=math.sqrt(E*E-M*M)
        for a,k in enumerate(Ee):
            for mu,wmu in zip(mus,mw):
                s=M*M+2*k*(E-p*mu)
                vflux=1-p*mu/E
                Ep,kp,t=cm_event(M,p,k,mu,U,Phi)
                delta=s-M*M;u=2*M*M-s-t
                amp=2*G*G*(delta*delta+(u-M*M)**2+2*M*M*t)
                dsdom=amp/(64*math.pi**2*s)
                # Incoming relative direction 1/2 dmu; outgoing 2pi du.
                # Vi and Ve have the full g/(2pi)^3 phase-space capacities.
                weight=Vchi[i]*Ve[a]*math.pi*wmu*vflux*dsdom*aw
                ep=Ep.ravel();w=weight.ravel()
                F=expit(-E/T)*expit(-k/T)*(1-expit(-ep/T))*(1-expit(-kp.ravel()/T))
                total+=float(np.sum(w));total_eq+=float(np.sum(w*F))
                Dcontinuous+=.5*float(np.dot(w*F,((ep-E)/T)**2))
                z=(ep-M)/dE-.5
                j=np.floor(z).astype(int)
                alpha=z-j
                for dest,fraction in ((j,1-alpha),(j+1,alpha)):
                    b=i+a-dest
                    valid=(dest>=0)&(dest<nchi)&(b>=0)&(b<ne)&(fraction>=0)&(fraction<=1)
                    bad=np.maximum(0.,fraction[~valid])
                    lost+=float(np.sum(w[~valid]*bad))
                    lost_eq+=float(np.sum(w[~valid]*bad*F[~valid]))
                    nn=int(np.sum(valid))
                    np.add.at(W,(np.full(nn,i),np.full(nn,a),
                                 dest[valid],b[valid]),w[valid]*fraction[valid])
    reversed_W=W.transpose(2,3,0,1)
    paired=(W+reversed_W)/2
    fc=expit(-Echi/T);fe=expit(-Ee/T)
    Feq=fc[:,None,None,None]*fe[None,:,None,None]*(1-fc)[None,None,:,None]*(1-fe)[None,None,None,:]
    d=Echi[None,None,:,None]-Echi[:,None,None,None]
    Dpair=.5*float(np.sum(paired*Feq*(d/T)**2))
    weight_norm=float(np.sum((W+reversed_W)*Feq))
    asym=float(np.sum(np.abs((W-reversed_W)*Feq)))/weight_norm if weight_norm>0 else 0.
    raw_asym=float(np.sum(np.abs(W-reversed_W)))/(2*np.sum(W)) if np.any(W) else 0.
    return dict(Echi=Echi,Ee=Ee,Vchi=Vchi,Ve=Ve,Wraw=W,W=paired,
        W_asym_L1=raw_asym,W_eq_Pauli_asym_L1=asym,
        D_quartet_discrete=Dpair,D_continuous_5D=Dcontinuous,
        omitted_raw_weight_fraction=lost/total if total else 0.,
        omitted_FD_weight_fraction=lost_eq/total_eq if total_eq else 0.,
        model=dict(M=M,T=T,G=G,gchi=gchi,ge=ge,dE=dE,nchi=nchi,ne=ne,
                   nmu=nmu,ncm=ncm,nphi=nphi),
        exact_unmodified_physical_collision=False)


def sources(g,fchi=None,fe=None,raw=False):
    """Both fermion Pauli directions, χ/e number/energy/entropy ledgers.

    With raw=True, show false thermal source of unpaired physical projection.
    """
    E=g["Echi"];K=g["Ee"];T=g["model"]["T"]
    W=g["W"]
    if fchi is None:fchi=expit(-E/T)
    if fe is None:fe=expit(-K/T)
    x=np.asarray(fchi);y=np.asarray(fe)
    if (x.shape!=E.shape or y.shape!=K.shape
        or np.any(~np.isfinite(x)) or np.any(~np.isfinite(y))
        or np.any((x<0)|(x>1)) or np.any((y<0)|(y>1))):
        raise ValueError("fermionic χ/e occupations required")
    F=x[:,None,None,None]*y[None,:,None,None]*(1-x)[None,None,:,None]*(1-y)[None,None,None,:]
    R=x[None,None,:,None]*y[None,None,None,:]*(1-x)[:,None,None,None]*(1-y)[None,:,None,None]
    if raw:
        L=g["Wraw"];J=.5*(L*F-L.transpose(2,3,0,1)*R)
    else:
        J=.5*W*(F-R)
    dchi=np.sum(J,axis=(0,1,3))-np.sum(J,axis=(1,2,3))
    de=np.sum(J,axis=(0,1,2))-np.sum(J,axis=(0,2,3))
    Qchi=float(np.dot(E,dchi));Qe=float(np.dot(K,de))
    mask=(F>0)&(R>0)&(W>0)
    full_entropy=(not raw and np.all((x>0)&(x<1)) and np.all((y>0)&(y<1)))
    entropy=float(np.sum(J[mask]*np.log(F[mask]/R[mask]))) if full_entropy else None
    return dict(dchi=dchi,de=de,Qchi=Qchi,Qe=Qe,
       number_chi=float(np.sum(dchi)),number_e=float(np.sum(de)),
       energy_total=Qchi+Qe,entropy=entropy,
       max_event_balance=float(np.max(np.abs(F-R))) if np.any(W) else 0.)


def implicit_one_quartet(g,i,a,j,b,fi,ga,fj,gb,dt):
    """Bounded implicit four-state full-Pauli step with PAIR-AVERAGED 5D rate.

    Each physical discrete quartet independently conserves χ/e number and
    total χ+e energy; entropy is nondecreasing, no clipping of occupations.
    The pair-averaged W differs from raw microphysical finite-bin rates.
    """
    from scipy.optimize import brentq
    from fermion_finite_bath_pairs import _fermi_entropy
    E=g["Echi"];K=g["Ee"];Vc=g["Vchi"];Ve=g["Ve"]
    if (not 0<=i<len(E) or not 0<=j<len(E) or not 0<=a<len(K)
        or not 0<=b<len(K) or (i,a)==(j,b)):
        raise ValueError("different in-bounds quartet required")
    if abs(E[i]+K[a]-E[j]-K[b])>1e-12:
        raise ValueError("quartet violates center-energy conservation")
    vals=(fi,ga,fj,gb,dt)
    if not all(math.isfinite(v) for v in vals) or not all(0<=v<=1 for v in vals[:4]) or dt<0:
        raise ValueError("bounded occupations and nonnegative dt required")
    w=float(g["W"][i,a,j,b])
    capacities=(Vc[i],Ve[a],Vc[j],Ve[b])
    S0=sum(v*_fermi_entropy(f) for v,f in zip(capacities,vals[:4]))
    lo=-min(Vc[j]*fj,Ve[b]*gb,Vc[i]*(1-fi),Ve[a]*(1-ga))
    hi=min(Vc[i]*fi,Ve[a]*ga,Vc[j]*(1-fj),Ve[b]*(1-gb))
    def state(z):
        return fi-z/Vc[i],ga-z/Ve[a],fj+z/Vc[j],gb+z/Ve[b]
    def J(z):
        x,y,u,v=state(z)
        return w*(x*y*(1-u)*(1-v)-u*v*(1-x)*(1-y))
    amount=0. if w==0 or dt==0 else float(brentq(
        lambda z:z-dt*J(z),lo,hi,xtol=5e-15,rtol=5e-15))
    x,y,u,v=state(amount)
    S1=sum(c*_fermi_entropy(f) for c,f in zip(capacities,(x,y,u,v)))
    energy_chi=float((E[j]-E[i])*amount)
    energy_e=float((K[b]-K[a])*amount)
    return dict(fi=x,ga=y,fj=u,gb=v,transfer_number_GeV3=amount,
        Q_dark_integrated_GeV4=energy_chi,Q_e_integrated_GeV4=energy_e,
        chi_charge_residual=float(Vc[i]*x+Vc[j]*u-Vc[i]*fi-Vc[j]*fj),
        electron_charge_residual=float(Ve[a]*y+Ve[b]*v-Ve[a]*ga-Ve[b]*gb),
        energy_residual_GeV4=energy_chi+energy_e,
        entropy_change_GeV3=float(S1-S0),
        pair_averaged_5D_W_GeV4=w,not_exact_unprojected_microphysics=True)
