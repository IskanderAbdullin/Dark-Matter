"""Full five-dimensional thermal e−e+ ↔ gamma gamma on-shell QED collision moments.

Independent from prior total-loss σv implementation: p,k,relative μ,CM photon
cos θ*, azimuth φ*; explicit Lorentz boosts of BOTH photons, physical QED
tree-level spin-averaged differential amplitude, initial Fermi and final
Bose for forward PLUS inverse. Thermal FD/BE subfamily ONLY; no general
nonthermal C[f], no time evolution, neutrino or dark coupling. Units GeV.
"""
from __future__ import annotations
import math
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import expit

ALPHA=1/137.035999084


def gauss(a,b,n):
    x,w=leggauss(n)
    return a+(b-a)*(x+1)/2,w*(b-a)/2


def thermal_collision(Te,Tgamma,mu_minus=0.,mu_plus=0.,*,m=1.,
                      alpha=ALPHA,nrad=24,nang=20,nphi=12,pcut_over_T=32):
    """Physical 5D homogeneous thermal QED reaction and energy moments.

    Tgamma=None implements *vacuum final-state LOSS-ONLY control*;
    with positive Tgamma, exact Fermi/Bose forward+inverse and entropy.
    e−,e+ each have two spins and 2γ are identical; photon polarizations
    are already summed into H(s,theta) with photon 1/2 in differential σ.
    """
    if not (math.isfinite(Te) and Te>0 and
            (Tgamma is None or (math.isfinite(Tgamma) and Tgamma>0)) and
            m>0 and alpha>0 and nrad>=8 and nang>=8 and nphi>=4 and
            pcut_over_T>0 and all(math.isfinite(x) for x in (mu_minus,mu_plus))):
        raise ValueError("positive thermal scales, physical chemical potentials")
    q,wq=gauss(0,pcut_over_T,nrad)
    mu,wm=gauss(-1,1,nang)
    c,wc=gauss(-1,1,nang)
    phi,wphi=gauss(0,2*math.pi,nphi)
    p=Te*q
    E=np.hypot(p,m)
    fm=expit((mu_minus-E)/Te)
    fp=expit((mu_plus-E)/Te)
    c3=c[None,:,None]
    sin3=np.sqrt(np.maximum(0,1-c3*c3))
    cosphi=np.cos(phi)[None,None,:]
    angle_weights=wm[:,None,None]*wc[None,:,None]*wphi[None,None,:]/(2*math.pi)
    tot_forward=tot_inverse=tot_net=tot_entropy=0.
    q_minus=q_plus=q_gamma=0.
    max_shell_error=0.
    for i,(pi,Ei) in enumerate(zip(p,E)):
      for j,(kj,wj) in enumerate(zip(p,E)):
        Et=Ei+wj
        P2=pi*pi+kj*kj+2*pi*kj*mu
        P=np.sqrt(np.maximum(0.,P2))
        s=Et*Et-P2
        beta2=np.maximum(0.,1-4*m*m/s)
        boost=P/Et
        gamma=Et/np.sqrt(s)
        pcm=.5*np.sqrt(np.maximum(0.,s-4*m*m))
        parallel=np.where(P>0,pi*(pi+kj*mu)/np.maximum(P,1e-100),0.)
        cosa=np.clip(gamma*(parallel-boost*Ei)/np.maximum(pcm,1e-100),-1.,1.)
        cosout=np.clip(cosa[:,None,None]*c3+
                       np.sqrt(np.maximum(0.,1-cosa[:,None,None]**2))*
                       sin3*cosphi,-1,1)
        photon1=Et/2*(1+boost[:,None,None]*cosout)
        photon2=Et-photon1
        max_shell_error=max(max_shell_error,
                            float(np.max(abs(photon1+photon2-Et))))
        D=1-beta2[:,None,None]*c3*c3
        H=(2-D)/D+2*(1-beta2[:,None,None])/D-2*(1-beta2[:,None,None])**2/D**2
        if np.any(H<0):raise ArithmeticError("negative QED spin amplitude")
        if Tgamma is not None:
            b1=1/np.expm1(photon1/Tgamma)
            b2=1/np.expm1(photon2/Tgamma)
            F=fm[i]*fp[j]*(1+b1)*(1+b2)
            R=(1-fm[i])*(1-fp[j])*b1*b2
            # Exact on-shell log(F/R): helps avoid subtraction at equilibrium.
            affinity=Et/Tgamma-(Et-mu_minus-mu_plus)/Te
            flow=F*(-np.expm1(-affinity))
        else:
            F=fm[i]*fp[j]
            R=np.zeros_like(F)
            flow=np.broadcast_to(F,H.shape)
            affinity=None
        # d(σ v_M,lab)/dc dφ/(2π) = πα² H/(2E−E+).
        # Incoming g−g+=2×2 and 4π×2π angular phase space give 4/(8π⁴).
        A=(Te**6*wq[i]*wq[j]*q[i]**2*q[j]**2 *
           alpha**2/(4*math.pi**3*Ei*wj))
        weight=A*angle_weights*H
        tot_forward+=float(np.sum(weight*F))
        tot_inverse+=float(np.sum(weight*R))
        dnet=float(np.sum(weight*flow))
        tot_net+=dnet
        q_minus-=Ei*dnet
        q_plus-=wj*dnet
        q_gamma+=float(np.sum(weight*(photon1+photon2)*flow))
        if Tgamma is not None:
            tot_entropy+=float(np.sum(weight*flow*affinity))
    return {
       "R_forward_GeV4":tot_forward,"R_inverse_GeV4":tot_inverse,
       "R_net_GeV4":tot_net,
       "N_e_minus_source_GeV4":-tot_net,
       "N_e_plus_source_GeV4":-tot_net,
       "N_photon_source_GeV4":2*tot_net,
       "electric_charge_source_GeV4":0.,
       "Q_e_minus_GeV5":q_minus,"Q_e_plus_GeV5":q_plus,
       "Q_photon_GeV5":q_gamma,
       "energy_residual_GeV5":q_minus+q_plus+q_gamma,
       "entropy_production_GeV4":tot_entropy if Tgamma is not None else None,
       "max_on_shell_energy_residual_GeV":max_shell_error,
       "photon_bose_and_inverse_included":Tgamma is not None,
       "complete_thermal_FD_BE_moment_not_general_collision":True,
       "nonlinear_time_evolution_not_implemented":True
    }


if __name__=="__main__":
    import argparse,json
    a=argparse.ArgumentParser()
    a.add_argument("--Te-over-me",type=float,default=.1)
    a.add_argument("--Tgamma-over-me",type=float,default=.1)
    a.add_argument("--pair-affinity-over-Te",type=float,default=.02)
    a.add_argument("--nrad",type=int,default=24)
    a.add_argument("--nang",type=int,default=20)
    o=a.parse_args()
    print(json.dumps(thermal_collision(
        o.Te_over_me,o.Tgamma_over_me,
        mu_minus=.5*o.pair_affinity_over_Te*o.Te_over_me,
        mu_plus=.5*o.pair_affinity_over_Te*o.Te_over_me,
        nrad=o.nrad,nang=o.nang),indent=2))
