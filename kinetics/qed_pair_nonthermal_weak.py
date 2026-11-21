"""Independent arbitrary isotropic nonthermal QED e−e+↔γγ five-dimensional WEAK moments.

Accepts physical fermion occupations (0<f<1) for electrons, positrons,
positive finite photon Bose-mode occupations b(E). Uses fully on-shell
CM photon boosts and tree-level spin/polarisation-averaged QED amplitude,
both e−e+ annihilation and inverse gamma gamma creation. Outputs only
INTEGRATED number, physical energy and entropy moments; not strong C[f]
on a physical phase-space mesh and NOT an evolving plasma.
"""
from __future__ import annotations
import math
import numpy as np
from numpy.polynomial.legendre import leggauss

ALPHA=1/137.035999084


def gauss(a,b,n):
    x,w=leggauss(n)
    return a+(b-a)*(x+1)/2,w*(b-a)/2


def weak_moments(fminus,fplus,bphoton,*,m=1.,Escale=.1,alpha=ALPHA,
                 nrad=24,nang=20,nphi=12,pcut_over_scale=32):
    if not (m>0 and Escale>0 and alpha>0 and nrad>=8 and
            nang>=8 and nphi>=4 and pcut_over_scale>0):
        raise ValueError("positive scales and sufficient quadrature")
    q,wq=gauss(0,pcut_over_scale,nrad)
    mu,wm=gauss(-1,1,nang)
    c,wc=gauss(-1,1,nang)
    phi,wphi=gauss(0,2*math.pi,nphi)
    p=Escale*q
    E=np.hypot(p,m)
    c3=c[None,:,None]
    sint=np.sqrt(np.maximum(0,1-c3*c3))
    cosphi=np.cos(phi)[None,None,:]
    angle_weights=wm[:,None,None]*wc[None,:,None]*wphi[None,None,:]/(2*math.pi)
    forward=inverse=net=entropy=Qminus=Qplus=Qphoton=0.
    max_shell=0.
    for i,(pi,Ei) in enumerate(zip(p,E)):
        fi=float(np.asarray(fminus(float(Ei))))
        if not (0<fi<1):
            raise ValueError("electron occupation strictly between 0 and 1")
        for j,(kj,wj) in enumerate(zip(p,E)):
            fj=float(np.asarray(fplus(float(wj))))
            if not (0<fj<1):
                raise ValueError("positron occupation strictly between 0 and 1")
            Et=Ei+wj
            P2=pi*pi+kj*kj+2*pi*kj*mu
            P=np.sqrt(np.maximum(0.,P2))
            s=Et*Et-P2
            beta2=np.maximum(0.,1-4*m*m/s)
            betaB=P/Et
            gam=Et/np.sqrt(s)
            pcm=.5*np.sqrt(np.maximum(0.,s-4*m*m))
            parallel=np.where(P>0,pi*(pi+kj*mu)/np.maximum(P,1e-100),0.)
            ci=np.clip(gam*(parallel-betaB*Ei)/np.maximum(pcm,1e-100),-1.,1.)
            cf=np.clip(ci[:,None,None]*c3+
                       np.sqrt(np.maximum(0.,1-ci[:,None,None]**2))*
                       sint*cosphi,-1.,1.)
            k1=Et/2*(1+betaB[:,None,None]*cf)
            k2=Et-k1
            max_shell=max(max_shell,float(np.max(abs(k1+k2-Et))))
            b1=np.asarray(bphoton(k1),float)
            b2=np.asarray(bphoton(k2),float)
            if (b1.shape!=k1.shape or b2.shape!=k2.shape or
                np.any(b1<=0) or np.any(b2<=0) or
                not np.all(np.isfinite(b1)) or not np.all(np.isfinite(b2))):
                raise ValueError("photon occupations finite and strictly positive")
            D=1-beta2[:,None,None]*c3*c3
            H=(2-D)/D+2*(1-beta2[:,None,None])/D-2*(1-beta2[:,None,None])**2/D**2
            if np.any(H<0):raise ArithmeticError("negative QED amplitude")
            F=fi*fj*(1+b1)*(1+b2)
            B=(1-fi)*(1-fj)*b1*b2
            flow=F-B
            # The same F,B determine both flow and log(F/B): ensures H theorem.
            logratio=np.log1p(flow/B)
            weight=(Escale**6*wq[i]*wq[j]*q[i]**2*q[j]**2 *
                    alpha**2/(4*math.pi**3*Ei*wj))*angle_weights*H
            forward+=float(np.sum(weight*F))
            inverse+=float(np.sum(weight*B))
            delta=float(np.sum(weight*flow))
            net+=delta
            entropy+=float(np.sum(weight*flow*logratio))
            Qminus-=Ei*delta
            Qplus-=wj*delta
            Qphoton+=float(np.sum(weight*(k1+k2)*flow))
    return {"R_forward_GeV4":forward,"R_inverse_GeV4":inverse,
            "R_net_GeV4":net,"N_e_minus_GeV4":-net,"N_e_plus_GeV4":-net,
            "N_photon_GeV4":2*net,"charge_residual_GeV4":0.,
            "Q_e_minus_GeV5":Qminus,"Q_e_plus_GeV5":Qplus,
            "Q_photon_GeV5":Qphoton,
            "Q_total_residual_GeV5":Qminus+Qplus+Qphoton,
            "entropy_production_GeV4":entropy,
            "max_shell_energy_error_GeV":max_shell,
            "arbitrary_nonthermal_weak_moments_only":True,
            "not_a_strong_mode_resolved_collision_or_time_integrator":True}


if __name__=="__main__":
    import argparse,json
    from scipy.special import expit
    p=argparse.ArgumentParser()
    p.add_argument("--T-over-me",type=float,default=.1)
    p.add_argument("--electron-logit-bump",type=float,default=.15)
    o=p.parse_args()
    T=o.T_over_me
    fminus=lambda E:expit(-E/T+o.electron_logit_bump*np.exp(-E/T))
    fplus=lambda E:expit(-E/T)
    bphoton=lambda E:1/np.expm1(E/T)
    print(json.dumps(weak_moments(fminus,fplus,bphoton,Escale=T),indent=2))
