"""5D deterministic Lorentz-invariant Pauli-blocked chi+e 2-to-2 weak source.

d rho_chi / dt = (gchi*ge/2) * integral[d^3p d^3k/(2pi)^6]
    v_M * integral dOmega_cm (d sigma/dOmega_cm) (F-R) (Echi'-Echi).
F=fchi fe (1-fchi') (1-fe'), R=fchi' fe' (1-fchi)(1-fe).
The 1/2 is ONLY for the symmetrised forward-minus-reverse weak form.
Incoming spins are averaged in |M|^2; incoming gchi and ge multiply once.
The electron is massless, and fchi and fe are isotropic FD trial shapes.

A separate unsymmetrised Cchi(p) is integrated over electron and outgoing
angles. Its number moment and energy moment must CONVERGE to their respective
weak-form invariant/energy limit, not be imposed by projection; finite
p/k cutoffs can produce real residuals. Results are in natural GeV units.

NOT a generic nonthermal evolving collision solver, not BEST gamma_el and
not a 86-eV DM relic.  The electron-vector EFT is broken-EW low-energy
only and lacks neutrino, electron mass, medium effects, UV completion.
"""
from __future__ import annotations
import math
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import expit


def _gauss_span(n, low, high):
    xx, ww = leggauss(n)
    return (low + (high-low)*(1+xx)/2), (high-low)*ww/2


def _input_checks(M,Tchi,Te,G,gchi,ge,n_p,n_k,n_mu,n_cm,n_phi,cutoff):
    if not all(math.isfinite(float(a)) for a in (M,Tchi,Te,G,gchi,ge,cutoff)):
        raise ValueError("nonfinite physical parameter")
    if M<=0 or Tchi<=0 or Te<=0 or G<0 or gchi<=0 or ge<=0 or cutoff<=2:
        raise ValueError("positive masses/temperatures/degeneracies, G>=0, cutoff>2")
    for n in (n_p,n_k,n_mu,n_cm,n_phi):
        if type(n) is not int or n<3:
            raise ValueError("mesh n >= 3 integer required")


def cm_event(M,p,k,mu,cos_star,phi):
    """Exact incoming lab to outgoing CM to lab scattering kinematics.

    Incoming chi=(E,0,0,p), e=(k,k*sqrt(1-mu^2),0,k*mu).
    cos_star,phi broadcast as numpy arrays. Returns Echi',Ee',t,shell.
    """
    E=math.hypot(M,p)
    rx=k*math.sqrt(max(0.,1-mu*mu)); rz=p+k*mu
    Etot=E+k
    s=M*M+2*k*(E-p*mu)
    root=math.sqrt(s)
    pstar=(s-M*M)/(2*root)
    Estar=(s+M*M)/(2*root)
    # Inverse Lorentz transform without dividing by small beta^2.
    factor=p*rz/(root*(Etot+root))-E/root
    cx=factor*rx; cz=p+factor*rz
    clen=math.hypot(cx,cz)
    if clen<=0:
        raise ArithmeticError("degenerate incoming CM direction")
    nx=cx/clen; nz=cz/clen
    v=np.asarray(cos_star);w=np.asarray(phi)
    st=np.sqrt(np.maximum(0.,1-v*v))
    xp=pstar*(v*nx-st*np.cos(w)*nz)
    yp=pstar*st*np.sin(w)
    zp=pstar*(v*nz+st*np.cos(w)*nx)
    dot=rx*xp+rz*zp
    Eprime=(Etot*Estar+dot)/root
    Eeprime=Etot-Eprime
    fac=dot/(root*(Etot+root))+Estar/root
    pxprime=xp+fac*rx
    pzprime=zp+fac*rz
    epx=rx-pxprime;epy=-yp;epz=rz-pzprime
    shell=Eprime*Eprime-(pxprime*pxprime+yp*yp+pzprime*pzprime)-M*M
    e_shell=Eeprime*Eeprime-(epx*epx+epy*epy+epz*epz)
    t=-2*pstar*pstar*(1-v)
    defect=np.maximum(np.abs(shell),np.abs(e_shell))/max(1,s)
    return Eprime,Eeprime,t,defect


def weak_energy_moments(*, M=.1,Tchi=.8,Te=1.,G=.001,
                        muchi=0.,mue=0.,gchi=4.,ge=4.,
                        n_p=10,n_k=10,n_mu=8,n_cm=10,n_phi=10,
                        cutoff=18.):
    """Pauli-inclusive isotropic thermal collision moments; Q in GeV^5.

    For nonzero chemical potentials with aggregated particle/antiparticle
    degeneracy, these are mathematical equal-mu trial populations, not
    a charge-conjugate early-universe chemical potential prescription.
    Physically charge-symmetric benchmark uses muchi=mue=0.
    """
    _input_checks(M,Tchi,Te,G,gchi,ge,n_p,n_k,n_mu,n_cm,n_phi,cutoff)
    if not all(math.isfinite(t) for t in (muchi,mue)):
        raise ValueError("chemical potentials must be finite")
    if G==0:
        return {"Q_chi_GeV5":0.,"Q_e_GeV5":0.,"entropy_GeV4":0.,
                "Q_fraction_dimensionless":0.,"max_mass_shell_defect":0.,
                "number_strong_GeV4":0.,"energy_strong_GeV5":0.,
                "weak_vs_strong_heat_relative_discrepancy":0.,
                "equilibrium_control_direct_difference":0.,
                "quadrature_points":0,"strong_modes":[],"no_gamma_el_or_relic":True}
    pgrid,pw=_gauss_span(n_p,0.,cutoff*Tchi)
    kgrid,kw=_gauss_span(n_k,0.,cutoff*Te)
    mugrid,mw=_gauss_span(n_mu,-1.,1.)
    ugrid,uw=_gauss_span(n_cm,-1.,1.)
    phigrid=2*np.pi*(np.arange(n_phi)+.5)/n_phi
    u=ugrid[:,None];phi=phigrid[None,:]
    one_w=uw[:,None]/n_phi
    total=0.;unpaired=0.;mass_shell=0.
    strong_number=0.;strong_energy=0.;p_rates=[]
    # After incident angular integrations: gchi*ge/(8*pi^4).
    # Outgoing azimuth integral 2*pi and symmetrising factor 1/2.
    pref=gchi*ge/(8*math.pi**3)
    for p,p_weight in zip(pgrid,pw):
        E=math.hypot(M,p)
        p_integral=0.
        fi=float(expit((muchi-E)/Tchi))
        for k,k_weight in zip(kgrid,kw):
            gi=float(expit((mue-k)/Te))
            for mu,mu_weight in zip(mugrid,mw):
                s=M*M+2*k*(E-p*mu)
                vflux=1-p*mu/E
                Eprime,Eeprime,t,defect=cm_event(M,p,k,mu,u,phi)
                D=s-M*M
                us=2*M*M-s-t
                # Same averaged vector-contact |M|^2 as microphysics_portal.
                amp=2*G*G*(D*D+(us-M*M)**2+2*M*M*t)
                dsdom=amp/(64*math.pi**2*s)
                f1=expit((muchi-Eprime)/Tchi)
                f2=expit((mue-Eeprime)/Te)
                F=fi*gi*(1-f1)*(1-f2)
                R=f1*f2*(1-fi)*(1-gi)
                A=(Eprime-E)*(1/Tchi-1/Te)
                # Equivalent to F-R; stable exact common-T cancellation.
                net=F*(-np.expm1(-A))
                if Tchi==Te:
                    unpaired=max(unpaired,float(np.max(np.abs(F-R))))
                delta=Eprime-E
                integral=float(np.sum(dsdom*net*delta*one_w))
                measure=k_weight*k*k*mu_weight*vflux
                total+=p_weight*p*p*measure*integral
                # Unprojected mode-resolved strong collision: R-F.
                # Finite p/k cutoffs can cause charge and energy leakage.
                p_integral+=measure*float(np.sum(dsdom*(-net)*one_w))
                mass_shell=max(mass_shell,float(np.max(defect)))
        C_at_p=ge/(2*math.pi)*p_integral
        p_rates.append((float(p),float(C_at_p)))
        strong_number+=gchi/(2*math.pi**2)*p_weight*p*p*C_at_p
        strong_energy+=gchi/(2*math.pi**2)*p_weight*p*p*E*C_at_p
    Q=pref*total
    return {
        "Q_chi_GeV5":float(Q),
        "number_strong_GeV4":float(strong_number),
        "energy_strong_GeV5":float(strong_energy),
        "weak_vs_strong_heat_relative_discrepancy":float((strong_energy-Q)/Q) if Q else float(strong_energy),
        "strong_modes":p_rates,
        "Q_e_GeV5":float(-Q),
        "entropy_GeV4":float(Q*(1/Tchi-1/Te)),
        "Q_fraction_dimensionless":float(Q/Te**5),
        "max_mass_shell_defect":mass_shell,
        "equilibrium_control_direct_difference":unpaired,
        "quadrature_points":n_p*n_k*n_mu*n_cm*n_phi,
        "parameters":dict(M=M,Tchi=Tchi,Te=Te,G=G,muchi=muchi,mue=mue,
                          gchi=gchi,ge=ge,cutoff=cutoff,n_p=n_p,n_k=n_k,
                          n_mu=n_mu,n_cm=n_cm,n_phi=n_phi),
        "physical_scope":"Electron-only vector EFT, FD isotropic trial shapes, integrated weak moment and sampled strong Cchi(p)",
        "no_gamma_el_or_relic":True,
    }


def rest_blocked_heat_independent(*,M=1.,Te=.1,G=.2,ge=4.,n_k=24,n_cm=24,cutoff=30.):
    """Independent 2D Gauss angular+k check vs existing t/quadrature rest Q.

    d sigma/d cos(theta*) = |M|^2/(32*pi*s). This is NOT gas Q.
    """
    if M<=0 or Te<=0 or G<0 or ge<=0:
        raise ValueError("physical inputs required")
    ks,ws=_gauss_span(n_k,0,cutoff*Te)
    us,uw=_gauss_span(n_cm,-1,1)
    total=0.
    for k,w in zip(ks,ws):
        s=M*M+2*M*k
        ps=(s-M*M)/(2*math.sqrt(s))
        t=-2*ps*ps*(1-us)
        D=s-M*M
        v=2*M*M-s-t
        amp=2*G*G*(D*D+(v-M*M)**2+2*M*M*t)
        dsdu=amp/(32*math.pi*s)
        delta=-t/(2*M)
        fe_prime=expit(-(k-delta)/Te)
        integral=float(np.dot(uw,dsdu*delta*(1-fe_prime)))
        total+=w*k*k*float(expit(-k/Te))*integral
    return ge/(2*math.pi**2)*total
