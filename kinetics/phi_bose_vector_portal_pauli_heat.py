"""On-shell Pauli-blocked FULL finite-vector-mediator χ1/χ2+e± heat sources.

Real φ gauge neutral and qχ1=qχ2 under vectorlike dark U(1); both Dirac
particles and BOTH antiparticles interact through a massive vector A' and
low-energy electromagnetic lepton current εe. This is an illustrative
ELECTRON/positron-only below-EWSB kinetic-mixing portal partial channel;
muons, charged quarks/hadrons, mediator self-energy, electron mass,
neutrino mixing and kinetic backreaction are NOT yet included.

For each χi+e -> χi'+e' on-shell physical event:
    dσ/dΩcm = [2 (gD εe)^2/(mV²-t)^2
                *((s-M²)^2+(u-M²)^2+2 M² t)]/(64π² s).
True Pauli-blocked quantum microscopic F,R and exact 4-momentum event
give pairwise Qchi+Qe=0 and H=Qchi(1/Tchi-1/Te) >=0.
This is an integrated WEAK heat source for fixed FD trial spectra,
NOT a strong evolving thermal equilibration eigenvalue or FRW dark relic.
"""
import math
import json
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import expit

from kinetics.relativistic_weak_form_2to2 import cm_event
from kinetics.phi_bose_vector_portal_finite_mediator import (
    MCHI1,MCHI2,vector_portal_setup,finite_vector_cross_section,
)


def _positive_gauss(n,low,high):
    z,w=leggauss(n)
    return low+(high-low)*(z+1)/2,w*(high-low)/2


def one_dirac_charge_sector_weak_heat(
    *,mchi,Tchi,Te=1.,gD=.1,epsilon=1e-5,mV=10.,
    mu_chi=0.,mu_e=0.,gchi=2.,ge=4.,
    n_p=8,n_k=8,n_mu=8,n_cm=8,n_az=8,cutoff=20.):
    """One χ or antiχ (spin gχ=2) elastic WEAK heat vs e± (g_e=4).

    Photon and electron SM bath are here a prescribed FD temperature Te;
    only this channel Qchi=-Qe has a physical event-level heat ledger.
    """
    if min(mchi,Tchi,Te,mV,gchi,ge)<=0 or cutoff<=2:
        raise ValueError("positive physical masses/temperatures/degeneracy")
    if min(gD,epsilon)<0 or min(n_p,n_k,n_mu,n_cm,n_az)<3:
        raise ValueError("nonnegative coupling and n>=3 required")
    if not (math.isfinite(mu_chi) and math.isfinite(mu_e)):
        raise ValueError("finite FD chemical potentials required")
    ref=vector_portal_setup(gD,epsilon,mV)
    pair=ref["vector_coupling_product"]
    p,pw=_positive_gauss(n_p,0.,cutoff*Tchi)
    k,kw=_positive_gauss(n_k,0.,cutoff*Te)
    mugrid,muw=leggauss(n_mu)
    costar,cw=leggauss(n_cm)
    phi=(np.arange(n_az)+.5)*(2*math.pi/n_az)
    uc=costar[:,None]
    az=phi[None,:]
    weight_cm=cw[:,None]/n_az
    # True weak forward/reverse paired normalization (not identical spins).
    pref=gchi*ge/(8*math.pi**3)
    integral=0.
    entropy=0.
    maximum_shell=0.
    nonnegative_flux_sign=True
    n_events=0
    event_weight_total=0.
    for pa,wpa in zip(p,pw):
        E=math.hypot(pa,mchi)
        fchi=float(expit((mu_chi-E)/Tchi))
        for ke,wke in zip(k,kw):
            fe=float(expit((mu_e-ke)/Te))
            for mui,wmui in zip(mugrid,muw):
                s=mchi*mchi+2*ke*(E-pa*mui)
                vM=1-pa*mui/E
                Ep,Kp,t,shell=cm_event(mchi,pa,ke,mui,uc,az)
                mandel_u=2*mchi*mchi-s-t
                D=s-mchi*mchi
                spintrace=2*(D*D+(mandel_u-mchi*mchi)**2+2*mchi*mchi*t)
                amp2=pair*pair*spintrace/(mV*mV-t)**2
                if np.min(amp2)<-1e-11*np.max(abs(amp2)):
                    raise ArithmeticError("negative physical spin-averaged amplitude")
                dsdo=amp2/(64*math.pi**2*s)
                assert np.all(dsdo>=0)
                fchip=expit((mu_chi-Ep)/Tchi)
                fep=expit((mu_e-Kp)/Te)
                F=fchi*fe*(1-fchip)*(1-fep)
                R=fchip*fep*(1-fchi)*(1-fe)
                affinity=(Ep-E)*(1/Tchi-1/Te)
                # stable on-shell thermodynamic detailed balance and sign
                signed=np.where(affinity>=0.,
                                F*(-np.expm1(-affinity)),
                                R*np.expm1(affinity))
                if np.min(signed*affinity)<-1e-25:
                    nonnegative_flux_sign=False
                delta=Ep-E
                w=wpa*pa*pa*wke*ke*ke*wmui*vM*weight_cm*dsdo
                integral+=float(np.sum(w*signed*delta))
                entropy+=float(np.sum(w*signed*affinity))
                event_weight_total+=float(np.sum(w))
                maximum_shell=max(maximum_shell,float(np.max(shell)))
                n_events+=len(np.broadcast_to(dsdo,(n_cm,n_az)).ravel())
    Q=pref*integral
    H=pref*entropy
    if not math.isclose(H,Q*(1/Tchi-1/Te),rel_tol=2e-11,abs_tol=3e-24):
        raise ArithmeticError("wrong Pauli/entropy source")
    return dict(mchi_GeV=mchi,Tchi_GeV=Tchi,Te_GeV=Te,
                vector_mediator_GeV=mV,gD=gD,epsilon=epsilon,
                Q_chi_GeV5=float(Q),Q_e_GeV5=float(-Q),
                quantum_H_GeV4=float(H),
                max_four_shell_relative_defect=float(maximum_shell),
                microscopic_nonnegative_H=nonnegative_flux_sign,
                total_positive_event_measure=float(event_weight_total),
                number_chi_source_GeV4=0.,number_e_source_GeV4=0.,
                quadrature_orders=[n_p,n_k,n_mu,n_cm,n_az],
                cutoff_T=cutoff,n_physical_sampled_events=n_events,
                no_strong_Pauli_or_BGK_or_SM_FWR_or_relic=True)


def four_dirac_sectors_electron_heat(
    *,Tc1=.8,Tc2=.9,Te=1.,gD=.1,epsilon=1e-5,mV=10.,
    n=8,cutoff=20.):
    """χ1, antiχ1, χ2, antiχ2 zero chemical asymmetry, common e± bath.

    All four independent elastic number sources vanish; total QSM=-ΣQdark
    is an INSTANTANEOUS weak event identity, not finite SM heat evolution.
    """
    charges=[]
    for name,m,T in (
            ("chi1",MCHI1,Tc1),("anti_chi1",MCHI1,Tc1),
            ("chi2",MCHI2,Tc2),("anti_chi2",MCHI2,Tc2)):
        q=one_dirac_charge_sector_weak_heat(
            mchi=m,Tchi=T,Te=Te,gD=gD,epsilon=epsilon,mV=mV,
            n_p=n,n_k=n,n_mu=n,n_cm=n,n_az=n,cutoff=cutoff)
        q["species"]=name
        charges.append(q)
    Qdark=sum(x["Q_chi_GeV5"] for x in charges)
    Qsm=sum(x["Q_e_GeV5"] for x in charges)
    H=sum(x["quantum_H_GeV4"] for x in charges)
    return dict(schema_version=1,scope="real-phi gauge-compatible both-Dirac finite-vector e± Pauli 5D *instantaneous weak heat* only, not full SM or FRW",
                mV_GeV=mV,Te_GeV=Te,
                Q_dark_from_e_GeV5=Qdark,Q_e_bath_GeV5=Qsm,
                exact_heat_ledger_Qdark_plus_Qe_GeV5=Qdark+Qsm,
                H_total_GeV4=H,components=charges,
                number_fluxes_all_four_zero=True,class_ready=False)


def audit():
    results=[]
    for n in (6,8,10,12):
        results.append(four_dirac_sectors_electron_heat(n=n,cutoff=20.))
    return dict(grid_rows=results,
                scope="physical finite-mV positive FD/Pauli vector electron-only weak heat, not SM+FRW/kinetic drag/BE relic",
                class_ready=False)


if __name__=="__main__":print(json.dumps(audit(),indent=2))
