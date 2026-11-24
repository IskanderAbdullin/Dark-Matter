"""Physical linear chemical-relaxation eigenmodes near the real-scalar
Dirac Yukawa FULL four-reaction FD/BE equilibrium, conditional fixed T.

A chemical fugacity perturbation around exact common-temperature equilibrium
with mu_chi1=mu_chi2=mu_U1, mu_anti=-mu_U1, mu_phi=0 has:
  net microscopic event rate R_r = - G_r (S_r . delta_mu) / T + O(mu²)
where G_r = ∫ physical positive common on-shell dW_r F_eq = ∫W_r R_eq,
and S_r is the stoichiometric row (products minus reactants).

For physical diagonal equilibrium susceptibilities
  chi_i=∂n_i/∂mu_i= g_i/T ∫p²dp/(2π²)f_i(1∓f_i),
the symmetric chemical-relaxation matrix
  L=chi^(-1/2) S^T diag(G/T) S chi^(-1/2)
is positive semidefinite with exactly the true vector Dirac U1 zero mode
on a chemically connected continuum network.

The zero eigenvalue is exact stoichiometry. POSITIVE finite eigenvalues
in a restricted thermal fugacity five-parameter ansatz are conditional
linear chemical equilibration diagnostics, NOT proof of full momentum
thermalization or Boltzmann+Friedmann freeze-out. Hard parent P cuts can
leave continuum scalar modes collision-disconnected even if L>0.
"""
import json
import math
import numpy as np
from scipy.integrate import quad
from scipy.special import expit
from kinetics.phi_bose_five_dirac_full_yukawa_weak import FullDiracYukawaWeak
from kinetics.phi_bose_yukawa_charge_audit import stoichiometry
from kinetics.phi_bose_physical_cutcell_weak import M1,M2,MPHI

MPLANCK_GEV=1.2209e19
Q_U1=np.array([1.,-1.,1.,-1.,0.])


def susceptibility(T,mass,g,statistics,mu=0.):
    """Physical thermal FD/BE per-species number susceptibility, GeV²."""
    if T<=0 or mass<0 or g<=0:raise ValueError("Nonpositive physical input")
    if statistics=="boson" and mu>=mass:
        raise ValueError('Subcritical Bose chemical potential required')
    m=mass/T
    muhat=mu/T
    def f(x):
        E=math.hypot(x,m)-muhat
        if statistics=="fermion":
            occ=float(expit(-E))
            var=occ*(1-occ)
        elif statistics=="boson":
            occ=1/math.expm1(E)
            var=occ*(1+occ)
        else:raise ValueError("Unknown physical statistics")
        return x*x*var
    return g*T*T/(2*math.pi**2)*quad(f,0.,45.,
          epsabs=2e-12,epsrel=3e-11,limit=240)[0]


def equilibrium_event_rate_weight(channel,temperature,muU=0.):
    """Physical on-shell equilibrium ∫W F, with particle/anti potentials.

    For each event F=R by χ1/χ2 common Dirac U1 and real φ mu=0.
    """
    T=temperature
    if T<=0:raise ValueError("Positive T")
    spec={
       "chi1_decay":None,"anti_chi1_decay":None,
       "chi1_pair":None,"chi2_pair":None}
    if channel.name=="chi1_decay":
        E1,E2,Ep=channel.energy
        W=channel.W
        f=expit(-(E1-muU)/T)
        g=expit(-(E2-muU)/T)
        b=1/np.expm1(Ep/T)
        return float(np.dot(W,f*(1-g)*(1+b)))
    if channel.name=="anti_chi1_decay":
        E1,E2,Ep=channel.energy
        W=channel.W
        f=expit(-(E1+muU)/T)
        g=expit(-(E2+muU)/T)
        b=1/np.expm1(Ep/T)
        return float(np.dot(W,f*(1-g)*(1+b)))
    if channel.name in ("chi1_pair","chi2_pair"):
        Echi,Eanti,E3,E4=channel.energy
        W=channel.W
        f=expit(-(Echi-muU)/T)
        g=expit(-(Eanti+muU)/T)
        b=1/np.expm1(E3/T)
        c=1/np.expm1(E4/T)
        return float(np.dot(W,f*g*(1+b)*(1+c)))
    raise ValueError("Unknown physical Yukawa channel "+channel.name)


def equilibrium_chemical_relaxation(model,T=1.,gstar=61.75,muU=0.):
    """Direct physical event weights and susceptibility matrices."""
    if T<=0 or gstar<=0:raise ValueError("Physical positive T and gstar")
    # Expected model is the fully physical five-field Bernstein source,
    # but only its channel exact event energies/W are used; source weak.
    ms=(M1,M1,M2,M2,MPHI)
    ss=("fermion","fermion","fermion","fermion","boson")
    gs=(2,2,2,2,1)
    chem=(muU,-muU,muU,-muU,0.)
    chi=np.array([susceptibility(T,m,g,k,mu) for m,g,k,mu in
                  zip(ms,gs,ss,chem)])
    gamma=np.array([equilibrium_event_rate_weight(c,T,muU)
                    for c in model.channels])
    if np.min(gamma)<=0 or np.min(chi)<=0:
        raise ArithmeticError("Physical rate/thermodynamic χ nonpositive")
    # Note: channel order of FullDiracPhysicalMoments equals
    # [χ1 decay, antiχ1 decay, χ1 pair, χ2 pair], matching exact S.
    S=stoichiometry(True).astype(float)
    K=S.T @ np.diag(gamma/T) @ S
    L=K/np.sqrt(np.outer(chi,chi))
    spectrum=np.linalg.eigvalsh((L+L.T)/2)
    H=1.66*math.sqrt(gstar)*T*T/MPLANCK_GEV
    uv=np.sqrt(chi)*Q_U1
    u1_rel_residual=float(np.linalg.norm(L@uv)/(max(spectrum)*np.linalg.norm(uv)))
    return dict(T_GeV=T,P_parent_GeV=model.P,P_pair_GeV=model.Ppair,
                gstar_rho=gstar,H_RD_GeV=H,
                susceptibility_GeV2=chi.tolist(),
                channel_equilibrium_F_integrated_GeV4=gamma.tolist(),
                chemical_relaxation_eigenvalues_GeV=np.sort(spectrum).tolist(),
                chemical_relaxation_over_H=(np.sort(spectrum)/H).tolist(),
                true_vector_U1_symmetric_null_residual=u1_rel_residual,
                scope="fixed-T, 5-species thermal FD/BE fugacity response; cannot establish thermalization of all momentum bins or actual χ-SM/FRW relic",
                class_ready=False)


def audit():
    rows=[]
    for P,nr,nang,pnr,prel,pcm,naz in (
        (2.,4,8,4,8,32,4),
        (8.,6,12,6,8,48,4),
        (16.,8,12,8,12,64,6)):
        from kinetics.phi_bose_five_dirac_full_yukawa_implicit import FullDiracPhysicalMoments
        m=FullDiracPhysicalMoments(P=P,Ppair=P,nrad=nr,
                  nangle=nang,pair_nrad=pnr,pair_nrelative=prel,
                  pair_ncm=pcm,pair_naz=naz)
        rows.append(equilibrium_chemical_relaxation(m,T=1.,gstar=61.75))
    return dict(schema_version=1,
                model="two Dirac+anti fields plus real φ, four Yukawa reactions, y=.001, vacuum masses .1/.02/86e−9 GeV, thermal FD/BE at zero μ",
                rows=rows,class_ready=False)

if __name__=="__main__":
    print(json.dumps(audit(),indent=2,allow_nan=False))
