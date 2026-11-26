"""Linearized χ thermal ENERGY-moment relaxation from full Pauli 5D 2->2.

Elastic scattering conserves Nχ. Along a restricted FD(Tχ,μχ) ansatz,
vary μχ(Tχ) to keep Nχ fixed, and calculate c_N=dρχ/dTχ|N.
The resulting gamma_T= -[dQχ/dTχ|N,Tχ=Te] / c_N is the relaxation rate
of ONE thermal-energy moment on the FD trial manifold.

This is a real collision-derived scalar under the specified electron
vector EFT, but emphatically NOT a momentum-dependent gamma_el(p,T),
not a true linearized Boltzmann spectral eigenvalue, and not a physical
freezeout rate for arbitrary nonequilibrium χ. A finite SM bath,
gauge-required neutrinos and plasma/UV corrections are not evolved.
Natural units with T,M,μ and gamma in GeV, G in GeV^-2.
"""
from __future__ import annotations
import math
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import expit
from relativistic_weak_form_2to2 import weak_energy_moments


def FD_thermo_moments(M,T,mu,gchi=4):
    if (not all(math.isfinite(v) for v in (M,T,mu,gchi)) or
            M<=0 or T<=0 or gchi<=0):
        raise ValueError("M,T,gchi>0 and finite mu")
    def vals(p):
        E=math.hypot(M,p)
        f=float(expit((mu-E)/T))
        z=f*(1-f)
        pref=gchi/(2*math.pi**2)*p*p
        return (pref*f,pref*E*f,pref*z*(E-mu)/(T*T),
                pref*z/T,pref*E*z*(E-mu)/(T*T),pref*E*z/T)
    v=[]
    for i in range(6):
        v.append(quad(lambda p:vals(p)[i],0,60*max(T,M),
                      epsabs=2e-11,epsrel=2e-10,limit=200)[0])
    N,rho,NT,Nmu,rhoT,rhomu=v
    slope=-NT/Nmu if Nmu>0 else float("nan")
    heat=rhoT-rhomu*NT/Nmu if Nmu>0 else float("nan")
    return {"N":N,"rho":rho,"N_T":NT,"N_mu":Nmu,
            "rho_T":rhoT,"rho_mu":rhomu,
            "mu_slope_at_fixed_N":slope,"c_N":heat}


def matched_mu(M,T,Ntarget,gchi=4):
    if not math.isfinite(Ntarget) or Ntarget<=0:
        raise ValueError("positive finite Ntarget required")
    def delta(mu):
        return FD_thermo_moments(M,T,mu,gchi)["N"]-Ntarget
    return brentq(delta,-40*T-M,M+40*T,xtol=1e-12,rtol=1e-12)


def thermal_relaxation_rate(M=.1,Te=1.,mu0=0.,G=.001,gchi=4.,ge=4.,
                             delta_T=.02,n=20,cutoff=32.):
    vals=(M,Te,mu0,G,gchi,ge,delta_T,cutoff)
    if (not all(math.isfinite(x) for x in vals) or M<=0 or Te<=0 or
        G<0 or gchi<=0 or ge<=0 or not 0<delta_T<Te*.25):
        raise ValueError("physical inputs and sufficiently small delta_T required")
    ref=FD_thermo_moments(M,Te,mu0,gchi)
    if ref["c_N"]<=0 or not math.isfinite(ref["c_N"]):
        raise ArithmeticError("nonpositive or invalid c_N")
    qs=[];mus=[];drifts=[]
    for sign in (-1,1):
        T=Te+sign*delta_T
        mu=matched_mu(M,T,ref["N"],gchi)
        m=FD_thermo_moments(M,T,mu,gchi)
        q=weak_energy_moments(
            M=M,Tchi=T,Te=Te,muchi=mu,mue=0.,G=G,
            gchi=gchi,ge=ge,n_p=n,n_k=n,n_mu=n,n_cm=n,n_phi=n,
            cutoff=cutoff)["Q_chi_GeV5"]
        qs.append(q);mus.append(mu)
        drifts.append(m["N"]/ref["N"]-1)
    # Thermal energy moment is cooling: Qχ(Tχ>Te) must be negative.
    K=-(qs[1]-qs[0])/(2*delta_T)
    return {
      "baseline_N_GeV3":ref["N"],
      "baseline_rho_GeV4":ref["rho"],
      "fixed_number_heat_capacity_cN_GeV3":ref["c_N"],
      "chemical_mu_slope_per_temperature":ref["mu_slope_at_fixed_N"],
      "mus_minus_plus_GeV":mus,
      "number_residual_minus_plus_fraction":drifts,
      "Q_chi_minus_plus_GeV5":qs,
      "K_linear_heat_conductance_GeV4":K,
      "gamma_T_energy_moment_GeV":K/ref["c_N"],
      "is_true_spectral_eigenvalue_or_BGK_gamma":False,
      "physical_scope":"linearized FD-family chi thermal energy moment for electron-only 5D Pauli EFT; fixed SM, not full kinetic spectrum or relic",
    }
