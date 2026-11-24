"""Conditional thermal 86-eV real scalar relic abundance and free-streaming gate.

A FREE collisionless real scalar with zero chemical potential, no condensate,
instant relativistic freeze-out and subsequent adiabatic SM entropy is used.
This is NOT the nonthermal self-consistent Yukawa+SM+FRW solution.

External inputs (not inferred from the Yukawa solver):
 - Planck 2018 VI Ω_c h²≈0.1200, Ω_m≈0.315, z_eq≈3402, h≈0.674.
 - PDG photon number ≈410.7/cm³ and rho_critical/h²≈1.053672e4 eV/cm³,
   Tcmb≈2.7255K; g*S,SM(today)≈3.909.
 - g*S,SM(1GeV)≈61.75 is a rough step-function illustrative value.

Compare only CHARACTERISTIC mean collisionless path lengths, not a CLASS
transfer-function half-mode mass or a direct Lyman-alpha 95% exclusion.
"""
import math
import json

import numpy as np
from scipy.integrate import quad_vec
from scipy.special import roots_laguerre,zeta

T_GAMMA_K=2.7255
EV_PER_K=8.617333262e-5
PHOTON_DENSITY_CM3=410.7
RHOCRIT_H2_EV_CM3=10536.72
GSTARS0=3.909
OMEGA_C_H2=0.1200
OMEGA_M=0.315
Z_EQ=3402.0
H0_KM_S_MPC=67.4
C_KM_S=299792.458
ZS3=float(zeta(3.,1.))


def thermal_number_factor(kind):
    if kind=="real_boson":
        return .5       # relative to photon n_gamma for one BE degree
    if kind=="two_fermion":
        return .75      # relative to photon n_gamma for two FD degrees
    raise ValueError("Only real Bose g=1 and Fermi g=2 trials supported")


def temperature_ratio_for_total_dark_matter(mass_eV,kind,omega_h2=OMEGA_C_H2):
    """Required T_relic/T_gamma for a completely thermal ALL-DM component."""
    if not (mass_eV>0 and omega_h2>0):
        raise ValueError("Positive mass and abundance required")
    return (omega_h2*RHOCRIT_H2_EV_CM3/
            (mass_eV*PHOTON_DENSITY_CM3*thermal_number_factor(kind)))**(1/3)


def thermal_abundance(mass_eV,gstar_entropy_dec,
                      gstar_entropy_today=GSTARS0):
    """g=1, μ=0, no condensate and zero later entropy transfer into scalar."""
    if mass_eV<=0 or gstar_entropy_dec<=gstar_entropy_today:
        raise ValueError("Relativistic hot decoupling assumed")
    t_ratio=(gstar_entropy_today/gstar_entropy_dec)**(1/3)
    number=PHOTON_DENSITY_CM3*.5*t_ratio**3
    omega=mass_eV*number/RHOCRIT_H2_EV_CM3
    t0=t_ratio*T_GAMMA_K*EV_PER_K
    nu_ratio=(4/11)**(1/3)
    delta_neff=(4/7)*(t_ratio/nu_ratio)**4
    return dict(mass_eV=mass_eV,gstar_entropy_dec=gstar_entropy_dec,
                Tphi_over_Tgamma=t_ratio,
                Tphi0_eV=t0,nphi_cm_minus3=number,
                Omega_phi_h2=omega,delta_Neff_early_thermal=delta_neff,
                gstar_entropy_dec_for_omega_c_h2=
                    gstar_entropy_today*(mass_eV*PHOTON_DENSITY_CM3*.5)/
                    (OMEGA_C_H2*RHOCRIT_H2_EV_CM3))


def thermal_momentum_weights(kind,order=96):
    """Independent Gauss-Laguerre evaluation of exact FD/BE number measure.

    Number-weighted normalized P(x)dx∝x²dx/(exp(x)∓1), x=p/T0.
    """
    x,w=roots_laguerre(order)
    if kind=="real_boson":
        weight=w*x*x/(-np.expm1(-x))
    elif kind=="two_fermion":
        weight=w*x*x/(1+np.exp(-x))
    else:
        raise ValueError("Unsupported statistics")
    weight/=np.sum(weight)
    return x,weight


def mean_collisionless_horizon(mass_eV,T0_eV,kind,
                               order=96,zeq=Z_EQ,Om=OMEGA_M,H0=H0_KM_S_MPC):
    """Number-weighted free-streaming COMOVING path until matter-radiation equality.

    Under a_dec≈0, homogeneous free-streaming and standard RD+MD:
      λ(q)= (c/H0)∫_0^a_eq da [q/sqrt(q²+(m a)²)]
                     /sqrt(Ωr+Ωm a).
    This is *not* the linear-theory matter power transfer function.
    """
    if not (mass_eV>0 and T0_eV>0 and 0<Om<1 and zeq>1):
        raise ValueError("Invalid RD/MD free-streaming parameters")
    x,w=thermal_momentum_weights(kind,order)
    r=x*T0_eV/mass_eV
    ae=1/(1+zeq)
    Or=Om*ae  # equality condition under this RD+MD approximation
    integrals,err=quad_vec(
        lambda a: r/np.hypot(r,a)/np.sqrt(Or+Om*a),
        0,ae,epsrel=2e-10)
    mean_speed_eq=float(np.dot(w,r/np.hypot(r,ae)))
    return dict(mass_eV=mass_eV,kind=kind,
                T0_eV=T0_eV,
                T_over_T_gamma=T0_eV/(T_GAMMA_K*EV_PER_K),
                number_weighted_mean_x=float(np.dot(w,x)),
                mean_velocity_eq_over_c=mean_speed_eq,
                mean_velocity_today_over_c=float(np.dot(w,r/np.hypot(r,1.))),
                mean_comoving_freestream_to_eq_Mpc=float(
                    C_KM_S/H0*np.dot(w,integrals)),
                velocity_quadrature_order=order,
                integral_error_estimate=float(np.max(np.abs(err))))


def audit():
    baseline=thermal_abundance(86.,61.75)
    v_baseline=mean_collisionless_horizon(
        86.,baseline["Tphi0_eV"],"real_boson",128)
    t_all=temperature_ratio_for_total_dark_matter(86.,"real_boson")
    v_all=mean_collisionless_horizon(
        86.,t_all*T_GAMMA_K*EV_PER_K,"real_boson",128)
    t_wdm=temperature_ratio_for_total_dark_matter(5700.,"two_fermion")
    v_wdm=mean_collisionless_horizon(
        5700.,t_wdm*T_GAMMA_K*EV_PER_K,"two_fermion",128)
    return dict(schema_version=1,
                scope="strict conditional relativistic thermal decoupled scalar and characteristic free-streaming; not solved Yukawa relic, not CLASS/likelihood",
                early_decoupling_real_boson_86eV_gstar61p75=baseline,
                early_decoupling_free_streaming=v_baseline,
                all_DM_real_boson_86eV=v_all,
                comparison_all_DM_thermal_FD_two_fermion_5p7keV=v_wdm,
                mean_comoving_path_ratio_scalar_to_WDM=
                   v_all["mean_comoving_freestream_to_eq_Mpc"]/
                   v_wdm["mean_comoving_freestream_to_eq_Mpc"],
                references=[
                    "https://doi.org/10.1051/0004-6361/201833910",
                    "https://pdg.lbl.gov/2024/download/db2024.pdf",
                    "https://doi.org/10.1103/PhysRevD.109.043511"
                ],class_ready=False)


if __name__=="__main__":
    print(json.dumps(audit(),indent=2,allow_nan=False))
