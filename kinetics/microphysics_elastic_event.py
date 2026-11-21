"""A microscopic fermion-fermion 2->2 EVENT ledger for the chosen portal.

The Pauli factors are both-species full gain/loss for chi+e <-> chi'+e'.
It checks equilibrium and exact energy-exchange signs per event. It does
NOT perform the full Lorentz-invariant phase-space integral, handle bath
backreaction numerically, or substitute for gamma_el.

Use natural units (GeV for energies and T) and strictly positive, finite
occupations. e is treated as massless and the initial chi is at rest in
the provided kinematics reference.
"""
from __future__ import annotations
import math
from microphysics_thermal_gain_loss import fermi


def rest_chi_event_kinematics(mchi,k,z):
    """Physical chi at rest, massless e: interpolate t from 0 to t_min.

    Returns lab incoming/outgoing energies and scattering angle cosine,
    plus exact energy and 3-momentum shell residuals.
    """
    if (not all(math.isfinite(v) for v in (mchi,k,z)) or
            mchi<=0 or k<=0 or not 0<=z<=1):
        raise ValueError("positive mchi,k and z in [0,1] required")
    s=mchi*mchi+2*mchi*k
    tmin=-(s-mchi*mchi)**2/s
    t=z*tmin
    kout=k+t/(2*mchi)
    if kout<=0: raise ArithmeticError("nonpositive final electron energy")
    echi_out=mchi-t/(2*mchi)
    costheta=1+t/(2*k*kout)
    if not -1-1e-12<=costheta<=1+1e-12:
        raise ArithmeticError("unphysical scattering angle")
    pout_squared=k*k+kout*kout-2*k*kout*costheta
    shell=echi_out**2-mchi*mchi-pout_squared
    return {
      "t":t,"s":s,"chi_before_E":mchi,"chi_after_E":echi_out,
      "electron_before_E":k,"electron_after_E":kout,
      "electron_cos_theta_lab":costheta,
      "chi_recoil_energy":echi_out-mchi,
      "energy_sum_residual":(echi_out+kout)-(mchi+k),
      "mother_mass_shell_residual":shell,
    }


def elastic_fermion_event(Echi,Ee,Echi_out,Ee_out,Tchi,Te,muchi=0.,mue=0.):
    """Unnormalised 2->2 quantum forward/inverse probabilities.

    fchi*f_e*(1-fchi')*(1-f_e') -
    fchi'*f_e'*(1-fchi)*(1-f_e)

    No amplitude, Jacobian, delta function or density-of-states factor.
    """
    vals=(Echi,Ee,Echi_out,Ee_out,Tchi,Te,muchi,mue)
    if not all(math.isfinite(v) for v in vals) or min(Echi,Ee,Echi_out,Ee_out)<0:
        raise ValueError("Finite nonnegative energies and finite parameters required")
    if Tchi<=0 or Te<=0:
        raise ValueError("Positive temperatures required")
    if not math.isclose(Echi+Ee,Echi_out+Ee_out,rel_tol=1e-12,abs_tol=1e-12):
        raise ValueError("Elastic 2->2 physical energy conservation required")
    fi=fermi(Echi,Tchi,muchi)
    gi=fermi(Ee,Te,mue)
    fo=fermi(Echi_out,Tchi,muchi)
    go=fermi(Ee_out,Te,mue)
    forward=fi*gi*(1-fo)*(1-go)
    backward=fo*go*(1-fi)*(1-gi)
    net=forward-backward
    delta=Echi_out-Echi
    # Equilibrium theorem from log(f/(1-f))=-(E-mu)/T.
    # For chi heating (delta>0), colder chi gains energy if Te>Tchi.
    logratio=delta*(1/Tchi-1/Te)
    # No log(0) required even if probabilities underflow.
    entropy_density=net*logratio
    return {
      "forward_factor":forward,"reverse_factor":backward,"net_forward_factor":net,
      "expected_log_forward_reverse":logratio,
      "chi_energy_transfer_unweighted":delta*net,
      "electron_energy_transfer_unweighted":-delta*net,
      "total_energy_transfer_unweighted":0.,
      "chi_number_transfer_unweighted":0.,
      "electron_number_transfer_unweighted":0.,
      "entropy_production_unweighted":entropy_density,
      "not_full_2to2_collision_integral":True
    }
