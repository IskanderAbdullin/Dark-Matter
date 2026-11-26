"""Conditional p-wave depletion bound from a finite comoving incident support."""
import argparse
import json
from pathlib import Path

import numpy as np

from archived_scalar_pair_reaction_reference import amplitude
from sm_dark_entropy_exchange_reference import MPL
from three_species_decay_onsager_reference import M, m, mp, y


def coefficient(P):
    # t,u<=0, thus |(t-M^2)(u-M^2)| >= M^4. In the expression
    # for the spin-averaged amplitude, |x sqrt(1-x^2)| <= 1/2.
    K2 = m*m + P*P
    A = M*M + m*m + 2*P*P
    d = M**4
    transverse_over_p = 2*K2/d
    longitudinal_over_p = 2*m*(A+2*K2)/d + 2*M*A/d
    return 2*y**4*(K2*transverse_over_p**2 + longitudinal_over_p**2)


def test(P, C):
    p = np.geomspace(max(P*1e-6,1e-9),P,91)[:,None]
    a = np.linspace(-1.,1.,301)
    s = 4*(p*p+m*m)
    # Use p reconstructed from s for the check at very small p.
    pr = np.sqrt(s/4-m*m)
    r = amplitude(s,a)/pr**2/C
    assert np.isfinite(r).all() and r.min() >= 0 and r.max() <= 1.
    return float(r.max())


def report(shape, stable_omega, exterior_omega, target=.120, photon_factor=1.):
    a = shape["scale_factor"]
    P = shape["lab_momentum_edges_GeV"][-1]
    assert 0 <= exterior_omega < stable_omega
    n_ret = shape["represented_species_count_a3"][1]/a**3
    n = n_ret/(1.-exterior_omega/stable_omega)
    T = shape["T_visible_GeV"]
    C = coefficient(P)
    H = np.sqrt(np.pi**2/45)/MPL*(T*photon_factor)**2
    assert 0 < photon_factor <= 1 and stable_omega > target > 0
    # For two distinguishable charge populations totaling N at zero
    # asymmetry, each forward reaction removes two particles. A charge
    # asymmetry only lowers the product at fixed N.
    # sigma*v_M <= C P^2 b^-2 (1+fmax)^2 / (64 pi m^2).
    # The incident density is <= current inclusive N/b^3, even when
    # some particles are exterior. Thus d(1/N)/dlog(b) <= A b^-3.
    # H >= H_photon b^-2 and integral b^-3 dlog b = 1/3.
    tau = n*C*P*P/(384*np.pi*m*m*H)
    necessary = stable_omega/target-1.
    caps = [0.,.001,1.,2.]
    return dict(temperature_GeV=T, support_upper_GeV=P, initial_inclusive_stable_number_GeV3=float(n),
            spin_averaged_amplitude2_over_pCM2_bound_GeV_minus2=C,
            sampled_max_amplitude_fraction_of_analytic_bound=test(P,C),
            assumed_minimum_comoving_photon_temperature_ratio=photon_factor,
            quadratic_loss_integral_without_Bose_enhancement=float(tau),
            stable_fixed_count_Omega_h2_at_anchor=stable_omega,
            target_fixed_count_Omega_h2=target,
            required_normalized_quadratic_loss_integral=float(necessary),
            necessary_supremum_of_future_scalar_occupation_to_allow_target=float(np.sqrt(necessary/tau)-1.),
            conditional_cases=[dict(scalar_occupation_upper_cap=z,quadratic_loss_integral_bound=float(tau*(1+z)**2),
                final_stable_fixed_count_Omega_h2_lower=float(stable_omega/(1.+tau*(1.+z)**2))) for z in caps])


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("shape",type=Path)
    p.add_argument("count_audit",type=Path)
    p.add_argument("--photon-factor",type=float,default=1.)
    a=p.parse_args()
    s=json.loads(a.shape.read_text())["snapshots"][-1]
    r=json.loads(a.count_audit.read_text())["chained_endpoint_cases"][-1]
    assert s["steps"]==r["n_steps"] and abs(s["T_visible_GeV"]/r["temperature_GeV"]-1.)<1e-12
    v=report(s,r["retained_and_all_particle_tail_Omega_h2"][1],r["all_tail_Omega_h2"][1],photon_factor=a.photon_factor)
    v["scope"]="Post-anchor stable pair removal only. No new incident populations or impacts from exterior re-entry; future comoving incident support stays inside the archived support. Photon floor and scalar occupation cap assumed, not established by the snapshot. Dark self-scattering, number-changing processes or portal depletion that invalidate these conditions are outside the bound."
    print(json.dumps(v,indent=2,allow_nan=False))
