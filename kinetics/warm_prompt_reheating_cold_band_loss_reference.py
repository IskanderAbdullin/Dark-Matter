"""Compare late stable count loss with separate soft and hard scalar occupation caps."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from three_species_decay_onsager_reference import m, mp
from warm_late_visible_recycling_envelope_reference import (
    COUP2, VECTOR_GAMMA, VECTOR_MASS,
)


def retention(a, b):
    exposure = -np.expm1(-b) / b if b else 1.0
    return np.exp(-b) / (1.0 + a * exposure)


def vector_energy_ceiling():
    p = VECTOR_MASS**2
    w2 = (VECTOR_MASS * VECTOR_GAMMA)**2
    s = (p + np.sqrt(4*p*p + 3*w2))/3
    coeff = 3*COUP2/(4*np.pi)

    def rate(inv):
        return coeff*np.sqrt(inv)/((inv-p)**2 + w2)

    mesh = np.r_[np.geomspace(4*m*m, 1e7, 5001),
                 p + np.sqrt(w2)*np.linspace(-5, 5, 5000)]
    assert np.max(rate(mesh)) <= rate(s)*(1+2e-14)
    return float(rate(s)), float(s)


def test_kinematic_gate():
    rng = np.random.default_rng(720004)
    s = 4*m*m*np.exp(rng.uniform(0., 12., 16000))
    boost = np.exp(rng.uniform(0., 10., s.shape))
    root = np.sqrt(s)
    beta = np.sqrt(1-boost**-2)
    speed = np.sqrt(1-4*mp*mp/s)
    angle = rng.uniform(-1., 1., s.shape)
    E = boost*root
    e1 = E/2*(1+beta*speed*angle)
    e2 = E/2*(1-beta*speed*angle)
    ratio = np.minimum(e1, e2)*4*E/s
    assert np.min(ratio) >= 1-2e-14
    assert np.all(np.maximum(e1, e2) >= m)
    # Nearly saturate the analytic low-daughter inequality without cancellation.
    gamma = 310.
    s0 = 4*m*m
    beta0 = np.sqrt(1-gamma**-2)
    v0 = np.sqrt(1-4*mp*mp/s0)
    # Exact rationalized ratio 2*(1+gamma² beta² 4mp²/s)/(1+beta*v).
    near = 2*(1+gamma**2*beta0**2*4*mp*mp/s0)/(1+beta0*v0)
    assert 1 < near < 1.00002
    return dict(sample_count=len(s), min_sampled_ratio=float(np.min(ratio)),
                near_forward_corner_ratio=float(near))


def run(floor, isolated, visible_input):
    fl = json.loads(floor.read_text())
    en = json.loads(isolated.read_text())
    vi = json.loads(visible_input.read_text())
    k0 = fl['invariant_vacuum_upper_GeVminus2']
    kv, peak_s = vector_energy_ceiling()
    check = test_kinematic_gate()
    rows = []
    for data, energy, budget in zip(fl['anchors'], en['anchors'], vi['anchors'],
                                    strict=True):
        T = data['T_visible_GeV']
        assert T == energy['visible_T'] == budget['T_visible_GeV']
        H = energy['photon_H_floor_at_anchor']
        rho0 = energy['initial_dark_energy_density']
        rho = budget['remaining_energy_ceiling']
        assert rho >= rho0
        C = data['exposure_coefficient_per_Bose_multiplier']
        raw = data['initial_Omega_h2_if_counts_freeze']
        Bv = 2*kv*rho/H
        old_Bv = data['uniform_vector_exposure_coefficient']
        assert Bv < old_Bv
        for hard in (0., 1., 10.):
            a = C*(1+hard)**2
            for split in (1e-6, 1e-5, 1e-4, 1e-3):
                assert mp < split < m
                threshold = m*m/split
                b_per_excess = 2*k0*rho*(1+hard)/(H*threshold)
                def bound(f):
                    b_s = b_per_excess*(f-hard)
                    return raw*max(retention(a, Bv+b_s),
                                   retention(a+old_Bv, b_s))
                necessary = (brentq(lambda f: bound(f)-.12,
                                    hard, hard+100./b_per_excess)
                             if bound(hard) > .12 else None)
                samples = []
                for soft in (hard, max(100.,hard), max(10000.,hard),
                             max(1000000.,hard)):
                    b = Bv+b_per_excess*(soft-hard)
                    retain = float(retention(a, b))
                    # Choose either available independent bound, never make it weaker.
                    old = raw/(1+old_Bv+C*(1+soft)**2)
                    samples.append(dict(soft_cap=soft, energy_gate_omega_floor=raw*retain,
                                        old_uniform_omega_floor=old,
                                        combined_omega_floor=max(old, bound(soft))))
                rows.append(dict(T_visible_GeV=T, split_GeV=split,
                                 necessary_pair_energy_GeV=threshold, assumed_hard_cap=hard,
                                 quadratic_exposure=a, linear_vector_exposure=Bv,
                                 old_uniform_vector_exposure=old_Bv,
                                 linear_scalar_exposure_per_excess_occupation=b_per_excess,
                                 necessary_soft_cap_for_this_floor_to_reach_point12=necessary,
                                 samples=samples))
        for a, b in ((0., 0.), (1., 0.), (.2, min(Bv,8.)), (.3, 5.), (2., 14.)):
            expected = retention(a,b)
            rr = solve_ivp(lambda u, y: -a*y*y-b*y, (0.,1.), [1.],
                           rtol=1e-11, atol=1e-14)
            assert rr.success and abs(rr.y[0,-1]/expected-1) < 1e-8
    return dict(scope="Modeled anchors and stated charged-lepton plus scalar-pair stable losses only. "
                      "Dark energy bounded by the separate selected vacuum lepton energy allowance; "
                      "extra sources, other stable losses and hadrons not covered. "
                      "Caps are assumptions, not derived occupation constraints.",
                finite_width_lepton_ceiling_per_pair_energy_GeVminus3=kv,
                finite_width_peak_s_GeV2=peak_s, kinematic_controls=check, results=rows)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('floor', type=Path)
    p.add_argument('isolated', type=Path)
    p.add_argument('visible_input', type=Path)
    a = p.parse_args()
    print(json.dumps(run(a.floor, a.isolated, a.visible_input), indent=2, allow_nan=False))
