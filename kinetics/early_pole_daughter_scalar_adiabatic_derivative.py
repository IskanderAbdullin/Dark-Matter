"""Homogeneous stress susceptibility at fixed bare mass and frozen comoving shape."""
import argparse
import json

import numpy as np

from early_pole_daughter_bare_scalar_stress_control import state
from early_pole_daughter_retained_kinetic_stress import build, r


def check(q, n, lam, bare):
    x = state(q, n, lam, bare)
    m = x["effective_gap"]
    E = np.hypot(q, m)
    w = x["pressure_to_energy"]
    I = float(n @ (1/E))
    A = float(n @ (1/E**3))
    dm = .5*lam*(2*I + m*m*A)/(2*m + .5*lam*m*A)
    e_log = (q*q + m*dm)/(E*E)
    rp = 3*(x["energy"] + x["pressure"])
    kinetic = q*q/(3*E)
    pp = float(n @ (kinetic*(5-e_log))) + m*I*dm
    i_log = 3 - float(n @ (e_log/E))/I
    trace_form = (1 - bare*bare*I*i_log/rp)/3
    lower, upper = 4*w/(3*(1+w)), min(1/3, 2*w/(1+w), (1+9*w)/(9*(1+w)))
    if not (0 < dm/m <= 1+2e-12 and 2-2e-12 <= i_log <= 3+2e-12
            and lower-2e-12 <= pp/rp <= upper+2e-12 and abs(trace_form-pp/rp)<2e-12):
        raise ValueError("Homogeneous trace derivative or analytic bounds did not close.")
    h = 1e-3
    scales = np.exp(h*np.array([-2, -1, 1, 2]))
    direct = [state(q, n, lam, bare, s) for s in scales]
    z = np.array([[v[i] for i in ("energy", "pressure")] for v in direct])
    finite = (z[0] - 8*z[1] + 8*z[2] - z[3])/(12*h)
    if abs(finite[0]-rp)>1e-8*rp or abs(finite[1]-pp)>1e-7*pp:
        raise ValueError("Analytic susceptibility did not match fixed-bare dilation.")
    # Independent nonlinear identity under a change of variables with b fixed.
    for s, d in zip(scales, direct, strict=True):
        scaled = state(q, n, lam, bare/s)
        for k, pow_s in (("effective_gap", 1), ("energy", 4), ("pressure", 4)):
            if abs(d[k]-s**pow_s*scaled[k])>2e-11*d[k]:
                raise ValueError("Finite dilation scaling did not reproduce.")
    return dict(bare_mass=float(bare), effective_gap=m, w=w, adiabatic_derivative=pp/rp,
                trace_form_derivative=trace_form, bound_lower=lower, bound_upper=upper,
                pressure_log_slope=pp/x["pressure"], gap_log_slope=dm/m, scalar_density_log_slope=i_log,
                finite_difference_adiabatic_derivative=float(finite[1]/finite[0]))


def fixture():
    values = []
    for b in (0., .3, .85, .99):
        # Set m=1 via the exact zero-momentum gap. lambda=2.
        actual = check(np.zeros(1), np.array([1-b*b]), 2., b)
        w = (1-b*b)/(3+b*b)
        c = (1-b*b)/(3-b*b)
        for key, expected in (("w", w), ("adiabatic_derivative", c)):
            if abs(actual[key]-expected)>2e-13:
                raise ValueError("Exact zero-momentum fixture disagreed.")
        values.append(dict(expected_w=w, expected_adiabatic_derivative=c, **actual))
    return values


def run(archive, arms, bare_receipt, angular=None):
    C, hard = build(archive, angular)
    if bare_receipt["hard_depth"] != hard or bare_receipt["cold_modes"] != len(C.q):
        raise ValueError("Mismatched bare-stress receipt.")
    out = []
    q = np.r_[C.q, C.spq]
    for t, m, tag in zip(arms, bare_receipt["arms"], ("native","local-daughter"), strict=True):
        if (t.get("prepared_charged_swap", "native") != tag or m["arm"] != tag
                or t.get("aperture_subdivision", 1) != hard
                or abs(t["energy_defect"]) >= 3e-12 or abs(t["entropy_defect"])>=3e-11):
            raise ValueError("Terminal history did not match its receipt and solver gate.")
        np.testing.assert_array_equal(t["cold_mode_count_history"][0], C.ns)
        np.testing.assert_array_equal(t["cold_momentum_nodes_over_T"], C.q)
        final = np.asarray(t["cold_mode_count_history"][-1])
        if not np.all(np.isfinite(final) & (final > 0)):
            raise ValueError("Nonpositive terminal scalar population.")
        n = np.r_[final, C.spect]
        baseline = check(q, n, C.lam, 0.)
        if abs(baseline["w"]-1/3)>2e-13 or abs(baseline["adiabatic_derivative"]-1/3)>2e-13:
            raise ValueError("Conformal limit did not close.")
        checks = [(str(row["target_pressure_to_energy"]), row["state"])
                  for row in m["thresholds"]]
        caps = ("local_aperture_cap","whole_represented_hard_cap")
        checks.extend((k, m[k]["stress_at_ceiling"]) for k in caps)
        cases = {}
        for k, z in checks:
            now = check(q, n, C.lam, z["bare_mass"])
            if abs(now["w"]-z["pressure_to_energy"])>2e-12:
                raise ValueError("Threshold stress from the earlier independent control did not replay.")
            cases[k] = now
        out.append(dict(arm=tag, baseline=baseline, cases=cases))
    return dict(scope="Homogeneous fixed-comoving-shape, fixed-bare-mass derivative, not perturbative propagation.",
                hard_depth=hard, cold_modes=len(C.q), zero_momentum_fixture=fixture(), arms=out)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for k in ("archive","native","local","bare_receipt"):
        p.add_argument(k)
    p.add_argument("--angular-mesh")
    a = p.parse_args()
    print(json.dumps(run(a.archive, [r(a.native), r(a.local)], r(a.bare_receipt),
                         r(a.angular_mesh) if a.angular_mesh else None), indent=2, allow_nan=False))
