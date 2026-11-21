"""Deterministic mathematics tests; passing does not verify astrophysical viability."""
from __future__ import annotations
import ast
import importlib.util
import math
import sys
from pathlib import Path
import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]

def load(relative):
    path = ROOT / relative
    sys.path.insert(0, str(path.parent))
    try:
        spec = importlib.util.spec_from_file_location(path.stem, path)
        assert spec is not None and spec.loader is not None
        mod = importlib.util.module_from_spec(spec)
        sys.modules[path.stem] = mod
        spec.loader.exec_module(mod)
        return mod
    finally:
        sys.path.pop(0)

@pytest.fixture(scope="module")
def ladder():
    return load("continuation-05/reproduce_cutoff_ladder.py")

def test_python_source_syntax():
    scripts = list(ROOT.rglob("*.py"))
    assert len(scripts) >= 25, "Missing research programs from migration"
    for src in scripts:
        ast.parse(src.read_text(encoding="utf-8"), filename=str(src))

def test_positive_bose_taylor_remainders(ladder):
    t = np.geomspace(1e-5, 3., 180)
    for k in range(1, 7):
        f = ladder.occupation(t, 4., k)
        assert np.all(np.isfinite(f))
        assert np.all(f > 0), f"nonpositive occupation k={k}"
        assert np.all(np.diff(f) > 0), f"not monotone k={k}"

def test_lowered_bose_edge_exponent(ladder):
    depths = np.array([.004, .008, .016, .032])
    for k in range(1, 5):
        rho = ladder.I(depths, 4., k)
        slope = np.polyfit(np.log(depths), np.log(rho), 1)[0]
        assert abs(slope - (k + 1.5)) < .035, (k, slope)

def test_self_consistent_subcritical_profiles(ladder):
    for k, radius, mass in ((1, 27.37249, 18.75118), (2, 78.42542, 17.90742)):
        ans = ladder.halo(k, psi0=3.8)
        assert ans["finite_edge_found"]
        assert abs(ans["radius_over_R_or_limit"] / radius - 1) < .002
        assert abs(ans["Q_at_end"] / mass - 1) < .002

def test_critical_lowered_bose_profile():
    theory = load("continuation-04/code/critical_lowered.py")
    Ic, rt, Q, _, _ = theory.solve(eta=4., nq=96)
    assert abs(Ic - 1.087257505) < .001
    assert abs(rt - 40.29879) < .15
    assert abs(Q - 26.35702) < .15
    c = 1 / np.expm1(4.)
    assert theory.casimir_second(0., 4.) > 0
    assert math.isclose(theory.casimir_second(1.,4.), 1/((1+c)*(2+c)), rel_tol=1e-12)

def test_quantum_detailed_balance_contrast(ladder):
    t = np.array([.05, .35, .20, .20])
    assert abs(t[0]+t[1]-t[2]-t[3]) < 1e-14
    thermal = 1 / np.expm1(4.-t)
    L = np.log1p(1/thermal)
    assert abs(L[0]+L[1]-L[2]-L[3]) < 1e-12
    for k in range(1, 7):
        f = ladder.occupation(t, 4., k)
        L = np.log1p(1/f)
        assert L[0]+L[1]-L[2]-L[3] > .5

def test_degeneracy_near_critical():
    import mpmath as mp
    mp.mp.dps=35
    nuc=mp.zeta(mp.mpf("1.5"))
    z=mp.exp(mp.mpf("-1e-8"))
    nu=mp.polylog(mp.mpf("1.5"),z)
    A=mp.polylog(mp.mpf(".5"),z)/nu
    assert A>1000
    assert abs(float((nuc-nu)*A-2*mp.pi/nuc)) < .002

def test_open_flux_exact_riccati():
    solver = load("continuation-05/reproduce_open_flux.py")
    answer = solver.run()
    assert answer["tests_passed"] >= 4
    assert math.isclose(answer["exact_finite_binding_pole_s"],4.67940655,abs_tol=1e-7)


def test_finite_escape_state_density():
    model=load("continuation-06/code/phase_density_edge.py")
    result=model.run()
    assert 18.5 < result["ratio_g0_quadratic_to_linear"] < 19.6
    assert abs(result["ratio_mass_quadratic_to_linear"]-1.) < .01
    assert all(x["g0"]>0 and x["gp0"]<0 for x in result["models"])

def test_shallow_cutoff_state_density_scaling():
    model=load("continuation-06/code/dos_scaling.py")
    result=model.run()
    assert abs(result["G0"] - 157.32910288) < .05
    assert abs(result["points"][-1]["relative_error"]) < .001
