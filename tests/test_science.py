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


def test_power_law_open_diffusion_selection():
    solver=load("continuation-06/code/generalized_escape_flux.py")
    result=solver.run()
    assert result["checks_passed"] == 12
    for case in result["scenarios"]:
        assert abs(case["predicted_edge_power"]-case["numerical_log_slope"]) < .03


def test_CLASS_file_format_and_analytic_moments():
    import tempfile
    from pathlib import Path
    c = load("cosmology/thermal_bose_class.py")
    with tempfile.TemporaryDirectory() as folder:
        file = Path(folder)/"be-numeric-only.dat"
        test = c.write_psd(file)
        lines = file.read_text().splitlines()
        assert len(lines)>2000
        assert len(lines[0].split())==2  # CLASS fscanf cannot read an ASCII header
        float(lines[0].split()[0]);float(lines[0].split()[1])
        assert abs(test["relative_number_moment_error"])<1e-4
        assert abs(test["relative_radiation_moment_error"])<1e-4
    b = c.thermal_background()
    assert abs(b["omega_h2"]-.11912939486124709)<1e-6

def test_number_preserving_toy_cold_PSD_analytic():
    import tempfile
    from pathlib import Path
    toy = load("cosmology/cold_psd_scan.py")
    with tempfile.TemporaryDirectory() as folder:
        for alpha in (1.,.5,.3,.1):
            value=toy.cold_psd_file(Path(folder)/("alpha"+str(alpha)+".dat"),alpha)
            assert abs(value["number_ratio"]-1)<1e-4
            assert abs(value["early_radiation_ratio"]/alpha-1)<1e-4


def test_bimodal_bose_population_moments():
    import tempfile
    from pathlib import Path
    mod=load("cosmology/burst_residual_hot.py")
    with tempfile.TemporaryDirectory() as folder:
        for alpha in (.02,.01):
            for w in (0.,.01,.05,.1,.2):
                data=mod.write_bimodal(Path(folder)/("b"+str(alpha)+"-"+str(w)+".txt"),alpha,w)
                assert abs(data["number_moment_ratio"]-1)<1e-4
                assert abs(data["radiation_moment_ratio"]-1)<1e-4
                expect=w/((1-w)*alpha+w)
                assert abs(data["hot_fraction_of_early_relativistic_energy"]-expect)<1e-12


def test_quantum_1to2_event_level_kinetic_conservation():
    import numpy as np
    model = load("kinetics/one_to_two_fixed_bath.py")
    k = model.FlatThermalDecay(nq=50,np1=64,nmu=48)
    balance,energy=k.detailed_balance_error()
    assert balance<1e-10 and energy<1e-12
    rp=k.avg_q[k.ip]*(1-k.wp)+k.avg_q[k.ip+1]*k.wp
    rc=k.avg_q[k.ic]*(1-k.wc)+k.avg_q[k.ic+1]*k.wc
    assert max(np.max(np.abs(rp-k.qphi)),np.max(np.abs(rc-k.qchi))) <1e-12
    zero=np.zeros(2*k.nq)
    assert np.min(k.rhs(0.,zero,True))>=-1e-12
    saturated=zero.copy()
    saturated[:k.nq]=1.
    saturated[k.nq:]=1.
    assert np.max(k.rhs(0.,saturated,True)[k.nq:])<=1e-12

def test_CLASS_two_population_is_one_scalar_mass_invariant():
    import tempfile
    from pathlib import Path
    model = load("cosmology/burst_residual_hot.py")
    b = model.thermal_background()
    with tempfile.TemporaryDirectory() as d:
        cfg = model.class_mixture_params(b,.02,.05,d)
        assert cfg["N_ncdm"] == 2
        assert cfg["m_ncdm"] == "86,86"
        assert cfg["deg_ncdm"] == "1,1"
        assert cfg["use_ncdm_psd_files"] == "1,1"
        for file in cfg["ncdm_psd_filenames"].split(","):
            rows=Path(file).read_text().splitlines()
            assert len(rows)>2000 and len(rows[0].split())==2
        same = model.class_mixture_params(b,.02,0,d)
        assert same["N_ncdm"] == 1
        assert abs(same["T_ncdm"]/.02/b["T_ratio"]-1)<1e-12
