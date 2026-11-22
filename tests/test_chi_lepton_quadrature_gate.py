"""Independent quadrature negative convergence gates: synthetic model only."""
import math
from kinetics.independent_chi_lepton_5d_quadrature import evaluate


def test_independent_thermal_subset_matches_previously_published_reference():
    x=evaluate("fd_cold",nr=4,nmu=4,nc=12,nphi=12,cutoff=12)
    assert math.isclose(x["Qchi"],6.935904667187969e-14,rel_tol=1e-11)
    assert math.isclose(x["entropy"],7.706560741319965e-15,rel_tol=1e-11)
    eq=evaluate("fd_eq_chem",nr=4,nmu=4,nc=12,nphi=12)
    assert abs(eq["Qchi"])<5e-27 and abs(eq["entropy"])<1e-35
    assert eq["max_energy_defect"]<5e-14


def test_oscillatory_radial_and_uv_nonconvergence_including_heat_sign_changes():
    rows=[evaluate("oscillatory",nr=n,nmu=8,nc=24,nphi=16)
          for n in (4,6,8,10)]
    qs=[x["Qchi"] for x in rows]
    assert qs[0]>0>qs[1] and qs[2]>0>qs[3]
    assert all(x["entropy"]>0 and x["negative_entropy_events"]==0 for x in rows)
    cut=[evaluate("oscillatory",nr=8,nmu=8,nc=24,nphi=16,
                  cutoff=c)["Qchi"] for c in (8,12,16)]
    assert cut[0]<0<cut[1] and cut[2]<0


def test_smooth_outgoing_polar_grid_matters_even_when_entropy_positive():
    q8=evaluate("smooth",nr=8,nmu=8,nc=8,nphi=16)["Qchi"]
    q32=evaluate("smooth",nr=8,nmu=8,nc=32,nphi=16)["Qchi"]
    assert abs(q32/q8-1)>.15 and q32>0


def test_energy_and_entropy_are_exact_per_scattering_event():
    for shape in ("smooth","oscillatory"):
        x=evaluate(shape,nr=6,nmu=6,nc=12,nphi=12)
        assert x["Qlepton"]==-x["Qchi"]
        assert x["max_energy_defect"]<1e-13
        assert x["negative_entropy_events"]==0 and x["entropy"]>0
