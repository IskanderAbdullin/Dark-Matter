"""Independent closed ν+νbar+Zprime chemical equilibrium and charges."""
import math
import pytest
from kinetics.vector_neutrino_equilibrium_charges import (
 equilibrium_from_empty_vector,independent_massless_series_equilibrium)

@pytest.mark.parametrize("m",[.03,.1,.3])
def test_closed_quantum_equilibrium_conserves_Nnu_plus_2NV_and_total_energy(m):
    a=equilibrium_from_empty_vector(m,n=128)
    b=equilibrium_from_empty_vector(m,n=256)
    assert abs(a["relative_N_residual"])<1e-11
    assert abs(a["relative_E_residual"])<1e-11
    assert a["S_total_final_over_initial"]>1
    assert a["mu_nu_over_T_final"]<0
    assert 0<a["n_V_final_over_thermal_mu0_initial"]<.2
    assert math.isclose(a["T_final_over_T0"],b["T_final_over_T0"],rel_tol=1e-8)
    assert math.isclose(a["mu_nu_over_T_final"],b["mu_nu_over_T_final"],rel_tol=1e-8)

def test_independent_massless_polylog_fugacity_limit():
    a=independent_massless_series_equilibrium()
    b=equilibrium_from_empty_vector(.03)
    assert math.isclose(a["T_final_over_T0"],b["T_final_over_T0"],rel_tol=4e-5)
    assert math.isclose(a["mu_nu_over_T_final"],b["mu_nu_over_T_final"],rel_tol=2e-4)
    assert math.isclose(a["nV_over_initial_nnu"],b["nV_over_initial_nnu"],rel_tol=2e-3)
    assert abs(a["relative_number_error"])<1e-12
    assert abs(a["relative_energy_error"])<1e-12
