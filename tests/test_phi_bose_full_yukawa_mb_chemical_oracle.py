"""Independent MB oracle for the real-scalar four-Yukawa linear chemistry."""
import math
import numpy as np
from kinetics.phi_bose_full_yukawa_mb_chemical_oracle import (
    ideal_MB_number,equilibrium_MB_four_channel_rates,
    chemical_MB_eigenmodes
)
from kinetics.phi_bose_yukawa_charge_audit import stoichiometry

def test_MB_density_massless_analytic_and_separate_Dirac_spin_count():
    for T in (.1,1.,3.):
        assert math.isclose(ideal_MB_number(T,0.,1),
                            T**3/math.pi**2,rel_tol=2e-15)
        assert math.isclose(ideal_MB_number(T,0.,2),
                            2*T**3/math.pi**2,rel_tol=2e-15)
        assert ideal_MB_number(T,.1,2)<ideal_MB_number(T,0.,2)

def test_independent_MB_four_channel_rates_match_physical_normalization():
    v=equilibrium_MB_four_channel_rates(1.,128)
    assert v.shape==(4,) and np.min(v)>0
    assert v[0]==v[1]
    for v_actual,v_expect in zip(v,(2.7458059471132872e-11,
                                    2.7458059471132872e-11,
                                    3.4912014442108976e-16,
                                    2.0886349448791195e-16)):
        assert math.isclose(v_actual,v_expect,rel_tol=3e-5)

def test_chemical_relaxation_has_exact_one_vector_U1_zero_and_four_positive_modes():
    for T,g in ((.1,17.25),(1.,61.75)):
        a=chemical_MB_eigenmodes(T,g,96)
        lam=np.asarray(a["chemical_eigenvalues_GeV"])
        assert lam.shape==(5,)
        assert abs(lam[0]/lam[-1])<1e-12
        assert np.min(lam[1:])>0
        assert a["true_Dirac_U1_relative_matrix_zero_defect"]<1e-12
        S=stoichiometry(True)
        assert np.array_equal(S@np.array([1,-1,1,-1,0]),np.zeros(4))

def test_T1_MB_smallest_nonzero_chemical_mode_about_463_H_not_a_relic_proof():
    a=chemical_MB_eigenmodes(1.,61.75,128)
    eig=a["chemical_eigenvalues_over_H"]
    assert math.isclose(eig[1],462.952579,rel_tol=2e-5)
    assert math.isclose(eig[2],8979.28208,rel_tol=2e-5)
    assert math.isclose(eig[3],2.53971134e8,rel_tol=2e-5)
    assert math.isclose(eig[4],7.61273740e8,rel_tol=2e-5)
    assert a["class_ready"] is False
