"""Independent linearisation: FD susceptibility, full nonlinear kernel and Ritz spectrum."""
from __future__ import annotations
import math
import sys
from pathlib import Path
import numpy as np
import pytest
from scipy.special import expit
from scipy.integrate import quad
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"kinetics"))
from linearized_pauli_modes import linearized_modes,fd_logit_basis
from relativistic_nonlinear_fermion import collision_from_nonthermal_chi
from thermal_energy_relaxation import thermal_relaxation_rate


def test_number_orthogonal_fd_basis_and_thermodynamic_norm():
    M=.1;T=1.
    h1,h2,norm,_=fd_logit_basis(M,T)
    def integrate(fun):
        return quad(lambda p:p*p*expit(-math.hypot(p,M))*
          (1-expit(-math.hypot(p,M)))*fun(math.hypot(p,M)),
          0,80,epsrel=1e-11)[0]
    assert abs(integrate(h1))<1e-10
    assert abs(integrate(h2))<1e-10
    assert abs(integrate(lambda E:h1(E)*h2(E)))<1e-10
    for i,h in enumerate((h1,h2)):
        assert math.isclose(integrate(lambda E:h(E)**2),norm[i],rel_tol=1e-11)


def test_positive_selfadjoint_dirichlet_and_nontrivial_distinct_ritz_rates():
    a=linearized_modes(n=16,cutoff=28.)
    b=linearized_modes(n=24,cutoff=32.)
    D=np.array(b["D_matrix_GeV4"])
    S=np.array(b["M_matrix_GeV3"])
    assert np.array_equal(D,D.T)
    assert np.linalg.eigvalsh(D)[0]>0
    assert np.linalg.eigvalsh(S)[0]>0
    lslow,lfast=b["generalized_rates_GeV"]
    gamma=b["thermal_energy_mode_rate_GeV"]
    assert 0<lslow<gamma<lfast
    assert .20<1-lslow/gamma<.40
    for x,y in zip(a["generalized_rates_GeV"],b["generalized_rates_GeV"]):
        assert abs(x/y-1)<5e-4
    assert abs(a["thermal_energy_mode_rate_GeV"]/gamma-1)<5e-4
    assert .40<b["coupling_offdiagonal_fraction"]<.50


def test_independent_full_nonlinear_source_finite_difference_matches_D11():
    M=.1;T=1.;n=12;cut=24.
    r=linearized_modes(n=n,cutoff=cut)
    h1,_,_,_=fd_logit_basis(M,T)
    eps=.001
    def Q(sign):
        f=lambda p:expit(-np.hypot(p,M)/T+
                        sign*eps*h1(np.hypot(p,M)))
        return collision_from_nonthermal_chi(
          f,M=M,Te=T,G=.001,n_p=n,n_k=n,n_mu=n,n_cm=n,
          n_phi=n,cutoff=cut,p_max=cut)["Q_chi_GeV5"]
    numeric=(Q(+1)-Q(-1))/(2*eps)
    assert math.isclose(numeric,-T*r["D_matrix_GeV4"][0][0],
                        rel_tol=1e-5)


def test_thermal_FD_moment_gamma_equals_rayleigh_quotient_not_slowest_rate():
    g=thermal_relaxation_rate(n=24,cutoff=32.,delta_T=.01)[
       "gamma_T_energy_moment_GeV"]
    r=linearized_modes(n=24,cutoff=32.)
    assert math.isclose(g,r["thermal_energy_mode_rate_GeV"],rel_tol=2e-5)
    assert abs(g/r["generalized_rates_GeV"][0]-1)>.2


def test_zero_portal_and_bad_parameters():
    assert linearized_modes(G=0.,n=8)["generalized_rates_GeV"]==[0.,0.]
    with pytest.raises(ValueError):linearized_modes(M=0.)
    with pytest.raises(ValueError):linearized_modes(n=2)
