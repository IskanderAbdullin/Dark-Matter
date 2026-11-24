"""Independent physical quantum KL identity for coupled backward Euler.

Checks the entropy inequality is a thermodynamic consequence of the
*same physical DOS and event manifold*, not an assertion from successful CI.
"""
import math
import numpy as np
import pytest

from kinetics.phi_bose_coupled_bernstein_p2_finite import Coupled
from kinetics.phi_bose_coupled_bernstein_p2_implicit import implicit_step
from kinetics.phi_bose_implicit_entropy_kl import (
    scalar_quantum_kl,physical_relative_entropy,entropy_residual_audit,
)


def test_bernoulli_and_geometric_KL_nonnegative_over_independent_states():
    for kind in ("fermion","boson"):
        old=np.array([.00001,.01,.15,.6,1.,5.,15.])
        new=old*(1+.4*np.sin(np.arange(len(old))+.3))
        values=scalar_quantum_kl(old,new,kind)
        assert np.all(values>=-3e-13)
        assert max(values)>1e-4
        assert np.max(abs(scalar_quantum_kl(old,old,kind)))<1e-13
    with pytest.raises(ValueError):
        scalar_quantum_kl([.1,-.1],[.2,.2],"boson")


def test_physical_entropy_gap_matches_KL_plus_exact_implicit_residual_work():
    for P in (2.,8.,16.):
        model=Coupled(P=P,Nphi=8,nrad=4,nangle=8,nq=24)
        old=model.initial_etas()
        trial=[a*(1+.08*np.sin(np.arange(len(a))+.2)) for a in old]
        decomposition=entropy_residual_audit(model,old,trial,1e5)
        kl=decomposition["integrated_quantum_kl"]
        assert len(kl)==3
        assert all(v>=-1e-12 for v in kl)
        assert sum(kl)>1e-9
        assert abs(decomposition["decomposition_defect"])<3e-12
        assert math.isclose(
            decomposition["entropy_certificate"],
            decomposition["integrated_quantum_kl"][0]
            +decomposition["integrated_quantum_kl"][1]
            +decomposition["integrated_quantum_kl"][2]
            +decomposition["residual_work"],
            abs_tol=3e-12,rel_tol=5e-9
        )


def test_implicit_root_has_positive_KL_and_stronger_entropy_certificate():
    for P in (2.,8.,16.):
        model=Coupled(P=P,Nphi=12,nrad=6,nangle=12,nq=24)
        old=model.initial_etas()
        step=implicit_step(model,old,1e5)
        result=entropy_residual_audit(model,old,step["etas"],1e5)
        assert step["iterations"]>0
        assert result["new_microscopic_H"]>0
        assert result["actual_entropy_increase"]>0
        assert all(v>=-2e-12 for v in result["integrated_quantum_kl"])
        assert sum(result["integrated_quantum_kl"])>1e-12
        assert result["entropy_certificate"]>0
        assert abs(result["decomposition_defect"])<3e-12
        assert math.isclose(
            result["entropy_certificate"],
            step["entropy_certificate"],
            rel_tol=8e-7,abs_tol=3e-12
        )


def test_stationary_nonzero_mu_zero_KL_and_zero_entropy():
    model=Coupled(P=2.,Nphi=12,nrad=6,nangle=12)
    old=model.initial_etas(T=(1.,1.,1.),mu=(-.15,.05,-.2))
    step=implicit_step(model,old,1e5)
    r=entropy_residual_audit(model,old,step["etas"],1e5)
    assert abs(r["decomposition_defect"])<1e-16
    assert max(abs(v) for v in r["integrated_quantum_kl"])<1e-18
    assert abs(r["new_microscopic_H"])<1e-30
