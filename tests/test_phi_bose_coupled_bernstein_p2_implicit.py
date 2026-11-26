"""Independent analytical/numerical audits of the fully implicit physical Yukawa
Bernstein-P2 three-species moment solver.

The tested model is a closed *illustrative dark-sector* collision operator
at fixed physical particle masses, fixed vacuum Yukawa width and parent P.
It is NOT finite-SM+FRW, 86-eV relic freeze-out, a condensate or CLASS.
"""
import math

import numpy as np
from scipy.sparse import issparse

from kinetics.phi_bose_coupled_bernstein_p2_finite import Coupled
from kinetics.phi_bose_coupled_bernstein_p2_implicit import (
    asflat, unflat, physical_implicit_residual, implicit_step,
    implicit_adaptive_step, integrate_to,
)


def _state(model, variant):
    eta = model.initial_etas()
    if variant == "nonthermal":
        eta[0] = eta[0]*(1+.021*np.sin(np.arange(len(eta[0]))+.3))
        eta[1] = eta[1]*(1+.016*np.cos(np.arange(len(eta[1]))+.1))
        eta[2] = eta[2]*(1+.018*np.sin(np.arange(len(eta[2]))+.7))
    return eta


def test_sparse_full_species_jacobian_matches_independent_centered_variation():
    m = Coupled(P=2., Nphi=8, nrad=4, nangle=8, nq=24)
    old = asflat(_state(m, "nonthermal"))
    new = old*(1+.022*np.sin(np.arange(len(old))+.1))
    dt = 75000.
    r, J, info = physical_implicit_residual(m, new, old, dt)
    assert issparse(J)
    assert J.shape == (len(old), len(old))
    assert info["H"] >= 0 and info["weights_positive"]
    assert np.all(np.isfinite(r))
    assert np.all(np.isfinite(J.data))
    direction = new*np.cos(np.arange(len(old))*.37+.2)
    eps=3.e-5
    rp,_ = physical_implicit_residual(
        m,new+eps*direction,old,dt,with_jac=False)
    rm,_ = physical_implicit_residual(
        m,new-eps*direction,old,dt,with_jac=False)
    independent=(rp-rm)/(2*eps)
    analytic=J@direction
    relative=np.linalg.norm(independent-analytic)/np.linalg.norm(analytic)
    assert relative<6.e-6,relative
    # The Jacobian MUST couple all three fields; a diagonal-only
    # approximation would not be the physical implicit system.
    for i in range(3):
        for j in range(3):
            if i==j:continue
            sub=J[m.offsets[i]:m.offsets[i+1],
                  m.offsets[j]:m.offsets[j+1]]
            assert sub.nnz>0 and np.max(abs(sub.data))>0


def test_implicit_three_species_positive_entropy_and_exact_physical_charges():
    for P in (2.,8.,16.):
        m=Coupled(P=P,Nphi=12,nrad=6,nangle=12,nq=28)
        for variant in ("thermal","nonthermal"):
            initial=_state(m,variant)
            before=m.integrate_spectrum(initial)
            step=implicit_adaptive_step(m,initial,100000.)
            after=m.integrate_spectrum(step["etas"])
            assert step["accepted_dt"]==100000.
            assert step["min_bose_eta"]>0
            assert step["iterations"]>0
            assert step["scaled_residual"]<5e-8
            assert step["entropy_increase"]>0
            assert step["microscopic_H_new"]>=0
            scale=max(1e-25,100000.*abs(m.source(initial)["Rphi"]))
            assert step["entropy_certificate"] >= -1e-9*scale
            assert abs(after["S"]-before["S"]-step["entropy_increase"])<1e-11
            assert abs(step["energy_defect"])<1e-8*scale
            assert max(abs(x) for x in step["number_residuals"])<1e-8*scale
            for sp,eta in zip(m.species,step["etas"]):
                assert np.all(np.isfinite(eta))
                occ=sp.evaluate(eta,with_entropy=True)
                assert np.min(occ["M"])>=0
                assert math.isfinite(occ["S"])
            # Distinguishes continuum H>=0 (insufficient for finite Euler)
            # from physical implicit entropy concavity certificate.
            assert step["entropy_increase"] >= (
                step["accepted_dt"]*step["microscopic_H_new"]
                -1e-9*scale
            )


def test_implicit_outperforms_explicit_trial_that_failed_positive_bose_inverse():
    m=Coupled(P=2.,Nphi=16,nrad=8,nangle=12,nq=28)
    initial=m.initial_etas()
    try:
        m.step(initial,300000.)
    except (ValueError,ArithmeticError):
        pass
    else:
        raise AssertionError("Explicit failed trial regression changed; review whether issue was fixed")
    step=implicit_step(m,initial,300000.)
    assert step["accepted_dt"]==300000.
    assert step["iterations"]>0
    assert step["scaled_residual"]<5e-8
    assert step["entropy_increase"]>0
    assert step["entropy_certificate"]>=-1e-11
    assert step["min_bose_eta"]>0
    assert abs(step["energy_defect"])<1e-13
    assert max(abs(x) for x in step["number_residuals"])<1e-13


def test_common_temperature_nonzero_mu_is_stationary_in_implicit_equations():
    for P in (2.,8.,16.):
        m=Coupled(P=P,Nphi=12,nrad=6,nangle=12)
        old=m.initial_etas(T=(1.,1.,1.),mu=(-.15,.05,-.2))
        flat=asflat(old)
        residual,meta=physical_implicit_residual(m,flat,flat,100000.,with_jac=False)
        # Exact F=R, including physically nonzero chemical potentials.
        assert np.max(abs(residual))<1e-17
        assert abs(meta["H"])<1e-30
        r=implicit_step(m,old,100000.)
        assert r["iterations"]==0
        assert r["stalled_at_precision"]
        assert max(np.max(abs(a-b)) for a,b in zip(r["etas"],old))==0


def test_first_order_backward_euler_temporal_refinement_in_actual_three_species_solver():
    m=Coupled(P=2.,Nphi=8,nrad=4,nangle=8,nq=24)
    initial=_state(m,"nonthermal")
    baseline=m.integrate_spectrum(initial)["N"][2]
    shifts=[]
    for dt in (75000.,37500.,18750.):
        out=integrate_to(m,initial,300000.,dt)
        assert out["rejected"]==0
        end=m.integrate_spectrum(out["etas"])
        shifts.append(end["N"][2]-baseline)
        assert all(v["entropy_certificate"]>=-1e-12 for v in out["history"])
    e1=abs(shifts[0]-shifts[1])
    e2=abs(shifts[1]-shifts[2])
    assert e1>1e-10 and e2>1e-10
    ratio=e1/e2
    assert 1.7<ratio<2.5,(shifts,ratio)
