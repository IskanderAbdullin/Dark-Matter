"""Five physical FD+BE fields, FOUR Yukawa reactions: finite-step moment checks.

Uses REAL massive species DOS quadrature and common positive on-shell
1↔2/2↔2 event weights. Previous two positive 1↔2 counts must NOT be
assumed conserved after physically allowed pair annihilation.
"""
import math
import numpy as np
from scipy.sparse import issparse

from kinetics.phi_bose_five_dirac_full_yukawa_weak import FullDiracYukawaWeak
from kinetics.phi_bose_five_dirac_full_yukawa_implicit import (
    FullDiracPhysicalMoments,asflat,unflat,
)


def _m():
    return FullDiracPhysicalMoments(
        P=2.,Ppair=2.,N1=3,N2=3,Nphi=8,nq=18,
        nrad=4,nangle=8,pair_nrad=4,
        pair_nrelative=8,pair_ncm=32,pair_naz=4)


def test_physical_DOS_capacity_energy_bernstein_graph_exact_four_channels():
    m=_m()
    assert m.nvars==sum(s.nvars for s in m.species)
    assert len(m.channels)==4
    assert m.total_events==sum(len(c.W) for c in m.channels)
    assert all(np.min(c.W)>0 for c in m.channels)
    for ch in m.channels:
        assert ch.max_energy_defect<1e-12
        assert np.max(abs(ch.S@m.charge_nodal))<5e-12
        assert np.max(abs(ch.S@m.energy_nodal))<5e-12
    for s in m.species:
        p0=np.sqrt(np.maximum(0.,(s.edges[:-1]-s.mass)*
                                        (s.edges[:-1]+s.mass)))
        p1=np.sqrt(np.maximum(0.,(s.edges[1:]-s.mass)*
                                        (s.edges[1:]+s.mass)))
        V=s.g*(p1**3-p0**3)/(6*math.pi**2)
        got=np.sum(s.weights,axis=1)
        assert np.max(abs(got/V-1))<5e-13


def test_five_field_weak_collision_at_thermal_state_matches_independent_full_Dirac_weak_source():
    m=_m()
    T=(.9,.9,1.,1.,1.)
    start=m.initial_etas(T=T)
    full=m.source(start)
    oracle=m.weak.source(T=T)
    assert full["H"]>0 and oracle["H_GeV4"]>0
    for channel in full["channels"]:
        name=channel["name"]
        d=next(x for x in oracle["event_channels"].items() if x[0]==name)[1]
        assert math.isclose(channel["event_rate"],d["event_rate_GeV4"],
                            rel_tol=1e-10,abs_tol=1e-24)
    assert math.isclose(full["H"],oracle["H_GeV4"],rel_tol=1e-10)
    assert abs(full["exact_event_U1"])<1e-21
    assert abs(full["exact_event_energy"])<1e-20


def test_full_chemical_U1_equilibrium_including_both_pair_channels_stationary():
    m=_m()
    eq=m.initial_etas(T=(1.,)*5,mu=(.05,-.05,.05,-.05,0.))
    src=m.source(eq,with_jac=True)
    assert abs(src["H"])<1e-27
    assert max(abs(x["event_rate"]) for x in src["channels"])<1e-23
    assert issparse(src["dC"])
    assert np.max(abs(src["C"]))<1e-23
    out=m.implicit_step(eq,10000.)
    assert out["iterations"]==0 and out["stalled_at_precision"]
    assert max(np.max(abs(a-b)) for a,b in zip(eq,out["etas"]))==0.


def test_five_field_full_physical_implicit_finite_time_charge_energy_entropy_and_Bose_bound():
    m=_m()
    initial=m.initial_etas()
    before=m.integrate_spectrum(initial)
    o=m.adaptive_step(initial,10000.)
    assert 0<o["accepted_dt"]<=10000.
    assert o["iterations"]>0
    assert o["scaled_residual"]<2e-7
    assert o["min_Bose_eta"]>0
    assert o["physical_H"]>=0
    assert o["entropy_increase"]>0
    assert o["entropy_certificate"]>=-1e-10
    after=m.integrate_spectrum(o["etas"])
    assert abs(after["U1"]-before["U1"])<3e-11
    assert abs(sum(after["Q"]-before["Q"]))<3e-11
    for spec,eta in zip(m.species,o["etas"]):
        assert np.min(spec.evaluate(eta)["M"])>=0
        assert np.isfinite(spec.evaluate(eta,with_entropy=True)["S"])


def test_sparse_full_five_species_jacobian_covers_decay_and_pair_cross_derivatives():
    m=_m()
    old=asflat(m.initial_etas())
    new=old*(1+.002*np.cos(np.arange(len(old))*.1))
    dt=5000.
    r,j,_=m.implicit_residual(new,old,dt)
    assert issparse(j)
    assert j.shape==(len(old),len(old))
    assert np.all(np.isfinite(j.data))
    direction=old*np.cos(.21*np.arange(len(old))+.5)
    delta=1e-5
    rp,_=m.implicit_residual(new+delta*direction,old,dt,with_jac=False)
    rm,_=m.implicit_residual(new-delta*direction,old,dt,with_jac=False)
    numeric=(rp-rm)/(2*delta)
    analytic=j@direction
    discrepancy=np.linalg.norm(numeric-analytic)/np.linalg.norm(analytic)
    assert discrepancy<6e-5,discrepancy
    # Different particle-antiparticle blocks are COUPLED through pair
    # channels: isolated decay-only Jacobian would have no these blocks.
    chi1_anti=j[m.offsets[0]:m.offsets[1],m.offsets[1]:m.offsets[2]]
    chi2_anti=j[m.offsets[2]:m.offsets[3],m.offsets[3]:m.offsets[4]]
    assert chi1_anti.nnz>0 and chi2_anti.nnz>0
    assert max(abs(chi1_anti.data))>0 and max(abs(chi2_anti.data))>0


def test_nonzero_mu_restricted_old_three_species_fixed_point_violates_full_pair_equilibrium():
    m=_m()
    e=m.initial_etas(T=(1.,)*5,mu=(-.15,.15,.05,-.05,-.2))
    src=m.source(e)
    chan={x["name"]:x for x in src["channels"]}
    assert abs(chan["chi1_decay"]["event_rate"])<1e-23
    assert chan["anti_chi1_decay"]["event_rate"]>0
    assert chan["chi1_pair"]["event_rate"]>0
    assert chan["chi2_pair"]["event_rate"]>0
    assert src["H"]>0
    assert abs(src["exact_event_U1"])<1e-21
    assert abs(src["exact_event_energy"])<1e-20


def test_true_massive_pair_IR_grid_splits_at_all_exact_hard_parent_support_edges():
    from kinetics.phi_bose_pair_ir_parent_support import global_pair_scalar_support_floor
    from kinetics.phi_bose_physical_cutcell_weak import M1,M2,MPHI
    m=_m()
    e=m.species[4].edges
    for mass in (M1,M2):
        cutoff=global_pair_scalar_support_floor(m.Ppair,mass)
        assert np.min(abs(e-cutoff))<2e-14
    assert math.isclose(e[0],MPHI,rel_tol=1e-15)
    kmin=global_pair_scalar_support_floor(m.Ppair,M2)
    assert math.isclose(e[1],kmin,rel_tol=2e-13)
    src=m.source(m.initial_etas())
    # Below global pair2 kmin no physical full-network event reaches phi.
    assert np.max(abs(src["C"][m.offsets[4]:m.offsets[4]+3]))==0.


def test_float64_far_tail_root_precision_is_reported_not_misclassified_as_no_solution():
    m=_m()
    old=m.initial_etas()
    out=m.adaptive_step(old,10000.)
    assert out["iterations"]>0
    assert out["scaled_residual"]<2e-7
    assert out["max_abs_moment_residual"]<1e-15
    assert out["raw_scaled_residual"]>=out["scaled_residual"]
    assert 0<=out["precision_limited_dofs"]<=m.nvars
    assert abs(out["energy_defect"])<3e-11
    assert abs(out["U1_defect"])<3e-11
    assert out["entropy_certificate"]>=-1e-10
    # CI previously requested eta correction <1 ulp in a far-tail mode.
    far_eta=4.086731509081753
    impossible_delta=1.1567543973059635e-16
    assert far_eta+impossible_delta==far_eta
    assert np.nextafter(far_eta,np.inf)-far_eta>impossible_delta


def test_near_equilibrium_quantum_affinity_stable_flux_and_event_H_positive():
    m=_m()
    for mu_phi in (-1e-12,+1e-12):
        etas=m.initial_etas(T=(1.,)*5,
               mu=(.05,-.05,.05,-.05,mu_phi))
        flat=asflat(etas)
        src=m.source(etas)
        assert src["H"]>=0
        for channel in m.channels:
            result=channel.weights_and_derivatives(flat)
            affine=np.asarray(channel.S@flat)
            assert np.min(result["J"]*affine)>=-1e-26
            # F-R may lose relative accuracy by cancellation; its absolute
            # distance from stable expm1 flux remains at float64 level.
            old_diff=result["F"]-result["R"]
            err=np.max(abs(result["J"]-old_diff)/
                       np.maximum(result["F"]+result["R"],1e-300))
            assert err<4e-15
            assert np.all(result["J"]*affine>=-1e-26)
        assert abs(src["exact_event_U1"])<1e-21
        assert abs(src["exact_event_energy"])<1e-20


def test_full_dirac_quantum_implicit_across_three_physical_parent_cutoffs():
    from kinetics.phi_bose_five_dirac_kl_certificate import quantum_KL_certificate
    for P in (2.,8.,16.):
        m=FullDiracPhysicalMoments(
            P=P,Ppair=P,N1=3,N2=3,Nphi=8,nq=18,
            nrad=4,nangle=8,pair_nrad=4,
            pair_nrelative=8,pair_ncm=32,pair_naz=4)
        initial=m.initial_etas()
        all_channels=m.source(initial)
        assert len(all_channels["channels"])==4
        assert all_channels["H"]>0
        assert abs(all_channels["exact_event_U1"])<1e-18
        assert abs(all_channels["exact_event_energy"])<1e-17
        result=m.adaptive_step(initial,5000.)
        assert 0<result["accepted_dt"]<=5000.
        assert result["iterations"]>0
        assert result["scaled_residual"]<2e-7
        assert result["min_Bose_eta"]>0
        assert result["entropy_increase"]>0
        assert result["entropy_certificate"]>=-1e-10
        certificate=quantum_KL_certificate(
            m,initial,result["etas"],result["accepted_dt"])
        assert all(v>=-1e-11 for v in
                   certificate["physical_quantum_KL_by_species"])
        assert abs(certificate["exact_entropy_decomposition_defect"])<3e-10
        assert abs(result["energy_defect"])<3e-11
        assert abs(result["U1_defect"])<3e-11
        assert result["max_abs_moment_residual"]<3e-14
