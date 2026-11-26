"""Independent MB-Gondolo-Gelmini and exact quantum pair on-shell H checks."""
import math
import numpy as np
from numpy.polynomial.legendre import leggauss
from kinetics.phi_bose_pair_independent_vectorized import (
    sigma_independent,reaction_density_MB
)
from kinetics.phi_bose_pair_physical_onshell_quantum import (
    differential_pair_sigma,PhysicalPairEvents
)

def test_CM_differential_matrix_includes_one_real_scalar_identical_half():
    u,w=leggauss(96);u=14*u;c=np.tanh(u);weights=14*w*(1-c*c)
    for sqrt_s in (.05,.1,.5,1.,5.):
        s=sqrt_s**2
        actual=np.dot(weights,differential_pair_sigma(s,c))
        independent=sigma_independent(s,96)
        assert actual>0 and np.all(differential_pair_sigma(s,c)>=0)
        assert math.isclose(actual,independent,rel_tol=5e-11)

def test_positive_physical_real_mass_lab_events_and_exact_four_leg_energy():
    m=PhysicalPairEvents(P=20.,nrad=12,nrelative=12,ncm=32,naz=4)
    assert m.event_count==12**2*12*32*4
    assert min(m.W)>0
    assert min(m.e_phi3)>86e-9 and min(m.e_phi4)>86e-9
    assert m.max_kinematic_energy_defect<1e-13
    assert np.max(np.abs(m.e_phi3+m.e_phi4-m.e_chi-m.e_anti))<1e-13

def test_lab_5D_MB_weak_normalization_is_independent_from_Gondolo_Gelmini():
    g_ref=reaction_density_MB(1.,128)["gamma_forward_GeV4"]
    m=PhysicalPairEvents(P=20.,nrad=24,nrelative=32,ncm=96,naz=1)
    val=m.MB_forward(1.)
    assert math.isclose(val,g_ref,rel_tol=3e-4)
    assert math.isclose(val,2.0886565136041006e-16,rel_tol=2e-5)
    # Parent/cos grid is finite, so not a pointwise continuum theorem.

def test_pair_quantum_thermal_chemical_equilibrium_has_zero_net_source():
    m=PhysicalPairEvents(P=20.,nrad=12,nrelative=12,ncm=32,naz=4)
    g=m.MB_forward(1.)
    r=m.quantum_weak(Tchi=1.,Tanti=1.,Tphi=1.,
                     muchi=.05,muanti=-.05,muphi=0.)
    assert abs(r["R_events_GeV4"])<1e-14*g
    assert abs(r["quantum_H_GeV4"])<1e-15*g
    assert r["min_Bose_eta"]>0
    assert abs(r["total_energy_defect"])<1e-28
    assert r["conserved_U1_defect"]==0.

def test_pair_quantum_nonthermal_and_scalar_underpopulation_have_positive_H():
    m=PhysicalPairEvents(P=20.,nrad=16,nrelative=12,ncm=48,naz=8)
    for params in (dict(muchi=.05,muanti=-.05,muphi=-.2),
                   dict(Tchi=.9,Tanti=1.1,Tphi=1.,bose_curvature=.01)):
        s=m.quantum_weak(**params)
        assert s["quantum_H_GeV4"]>0
        assert s["min_micro_entropy_JA"]>=-1e-12
        assert abs(s["total_energy_defect"])<1e-28
        assert s["conserved_U1_defect"]==0.
        assert s["dN_phi"]+2*s["dN_chi"]==0.
    d=m.quantum_weak(muchi=.05,muanti=-.05,muphi=-.2)
    assert math.isclose(d["quantum_H_GeV4"],
                        .4*d["R_events_GeV4"],rel_tol=2e-13)
    assert math.isclose(d["R_events_GeV4"],
                        7.547876890384561e-17,rel_tol=5e-4)
