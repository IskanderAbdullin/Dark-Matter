"""Physical five-field counterexample to an assumed globally BE thermal shape.

Cannot infer chemical kinetics from a single homogeneous mu(T) when the
actual physical pair+decay operator has frozen massive-Bose IR modes.
"""
import math
import numpy as np
from kinetics.phi_bose_frozen_ir_thermal_manifold_no_go import (
    model, frozen_thermal_susceptibility, audit,
)
from kinetics.phi_bose_physical_cutcell_weak import MPHI


def test_exact_frozen_bose_cell_and_active_source_with_nonzero_affinity():
    for P in (2.,8.,16.):
        m=model(P,False)
        e=m.initial_etas(T=(1.,)*5,mu=(0.,0.,0.,0.,-.2))
        src=m.source(e)
        a=m.offsets[4]
        frozen=src["C"][a:a+3]
        active=src["C"][a+3:]
        assert np.max(abs(frozen))==0
        assert np.max(abs(active))>1e-17
        assert src["H"]>0
        assert abs(src["exact_event_energy"])<1e-17
        assert abs(src["exact_event_U1"])<1e-18


def test_two_parameter_Bose_T_mu_frozen_moment_jacobian_full_rank():
    m=model(2.,True)
    e=m.initial_etas(T=(1.,)*5,mu=(0.,0.,0.,0.,-.2))
    jac=frozen_thermal_susceptibility(m,e,1.,-.2)
    assert jac["rank2"]
    assert np.linalg.matrix_rank(jac["first_frozen_dM_dmu_dT"])==2
    # Both physical Bose chemical and thermal shape sensitivities are
    # linearly independent in *any* finite-width occupied energy cell.
    # If all 3 frozen moment derivatives vanish, dmu=0 and dT=0.
    for v in (np.array([1.,0.]),np.array([0.,1.])):
        assert np.max(abs(jac["first_frozen_dM_dmu_dT"]@v))>0


def test_IR_chemical_susceptibility_rises_in_nearcritical_limit_while_source_zero():
    m=model(2.,True)
    values=[]
    for delta in (1e-12,1e-14,1e-16):
        e=m.initial_etas(T=(1.,)*5)
        e[4]=m.species[4].nodal_energies-MPHI+delta
        source=m.source(e)
        assert max(abs(source["C"][m.offsets[4]:m.offsets[4]+3]))==0
        jac=frozen_thermal_susceptibility(m,e,1.,MPHI-delta)
        chi=float(jac["first_frozen_dM_dmu_dT"][0,0])
        assert chi>0
        values.append(chi)
        assert source["H"]>=0
    assert values[2]/values[1]>5
    assert values[1]/values[0]>4


def test_global_BE_thermal_closure_strictly_inconsistent_with_full_active_source():
    v=audit()
    assert v["rank_two_for_frozen_BE_mu_T"]
    assert v["supported_phi_max_abs_collision_rate"]>0
    assert max(abs(x) for x in v["exact_first_cell_microcollision_moment_rate"])==0
    assert abs(v["true_net_dirac_U1_defect"])<1e-18
    assert abs(v["exact_on_shell_energy_defect"])<1e-17
    print("PHYSICAL_5FIELD_GLOBAL_BE_NOT_INVARIANT",repr(dict(
        P=v["P_GeV"],mu_phi=v["mu_phi_GeV"],
        frozen_jac_singular=v["independent_singular_values"],
        active_phi_source=v["supported_phi_total_source"],
        entropy_H=v["microscopic_H"])),flush=True)
