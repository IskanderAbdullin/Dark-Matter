"""A global one-temperature/chemical-potential Bose ansatz is NOT invariant
under the physically on-shell finite-parent four-channel Dirac Yukawa operator.

The FIRST real-Bose cell below the exact light-fermion pair kinematic
support is entirely collision-frozen. Its three nondegenerate physical
Bernstein moments constrain both dT_phi and dmu_phi to ZERO under any
global BE(E;T(t),mu(t)) shape. At least one higher, on-shell active
scalar cell has nonzero quantum collision source whenever the chemical
affinity is nonzero. Hence no global BE-only T(t),mu(t) closure can
reproduce this exact nonthermal collision source at finite parent cutoff.

This is a strict, scoped *closure non-invariance theorem*, not a
physical no-go of Bose DM, nor an actual infinite thermal tail model.
"""
import math
import json
import numpy as np
from kinetics.phi_bose_five_dirac_full_yukawa_implicit import FullDiracPhysicalMoments


def model(P=2.,phi_IR_composite=True):
    return FullDiracPhysicalMoments(
        P=P,Ppair=P,N1=3,N2=3,Nphi=8,nq=18,nrad=4,nangle=8,
        pair_nrad=4,pair_nrelative=8,pair_ncm=32,pair_naz=4,
        phi_IR_composite=phi_IR_composite)


def frozen_thermal_susceptibility(m,eta5,T=1.,mu=-.2):
    """First physical scalar cell dM_a/d(mu,T) at current BE shape."""
    sp=m.species[4]
    e=np.asarray(eta5[4])
    G=sp.evaluate(e,with_G=True)["G"][0]
    # m(E;T,mu) and η(E) are affine in physical energy and exactly
    # reconstructed by Bernstein P2 at every physical p.
    eta_nodes=sp.nodal_energies.reshape((-1,3))[0]
    A=np.column_stack((G@np.ones(3)/T,
                       G@(eta_nodes-mu)/T**2))
    singular=np.linalg.svd(A,compute_uv=False)
    return dict(first_frozen_dM_dmu_dT=A,
                singular_values=singular,
                condition=float(singular[0]/singular[1]),
                rank2=bool(singular[1]>1e-10*singular[0]))


def audit(T=1.,mu_phi=-.2,P=2.):
    m=model(P,True)
    etas=m.initial_etas(T=(T,)*5,mu=(0.,0.,0.,0.,mu_phi))
    source=m.source(etas)
    Cphi=source["C"][m.offsets[4]:m.offsets[5]]
    A=frozen_thermal_susceptibility(m,etas,T,mu_phi)
    return dict(schema_version=1,P_GeV=P,T_GeV=T,
                mu_phi_GeV=mu_phi,
                scalar_first_energy_cell=m.species[4].edges[:2].tolist(),
                exact_first_cell_microcollision_moment_rate=Cphi[:3].tolist(),
                supported_phi_max_abs_collision_rate=float(max(abs(Cphi[3:]))),
                supported_phi_total_source=float(sum(Cphi[3:])),
                frozen_first_cell_BE_mu_T_jacobian=A["first_frozen_dM_dmu_dT"].tolist(),
                independent_singular_values=A["singular_values"].tolist(),
                rank_two_for_frozen_BE_mu_T=A["rank2"],
                condition_number=A["condition"],
                microscopic_H=source["H"],
                true_net_dirac_U1_defect=source["exact_event_U1"],
                exact_on_shell_energy_defect=source["exact_event_energy"],
                theorem="Nonzero active source + zero all frozen moment rates + full rank frozen (mu,T) Jacobian => no one-global-thermal-BE (mu(t),T(t)) evolution can satisfy all physical collision moments",
                limits="Hard finite per-parent P and restricted vacuum four-channel model; no full thermal parent UV, no SM/FRW/thermal 86eV DM, no zero-mode dynamics",
                class_ready=False)


if __name__=="__main__":
    print(json.dumps(audit(),indent=2,allow_nan=False))
