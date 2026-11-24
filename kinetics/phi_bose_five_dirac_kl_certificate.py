"""Independent finite-time physical FD/BE KL entropy check for the
combined FIVE-population Dirac/antiχ plus real φ Yukawa four-channel solver.

The relative quantum entropy is evaluated independently on actual
g*p²dp/(2π²) species quadrature and is NONNEGATIVE mode by mode.
For ANY trial new state, even off the nonlinear root:
  ΔS − dt*H(new) = Σ D_quantum(old||new) + η_new · residual.
At an admissible EXACT positive backward-Euler moment root R=0,
  ΔS >= dt*H >= 0.
This certifies physical entropy only for the DISCRETIZED vacuum dark
collision model and the chosen quadrature, not an unconditional root
existence theorem, nor cosmological freeze-out/BEC.
"""
import json
import numpy as np
from kinetics.phi_bose_implicit_entropy_kl import physical_relative_entropy
from kinetics.phi_bose_five_dirac_full_yukawa_implicit import (
    FullDiracPhysicalMoments,asflat
)

def quantum_KL_certificate(model,old_etas,new_etas,dt):
    if dt<=0:raise ValueError("Positive dt required")
    before=model.integrate_spectrum(old_etas)
    after=model.integrate_spectrum(new_etas)
    residual,source=model.implicit_residual(
        asflat(new_etas),asflat(old_etas),dt,with_jac=False)
    K=physical_relative_entropy(model,old_etas,new_etas)
    lhs=float(after["S"]-before["S"]-dt*source["H"])
    numerical_work=float(np.dot(asflat(new_etas),residual))
    rhs=sum(K)+numerical_work
    return dict(entropy_gain=float(after["S"]-before["S"]),
                microscopic_H_new=source["H"],
                physical_quantum_KL_by_species=K,
                root_residual_work=numerical_work,
                exact_entropy_decomposition_defect=lhs-rhs,
                entropy_certificate=lhs,
                reconstructed_certificate=rhs,
                U1_defect=float(after["U1"]-before["U1"]),
                energy_defect=float(sum(after["Q"]-before["Q"])),
                class_ready=False)

def audit():
    m=FullDiracPhysicalMoments()
    initial=m.initial_etas()
    out=m.adaptive_step(initial,10000.)
    return dict(schema_version=1,model="5 physical species 4 Yukawa channels",
                conditional_root=quantum_KL_certificate(
                    m,initial,out["etas"],out["accepted_dt"]),
                iterations=out["iterations"],accepted_dt=out["accepted_dt"],
                class_ready=False)

if __name__=="__main__":
    print(json.dumps(audit(),indent=2))
