# Discrete quantum KL certificate for fully coupled implicit χ1↔χ2+φ time evolution (2026-11-24 UTC)

## An independently checkable physical finite-time identity

The newly committed physically normalized quantum Bernstein-P2 three-species solver solves a *common-event* backward-Euler moment root

    R = M(η_new) - M(η_old) - Δt C(η_new) = 0,

where each species' physical moments are M_a=∫ (g p²dp/2π²) B_a(E) f(E;η(E)), with its own massive on-shell density of states. The entropy-dual η(E)=ΣB_a(E)η_a is finite for fermions and **strictly positive for the Bose scalar**. Quantum forward/reverse event weights are F=f1(1−f2)(1+bφ), R_event=(1−f1)f2 bφ, shared physical on-shell event weights W>0, J=F−R_event. All statements here are for the illustrative vacuum Yukawa dark sector, no FRW/SM.

For fermions ψ_F(η)=log(1+exp(−η)); for bosons ψ_B(η)=−log(1−exp(−η)) (valid only η>0). Both are **strictly convex**, with f=−ψ' and positive ψ''=f(1∓f). Local entropy s(f)=ψ(η)+η f. For any two positive admissible η_old,η_new, the exact modewise Bregman/Kullback-Leibler gap is

    D_quantum(old||new)
      = ψ(η_new)-ψ(η_old) +(η_new−η_old) f(η_old) >= 0.

For FD this is Bernoulli KL:
    f_old ln(f_old/f_new)+(1−f_old)ln[(1−f_old)/(1−f_new)].

For BE this is geometric KL:
    b_old ln(b_old/b_new)−(1+b_old)ln[(1+b_old)/(1+b_new)].

Integrate this gap with the **same physical DOS quadrature weights** as the entropy and moments, for every energy cell of all three species. The exact finite-dimensional identity is

    S_new - S_old - Δt H_new
       = Σ_species ∫physDOS D_quantum(old||new)
         + η_new · [M_new−M_old−Δt C_new],

because η_new·C_new=Σ_events W(F−R_event)(−η1+η2+ηφ)
=ΣW(F−R_event)ln(F/R_event)=H_new≥0.

Thus for an EXACT **positive** simultaneous backward-Euler moment root the residual term vanishes and

    S_new − S_old >= Δt H_new >= 0.

This is a **conditional theorem**, not a guarantee that a positive nonlinear root exists for arbitrary Δt, nor that Newton converges to it. If the numerical root has nonzero residual, the sign of the **residual work** η_new·R matters: checking a small Euclidean norm alone is not a proof of a precise entropy inequality. The independent KL decomposition provides a sharper numerical certificate and explicitly diagnoses errors from finite nonlinear solver tolerance.

## Independent implementation and tests

The new [separate physical FD/BE KL evaluator](phi_bose_implicit_entropy_kl.py) constructs the relative entropy **independently of the collision event list**, using the physical energy-mode DOS and the canonical quantum partition functions rather than reusing the solver's entropy finite difference. It then compares the full identity against the independently assembled sparse implicit residual. It tests both **off-root perturbed spectra** (where residual work is nonzero) and converged physically admissible roots, plus exact common-T nonzero-chemical-potential equilibrium. The separate [five focused implicit regression tests](../tests/test_phi_bose_coupled_bernstein_p2_implicit.py) also independently finite-difference the full interspecies sparse Jacobian, test large time steps, P2/8/16 thermally and nonthermally, and first-order dt refinement. [Four specific KL tests](../tests/test_phi_bose_implicit_entropy_kl.py) check this stronger identity.

Reproduction:
- OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_phi_bose_coupled_bernstein_p2_implicit.py tests/test_phi_bose_implicit_entropy_kl.py
- OPENBLAS_NUM_THREADS=1 python -m kinetics.phi_bose_implicit_entropy_kl

The previous implicit solver commit 6dfdea2a87b7406d692a914d8aad4cebfca57b14 has public Scientific verification **success**. New independent tests and KL evaluator were committed later and need their own CI results to be verified; do not assume passing before inspection.

## Context and limits

General Boltzmann/quantum Boltzmann positivity, entropy, and conserved moments are nontrivial independent numerical requirements; compare [Pareschi–Rey, SIAM J. Numerical Analysis, 2022, DOI 10.1137/21M1423452](https://doi.org/10.1137/21M1423452), [Zhang–Shen–Hu, 2024 SAV entropy-preserving schemes, arXiv:2408.16105](https://arxiv.org/abs/2408.16105), and [Discrete Quantum Kinetic Equation, 2023, DOI 10.1007/s44007-023-00070-1](https://doi.org/10.1007/s44007-023-00070-1). Those works concern different operators/schemes and do not independently certify this Yukawa kernel.

The finite-event H theorem and root KL inequality are exact on the **represented physical discrete state**, but remain conditional on the chosen internal-state degeneracies, vacuum Γ, finite parent cutoff, positive Bose entropy-dual family, source quadrature and physically valid root. No SM heat, expansion, actual late-time freeze-out, 86-eV condensate, Ω, CLASS transfer function, Lyman-alpha or halo likelihood is derived.

**NEXT_ACTION:** independently verify new CI, preserve concrete root residual/KL values, and test implicit solver at physically challenging Bose near-critical states, dt and q-grid/IR/UV convergence. Do not call the scheme unconditionally positive or cosmologically validated.
