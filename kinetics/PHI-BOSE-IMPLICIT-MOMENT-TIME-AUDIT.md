# Entropy-stable finite-step χ1↔χ2+real φ nonlinear moment time: 86-eV IR test and physical vacuum Yukawa weak matching (2026-11-22)

## Falsifiable question, scope and mathematical method

The preceding [physical φ Bose moment representation](PHI-BOSE-MOMENT-GALERKIN-AUDIT.md) avoids a false 86-eV scalar rest-node abundance, and proves on-shell *instantaneous* entropy, two conserved dark charges, and energy. Yet this did not provide finite-time evolution: the nonlinear physical momentum-integrated moments M_i(η)=∫dΠ h_i(E)n(E;η) must be inverted back to their Fermi/Bose inverse-fugacity nodes η, without making Bose ηφ nonpositive or silently clipping f/b. The next scientific gate was to derive and actually solve a finite-step implicit microcollision for χ1↔χ2+φ, with independent source Jacobian and time-order checks, and then normalize the parent event weight against an existing independent Yukawa weak reference.

All results here are **isolated 1↔2 on-shell Yukawa only**. No χ1/χ2↔SM elastic bath in the same solver, finite ν/e±γ/Z′ plasma, FRW/H, physically fitted Yukawa coupling, zero-mode condensate, Ω or frozen f(q).

The energy hats h_i(E) obey Σh_i=1 and ΣE_i h_i(E)=E, with positive Bose ηφ(E) and unconstrained finite fermion η. The full quantum F/R factors are

    F = f1(E1)(1−f2(E2))(1+bφ(Eφ))
    R = (1−f1(E1))f2(E2)bφ(Eφ)
    J = F−R,

with exact boosted shell E1=E2+Eφ. The same positive hats construct distributions and signed collision-source moments C_i=Σevents WJ[−h1_i,+h2_i,+hφ_i]. Thus C conserves N1+N2 and N1+Nφ and the physical total energy. The instant entropy production is η·C=Σevents W(F−R)ln(F/R)≥0.

**Finite-step method:** solve the fully coupled nonlinear backward Euler moment equation

    M(η_new) − M(η_old) − Δt C(η_new) = 0.

The exact analytic physical susceptibility ∂M_i/∂η_j=−∫dΠ h_i h_j n(1±n), plus a separately derived analytic full Pauli/Bose event Jacobian, form a dense nonlinear Newton system. Each Newton line search accepts a candidate only when all Bose ηφ_i>0 and the full moment-space residual decreases; fermion f_i are automatically in (0,1). Invalid candidates fail closed. There is **no occupation clipping, FD thermal reset, post-hoc charge correction or entropy adjustment**.

For the entropy maximum-closure manifold, S(M) is concave and ∇_M S=η. Therefore, for an admissible root of backward Euler,

    S_new−S_old ≥ η_new·(M_new−M_old)
                = Δt Σevents W(F−R)ln(F/R) ≥ 0.

This is a mathematical finite-step quantum H theorem up to the reported nonlinear solve and physical quadrature errors. The implementation evaluates actual integrated entropy and the lower bound independently. It is **first-order in time**; finite-η treatment does not solve bφ=0 (η→∞) or true Bose saturation/δ³(p).

## A. Synthetic positive event-weight finite-step verification

Parameters (illustrative common energy scale): mχ1=.1, mχ2=.02, mφ=86e−9, common energy upper4, 16 hat intervals; chemically consistent μ1=0,μ2=.01,μφ=−.01 and a nonthermal χ1 inverse-fugacity distortion +.12sin(3E). Sixteen on-shell moving events carry *arbitrarily chosen* positive W=10⁻⁴, not a physical Yukawa width. Directly integrated physical Bose/Fermi N and S were checked for every accepted implicit step.

| One implicit Δt | Newton iterations | ΔNφ | Physical S increase | Concavity lower bound Δt·H |
|---:|---:|---:|---:|---:|
| .25 | 4 | −1.4798331×10⁻⁵ | +1.2816095×10⁻⁶ | +1.1503269×10⁻⁶ |
| .5 | 4 | −2.4300591×10⁻⁵ | +1.9998937×10⁻⁶ | +1.6596990×10⁻⁶ |
| 1 | 4 | −3.6707731×10⁻⁵ | +2.7933229×10⁻⁶ | +2.0318554×10⁻⁶ |
| 2 | 5 | −5.0417717×10⁻⁵ | +3.4810446×10⁻⁶ | +2.0604255×10⁻⁶ |

The physical entropy change exceeds the positive lower bound; N1+N2, N1+Nφ and total energy agree to ~10⁻¹⁶; minimum Bose ηφ stays >.0141 and nonlinear residual ≤5×10⁻¹⁶. A common-temperature nonzero chemical potential Yukawa stationary point stays stationary. Three independently selected analytic source-Jacobian columns agree with finite differences; an explicitly negative Bose inverse fugacity is **rejected** rather than silently clipped.

At fixed total synthetic τ=2, successive 2/4/8/16/32/64 fully implicit steps produce ΔNφ=−5.6540/−6.0413/−6.2621/−6.3805/−6.4419/−6.4731 ×10⁻⁵. Adjacent self-differences ≈(3.8726,2.2087,1.1839,.6134,.3123)×10⁻⁶ have last halving ratio **1.964**, consistent with first-order temporal convergence on this fixed finite event graph.

**Negative mesh result:** for the same τ=2 and 8 time steps, energy-hat intervals N=8/16/32/64 yield ΔNφ=−5.9929/−6.2621/−2.8300/−1.6445 ×10⁻⁵. The sparse 16-event quadrature is **not converged in φ production under mode refinement** even though its H/charges/energy and temporal consistency tests pass. A physical collision cannot be inferred from this graph.

## B. Actual tree-level vacuum Yukawa normalization, matched independent reference

The scalar Yukawa vertex L=−y φ(barχ2χ1+h.c.) gives

    Γ_rest = y²[(m1+m2)²−mφ²]√λ(m1²,m2²,mφ²)/(16πm1³).

For the **hypothetically chosen** y=.001, m1=.1,m2=.02,mφ=86e−9 in GeV, Γrest=2.7501974166242436×10⁻⁹ GeV. The parent integration includes physical number-density measure g1 p²dp/(2π²) with g1=4, Lorentz time dilation Γrest m1/E1 and isotropic rest-frame decay dc/2. Physical daughter energies are boosted, and full Fermi Pauli plus Bose stimulation/inverse process remain present. This is an actual conditional tree-level source but **not a cosmologically viable parameter selection**.

At synthetic Tχ1=.9,Tχ2=Tφ=1, pmax=2, nrad12,nangle12, the moment-projected weak net scalar rate is **Rφ=−1.16493012862968×10⁻¹¹**, microscopic entropy +1.64846911073637×10⁻¹² (natural units, not Γelastic). An independently coded p/c integral AND the existing project [moving Yukawa 5D weak reference](moving_yukawa_bose_pauli_weak.py) reproduce this Rφ to better than 2×10⁻¹² relative. For mφ=.001, Γrest=2.749696097519×10⁻⁹ and Rφ=−1.1603015696×10⁻¹¹. Source radial nrad8/12/16/24/32 yields Rφ=−1.16494065/−1.16493013/−1.16492995/−1.16492997/−1.16492997 ×10⁻¹¹, a scoped smooth FD benchmark. Full physical plasma and UV validity are not verified by this normalization.

One physically weighted implicit Yukawa step, at illustrative durations 10⁶/10⁷/10⁸ GeV⁻¹ and 8/16/32 scalar energy intervals, preserves physical N/energy/entropy and Bose admissibility. For 10⁷ GeV⁻¹, the scalar number change is respectively **−8.5256899×10⁻⁵, −8.0680326×10⁻⁵, −7.7733958×10⁻⁵**. Thus even with independently matched physical parent weak Γrest, N8→32 time drift differs **≈9.7%**, N16→32 **≈3.8%**. Physical rate normalization is **not** proof of convergence of the nonlinear mode dynamics. No hidden freezing or artificial condensate occurs.

## Permanent evidence, reproduction, literature and NEXT_ACTION

[Executable code](phi_bose_implicit_moment_step.py), [9 focused tests](../tests/test_phi_bose_implicit_moment_step.py), [machine-readable primary result](results/phi-bose-implicit-moment-step-2026-11-22.json), [dated stage](../archive/stages/2026/11/2026-11-22-phi-bose-implicit-moment-step.md).

    OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_phi_bose_implicit_moment_step.py tests/test_phi_bose_moment_galerkin.py
    OPENBLAS_NUM_THREADS=1 python -m kinetics.phi_bose_implicit_moment_step --audit

Nine local targeted tests (plus eight from the preceding representation) passed. New GitHub Scientific CI must be checked after publication. Methodological literature on quantum Bose/Fermi collision inversion: [stiff quantum Boltzmann](https://numdam.org/articles/10.1051/m2an/2011051/) and [conservative/entropic boson methods](https://arxiv.org/abs/1009.2748); no external implementation was copied.

**Missing physics/numerics:** The current solver evolves χ1,χ2,real φ but **no shared SM lepton bath in the same update**, no physical e±/ν weak chemistry, Z′ mediator distribution, cosmological expansion, real thermal masses or independently constrained gauge/Yukawa parameters. The zero-occupation bφ=0 and Bose-critical ηφ→0 boundaries need explicit separate mathematical and delta-mode treatment. A full angular/radial/mode/time uncertainty budget for *coupled* χ+SM+φ is still absent. Therefore **no** physical Γelastic, 86-eV ΩDM/BEC, frozen f(q) or CLASS-ready result.

**NEXT_ACTION:** replace the sparse finite event graph by a 3+ resolved moving parent p/c quadrature and independently verify strong-mode q/IR/dt convergence at fixed physical Γrest; then add χ1/χ2 elastic Pauli sharing the SAME physical moment/entropy variables, with real SM ν/e±γ/Z′ heat+FRW, testing total energy and dynamically produced χ2/φ without any thermal resets or false BEC.

## 2026-11-23 scientific correction: 12×12 Yukawa source is not angle-converged

The quoted physical tree-level Rφ=−1.16493012862968×10⁻¹¹ at pmax2 is **the correct value for the specific 12×12 Gauss event graph**, but it is not a converged continuum Yukawa rate. A new [independent cut-cell angular weak calculation](PHI-BOSE-CUTCELL-WEAK-IR-AUDIT.md) finds total Rφ=−1.29545395549×10⁻¹¹ for exactly the same thermal state, a ~10.08% magnitude underestimate in the original sparse graph. The N16 quadratic first scalar cell has a physically allowed source Rφ=−2.76214936×10⁻¹² even though the sparse 12×12 graph deposits zero there. Consequently the reported independent reference agreement at 12×12 validates implementation consistency at **fixed quadrature**, not continuum convergence, and finite-time N/time results must be rechecked with cut-cell resolved angular events. No physical relic or BEC inference is warranted.
