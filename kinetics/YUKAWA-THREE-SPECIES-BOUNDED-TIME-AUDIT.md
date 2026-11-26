# Moving χ1↔χ2+real φ: a nonlinearly evolving three-species Pauli/Bose energy-matched source (2026-11-22)

## Scientific question and prior evidence

The earlier [Lorentz-boosted Yukawa weak rate](MOVING-YUKAWA-BOSE-PAULI-WEAK-AUDIT.md) independently established the vacuum scalar width, kinematics, full inverse decay, exact chemical detailed balance and energy/charge/entropy ledgers—but evaluated arbitrarily *prescribed* distributions, not a dynamical χ1/χ2/φ system. With χ2 and φ energies varying under moving decays, an arbitrary assignment of daughter production to finite mode bins can break the first energy moment or the quantum H theorem. The aim here is a falsifiable **three-species nonlinear quantum finite-volume pilot** with the original microscopic scalar decay and no phenomenological BGK term, not a claim of physical dark matter production.

Predeclared tests: (1) preserve Nχ1+Nχ2 and Nχ1+Nφ, (2) preserve the sum of three physical energy moments, (3) produce nonnegative microscopic and actual nodal entropy, (4) maintain 0<fχi<1 and bφ>0 through accepted finite steps, (5) reproduce the earlier independent **moving weak** rate on arbitrary 16/32/64 bin thermal shapes, (6) pass nonzero independent chemical μ equilibrium, and (7) demonstrate both stimulated production and net inverse absorption while retaining **numerical nonconvergence** with energy bins. These tests are strictly scoped to a synthetic isolated Yukawa reaction, not a full SM+FRW simulation.

## Derivation: one common on-shell event, three densities and discrete H

For the renormalizable real-scalar Yukawa −yφ(barχ2χ1+barχ1χ2) with m1>m2+mφ, physical Γrest is the tree-level width already independently checked in the prior audit. For parent momentum p and rest decay angle c, E1=√(p²+m1²) and Lorentz-boosted E2+Eφ=E1. Each event's **common positive phase-space coefficient** is

    dW = [g1 p² dp/(2π²)] [Γrest m1/E1] dc/2.

Its Pauli/Bose *forward* and *inverse* statistical factors are

    F=f1(E1)[1−f2(E2)][1+bφ(Eφ)],
    R=[1−f1(E1)]f2(E2)bφ(Eφ),
    J=F−R.

For each of three independent uniform physical energy grids E_{s,i}, use nonnegative barycentric α_{s,i}(E) with sum α=1 and sum α E_i=E for ANY on-shell incoming/daughter energy inside the grid. Deposit each event exactly once: dotN1_i = −dW J α1_i(E1), dotN2_j = +dW J α2_j(E2), dotNφ_k = +dW J αφ_k(Eφ). Thus the **same** on-shell physical event implies dotN1+dotN2=dotN1+dotNφ=0 and dotE1+dotE2+dotEφ=0 up to roundoff; a mere weak-moment-only check would not guarantee this.

To ensure that the quantum source and the projected *physical* entropy are consistent, interpolate Fermi LOG-ODDS and Bose inverse LOG-FUGACITY, not raw arbitrary nonlinear distributions:

    h1(E)=Σ_i α1_i(E) log[f1_i/(1−f1_i)]
    h2(E)=Σ_j α2_j(E) log[f2_j/(1−f2_j)]
    ηφ(E)=Σ_k αφ_k(E) log[(1+bφ_k)/bφ_k]
    f1(E)=expit(h1(E)), f2(E)=expit(h2(E)),
    bφ(E)=1/expm1(ηφ(E)).

Then nodewise quantum entropy S=−Σ_{a=1,2} Σ_i w_{a,i}[f ln f+(1−f)ln(1−f)] + Σ_k wφ,k[(1+b)ln(1+b)−b ln b] has

    dS/dt = −Σ_i dotN1_i h1_i − Σ_j dotN2_j h2_j + Σ_k dotNφ_k ηφ_k
          = ∫ dW (F−R) ln(F/R) >=0.

Here each physical **zeroth** radial density-of-states cell weight is

    w_{s,i}=g_s/(6π²){[U_i²−m_s²]^(3/2)−[L_i²−m_s²]^(3/2)}.

The semidiscrete nonlinear ODE is df1_i/dt=dotN1_i/w1_i, df2_j/dt=dotN2_j/w2_j, dbφ_k/dt=dotNφ_k/wφ_k. Chemical equilibrium for the CLOSED 1↔2+φ network requires μ1=μ2+μφ (μ need not vanish until other processes are included). For exact FD/BE at common T, both h and η are affine in energy; their barycentric reconstruction preserves microscopic F=R for any bin spacing.

Finite step: recompute all three species source after each provisional Euler half step. Accept two such half steps only when both intermediates and final node states remain strictly Pauli/Bose admissible, *actual* entropy does not decrease and the infinity norm of differences between one full versus two half steps for f1, f2 and log1p(bφ) is below atol. Rejected candidate halves dt, no clipping or thermal reset. **The accepted method is FIRST ORDER; local step-doubling is not a rigorous global error bound.**

## Real numerical results: exact moments but significant finite-bin uncertainty

All numbers are synthetic arbitrary common energy units, not an 86-eV physical mass: m1=.1,m2=.02,mφ=.001,y=.001, parent pmax12, Gauss radial 12, angular12, g1=g2=4 numerical Dirac+antiparticle convention, gφ=1 real boson, end time 5×10⁷ in inverse *toy* units, initial dt2×10⁷, atol=10⁻⁶. A proper physical charge/spin and mediator/thermal mass normalization is still needed.

**Weak/strong independent comparison:** with FD(T1=1.1), FD(T2=1), BE(Tφ=1), the three-species deposited full R,Q1,Q2,Qφ,H coincide with [previous independent moving weak integrator](moving_yukawa_bose_pauli_weak.py) to relative 10⁻¹⁰ for **16,32,64 energy bins** when using identical 12×16 p/CM integration; per-event on-shell E defect ~2.2×10⁻¹⁵. With T=.8, μ1=.03, μ2=.05, μφ=−.02, the separate arbitrary-number equilibrium is stationary to floating accuracy. At bins32, the forward initial dotNφ=+2.73784212332×10⁻¹¹; with colder T1=.9 the initial dotNφ=−2.65792481244×10⁻¹¹, despite quantum H being positive in BOTH cases.

Actual three-distribution nonlinear **time** evolution with Tχ2=Tφ=1 initially:

| bins | mother Tχ1 | accepted/rejected | ΔNφ=ΔNχ2 | ΔS (actual) |
|---:|---:|---:|---:|---:|
| 16 | 0.9 | 79 / 31 | −0.001219814 | +0.000346405 |
| 16 | 1.1 | 75 / 30 | +0.001275659 | +0.000328315 |
| 32 | 0.9 | 133 / 50 | −0.001194165 | +0.000337165 |
| 32 | 1.1 | 127 / 48 | +0.001268518 | +0.000324541 |
| 64 | 0.9 | 206 / 76 | −0.001111587 | +0.000310991 |
| 64 | 1.1 | 174 / 64 | +0.001241995 | +0.000314231 |

For all six tested histories, weighted N1+N2, N1+Nφ and total energy differ from their initial values by less than 2×10⁻¹⁵; accepted steps have positive entropy gain, Fermi occupations stay strictly between 0 and 1 and bosonic occupations remain >0. For bins32,T1=1.1 the entropy increase is +0.0003245412; mother N declines by 0.00126852 while BOTH χ2 and φ numbers increase by this value. **This is actual evolving φ and χ2 in a closed toy collision system**, not fixed prescribed daughters.

**Negative physical/numerical evidence:** at fixed toy time atol, 16→64 bins changes *cold-mother* ΔNφ magnitude from 0.00121981 to 0.00111159 (**~8.9%**), and hot-mother from 0.00127566 to 0.00124200 (~2.6%). The current 12-point parent radial quadrature is NOT independently converged for arbitrary nonthermal modes; the rapidly occupied very-low-energy Bose bins can require much finer IR resolution. Thus a perfect N/E/H ledger gives no trustworthy physical absolute relic abundance.

**Temporal tolerance** (bins32,T1=1.1): atol 10⁻⁵,10⁻⁶,10⁻⁷ gives accepted/rejected 40/15,127/48,396/140 and net ΔNφ=+0.00126915347,+0.00126851763,+0.00126831775. The changing result demonstrates actual nonzero first-order integration bias; local tolerance alone cannot be called global convergence.

## Scope, independent literature and next experiment

Scientific context: Wen Yin, [arXiv:2301.08735](https://arxiv.org/abs/2301.08735) (Bose enhanced eV-scale thermal DM), and Sakurai & Yin, [arXiv:2410.18968](https://arxiv.org/abs/2410.18968) (stimulated emission and freeze-in) describe why one must preserve quantum stimulation and inverse processes. The numerical implementation and tests here are original project calculations, not copied article code.

[Source](yukawa_bose_pauli_bounded_time.py), [eight focused parametrized tests](../tests/test_yukawa_bose_pauli_bounded_time.py), [six history JSON](results/yukawa-three-species-pauli-bose-time-2026-11-22.json), [three tolerance JSON](results/yukawa-bose-pauli-temporal-tolerance-2026-11-22.json), [dated stage](../archive/stages/2026/11/2026-11-22-yukawa-three-species-bounded-time.md).

    OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_yukawa_bose_pauli_bounded_time.py tests/test_moving_yukawa_bose_pauli_weak.py
    OPENBLAS_NUM_THREADS=1 python -m kinetics.yukawa_bose_pauli_bounded_time

**Important omissions:** Bose vacuum bφ=0 endpoint is *fail closed* (entropy derivative singular), not safely evolved. No zero-mode/IR limit, no FRW/Hubble/redshift or finite SM e±γν/Z′ bath and physical lepton/χ chemical heat exchange, no thermal masses or realistic couplings, and no observational BBN/ΔNeff/CMB/structure tests. Previous one-shared-lepton χ1χ2 Pauli elastic pilot is *not yet coupled in the same ODE*. mφ=.001 arbitrary numerical unit is NOT the historic 86 eV particle. CLASS-ready=false; NO physical ΩDM, BEC, freezeout, Γelastic or cosmological dark matter result. Priority: integrate the two collision generators for all four χ1,χ2,φ,ℓ dynamical distributions and independently verify combined charges/energy/H; quantify radial/IR/bin/time/cutoff before any physical inference.
