# Local detailed balance is the necessary bridge from 5D χe scattering to finite bins

2026-11-21. **Analytic consistency result, not a calibrated 5D collision kernel.** Baseline branch HEAD before this note: `8e0d5e63fda79331ea1270e88f10b2e261078b7d`. The actual GitHub Scientific verification runs [44938150282](https://github.com/IskanderAbdullin/Dark-Matter/actions/runs/44938150282) and [44938076520](https://github.com/IskanderAbdullin/Dark-Matter/actions/runs/44938076520) succeeded on the preceding commits; they **do not** test the new derivation.

## Hypothesis and derivation

For a number-conserving fermionic χ scattering from a thermal, chemically equilibrated fermionic electron bath at temperature T, discretize χ states into energies E_i and capacities g_i (number per physical volume, including phase-space measure). Let f_i=N_i/g_i, 0≤f_i≤1. Define K_ij(T)≥0 as the **per-state** forward coefficient, inclusive of the integrated electron incoming occupancy, outgoing Pauli blocking, Lorentz-invariant amplitude, flux and outgoing phase space, but *excluding* χ occupancy/blocking. If χ i→j transfers Δ=E_i−E_j to the bath, bath FD factors obey

f_e(k)[1−f_e(k′)] / {f_e(k′)[1−f_e(k)]} = exp[(E_i−E_j)/T]

when electron chemical potentials match on both sides and E_i+ε_k=E_j+ε_k′. Time-reversal-invariant |M|² and a **paired** identical quadrature imply, for a consistent per-state normalization,

K_ij / K_ji = exp[(E_i−E_j)/T].

The finite-bin net flow (number per volume per time) is

J_ij = g_i g_j [K_ij f_i(1−f_j) − K_ji f_j(1−f_i)],  J_ji=−J_ij.

A convention absorbing one g factor into K is equally valid, but then the ratio must carry the corresponding density-of-states factor. The common-T Fermi distribution f_i=[1+exp((E_i−μχ)/T)]⁻¹ makes each J_ij zero, for arbitrary conserved μχ. This is stronger than checking a summed heat source alone.

With dN_i/dt=−Σ_j J_ij and dE_bath/dt=Σ_{i<j}(E_i−E_j)J_ij, Nχ and Eχ+E_bath are conserved exactly at the semidiscrete ODE level. For the bath, dS_bath=dE_bath/T; the fermion entropy derivative is dSχ/dN_i=ln[(1−f_i)/f_i]. Thus the pair entropy production is

σ_ij = J_ij ln{K_ij f_i(1−f_j) / [K_ji f_j(1−f_i)]} ≥ 0,

because (x−y)ln(x/y)≥0 for positive x,y. The identity includes the bath energy term **only** when K_ij/K_ji=exp(Δ/T). It remains valid for a changing finite bath T(t) if coefficients are recomputed at its instantaneous T and energy is transferred, not reset. At f=0 or 1 use continuous endpoint limits.

## Crucial obstruction for real 5D-to-bin projection

The identity above is **not automatic** after mapping continuously scattered χ energies E and E′ to finite cell centers E_i,E_j. Exact microscopic electron recoil is ε′−ε=E−E′, whereas the finite bath ledger assigns E_i−E_j. Therefore naive independent interpolation of outgoing states or averaging K_ij and K_ji can break local detailed balance and produce spurious equilibrium heating, even if a weak-form energy deposition conserves an independently calculated continuous energy moment. A fixed finite cell with nonzero energy width cannot generally represent both exact microscopic energy and arbitrary nonthermal within-cell χ occupations. A thermodynamically consistent coarse model needs either (a) a declared within-cell equilibrium closure and paired projected transition kernel satisfying the bin-center ratio, with quantified projection bias, or (b) a higher-dimensional energy-resolved conservative discretization. Merely symmetrizing rates after the fact changes the underlying collision physics and requires an explicit consistency/convergence test.

## Independent minimal algebraic test

A 10,000-case seeded random numerical check of the dimensionless pair algebra (E_i,E_j∈[0.1,4], T∈[0.2,3], f_i,f_j∈[10⁻⁵,1−10⁻⁵]) using K_ij=exp[(E_i−E_j)/(2T)], K_ji=1/K_ij yielded minimum σ=0 (to floating precision); this is an **algebraic synthetic control**, not a numerical evaluation of the physical χe amplitude. Reproduce with `python kinetics/test_bin_detailed_balance_identity.py`. Additional decisive test: compute K_ij,K_ji by two *independent* physical 5D quadratures, check the ratio and normalization on adapted IR grids, and compare bin heat to the strong and weak continuous 5D operators at increasing resolution.

## Scope and negative result

No claim of actual K_ij from χe scattering, physical γ_el, convergence of the slow spectrum, BEC, ΩDM, frozen f(q) or CLASS readiness. This establishes a necessary mathematical acceptance gate for the next implementation, not a solution of the microphysical bottleneck.
