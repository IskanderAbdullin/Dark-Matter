# Fully coupled **continuum-time quantum mode generators** for nonthermal physical Yukawa 1↔2 (2026-11-24)

## The precise scientific gate

Prior progress included physically normalized Yukawa **thermal weak** cut-cell sources, two independent positive conditional strong subflows (scalar with frozen fermions; χ2 with frozen parent+scalar), and a globally exact two-sided inverse parent rapidity theorem. The outstanding question is whether the **same physical invariant microevents** define *simultaneous* admissible continuum-time gains/losses for all three species for **arbitrary nonthermal occupations** with physical p²dp mode capacities, two charges, total energy and quantum H theorem. This stage derives them and checks independent scalar-energy-first weak moments. **A finite-time simultaneous nonlinear solver is NOT constructed here.** No subflow is silently counted as energy conserving over time when the other species is frozen.

Model: χ1 (M=.1 GeV, g1=4) ↔ χ2 (m2=.02 GeV, g2=4) + real φ (mφ=86e−9 GeV, gφ=1), unpolarized isotropic inclusive Yukawa rest width Γrest from coupling y=.001. Internal-state factors are a conditional convention that must be checked against a microscopic spin/channel realization before physical predictions. On-shell events satisfy E1=E2+k exactly.

## Single common invariant event and physical strong mode equations

Parent p, angle c and boosted scalar energy k define exact event energies E1, E2=E1−k and k. The positive event rate per volume is

    dW = g1 Γrest M/(2π² E1) p² dp (dc/2).

Quantum gain and reverse weights are respectively

    F=f1(E1)(1−f2(E2))(1+b(k)),
    R=(1−f1(E1)) f2(E2) b(k), J=F−R.

For each actual event deposit one common signed extent dW J: dN1=−dW J; dN2=dNφ=+dW J; dQ1=−E1 dW J, dQ2=+E2 dW J, dQφ=+k dW J. Both dark number charges and energy are exact at the **continuum event level**. The microscopic entropy term dW(F−R)ln(F/R)≥0 if all occupations are interior and positive. No arbitrary daughter capacity is allowed.

The **modewise** continuum gain/loss equations arise from this same event measure after dividing by the species' physical density of states. Define per species N-mode DOS: parent g1 p²dp/(2π²); daughter g2 E2 sqrt(E2²−m2²)dE2/(2π²); scalar gφ k sqrt(k²−mφ²)dk/(2π²).

1. Parent: dot f1(p)=G1(p)(1−f1)−L1(p)f1, with G1=Γrest M/E1 ⟨f2 bφ⟩_c, L1=Γrest M/E1 ⟨(1−f2)(1+bφ)⟩_c (angle average on true rest shell). G1,L1≥0, hence both Pauli boundaries point inward.
2. Daughter χ2: dot f2(E2)=G2(E2)(1−f2)−L2(E2)f2. Both are positive physical integrals over the **exact two-sided inverse parent rapidity** interval p∈[p_low(E2),p_high(E2)] with proportional factor g1 Γrest M²/[2g2 kstar E2 sqrt(E2²−m2²)] and integrand p dp/E1 times respectively f1(1+bφ), (1−f1)bφ. At fixed E2, p dp/E1=dk for k=E1−E2; evaluate with positive log-k weights. G2,L2≥0 and the Fermi bounds are inward.
3. Scalar: dot bφ(k)=Sφ(k)−Dφ(k)bφ(k), where Sφ≥0 integrates f1(1−f2) over exact **two-sided** parent interval; Dφ integrates f2−f1 and may have either sign. The scalar prefactor is g1 Γrest M²/[2gφ kstar k sqrt(k²−mφ²)] and the integration weight p dp/E1. The physical Bose lower boundary b=0 is inward. D<0 is permitted and corresponds to genuine stimulated positive growth, not a negative occupation.

At common T with μ1=μ2+μφ, all three per-mode gain/loss ratios reproduce the correct FD1, FD2 or BEφ equilibrium independently. There is no hidden fermion or scalar multiplier introduced beyond the explicit physical DOS ratios.

## Independent physically normalized nonthermal N/Q/H

Use an intentionally nonthermal, positive occupation trio:

    f1(E)=expit(−E/.9+0.1sin(3E)),
    f2(E)=expit(−E+0.07cos(2E)),
    bφ(k)=[1+0.03sin(3k)]/expm1(k+0.005).

All energies in GeV. These nonthermal functions are NOT cosmological solutions or observations. The parent-first source integrates the **single same-event J** with 96×96 positive Gauss quadrature and yields:

| Parent P GeV | Rφ (GeV⁴) | Qφ (GeV⁵) | Q1 (GeV⁵) | Q2 (GeV⁵) | H≥0 |
|---:|---:|---:|---:|---:|---:|
| 2 | −9.667552981841565e−12 | −3.3328861369172934e−12 | +1.372375649268114e−11 | −1.0390870355763847e−11 | +1.794977527897187e−12 |
| 8 | −2.372668092370616e−11 | −1.369125131038746e−11 | +6.386183437836145e−11 | −5.0170583067973955e−11 | +7.414064029620446e−12 |
| 16 | −2.3888382051541185e−11 | −1.412416887039840e−11 | +6.534061449282427e−11 | −5.121644562242588e−11 | +7.576898137925836e−12 |

The Q1+Q2+Qφ numerical residual is 0/3.23e−26/−1.29e−26 GeV⁵ respectively; this is a *weak event-ledger* conservation check, not numerical finite-time conservation of a discretized state. Internal number charges vanish algebraically on each event.

A **separately implemented order-reversed** scalar-energy-first positive φ DOS integral, with physical two-sided parent support, was compared at k/parent Gauss levels 16,32,48,64,96. Worst relative discrepancy over Rφ,Qφ:

| P GeV | n16 | n32 | n48 | n64 | n96 |
|---:|---:|---:|---:|---:|---:|
| 2 | 1.014e−7 | 5.91e−13 | 6.22e−15 | 1.24e−14 | 1.49e−14 |
| 8 | 5.080e−4 | 9.364e−7 | 2.414e−10 | 1.37e−12 | 1.44e−14 |
| 16 | **3.665e−3** | 5.93e−5 | 5.412e−7 | 6.72e−8 | **2.728e−11** |

The coarse P16 error of 0.3665% is an explicit **negative grid-convergence control**, not a physics effect. This third-party mode order reversal uses only the same derived microscopic model and physical DOS, not a fit to the parent Gauss output. Accuracy requires convergence of the parent outer grid too: the quoted parent-first 96×96 benchmark is a numerical reference, not a proven arbitrary-nonthermal continuous exact integral.

All four local focused tests passed in 2.13s, independently covering Pauli inward directions at f=0,1, Bose inward b=0 with D<0 Bose-stimulated case, the 3-cutoff shared Q/H ledger, third-order scalar mode integral N/Q mesh convergence, and common-T FD/BE per-mode equilibrium. Public Scientific verification CI must be checked separately after commit.

## Scientific limits and NEXT_ACTION

Sources: [continuum generator](phi_bose_coupled_continuum_modes.py), [tests](../tests/test_phi_bose_coupled_continuum_modes.py).

**This is a genuine joint continuum collision vector-field derivation but NOT a solved finite-time strong physical scheme:** the parent/daughter/scalar spectra remain independent functions with no numerically integrated simultaneous evolution. The displayed equalities do not prove an arbitrary discretization will maintain physical mode capacities, Bose positivity, Pauli bounds or entropy at finite time. One-point/one-bin energy closure is known to be wrong by 45.16% in the first P2,N16 benchmark, and the affine DG-P1 Bose ansatz fails for a positive physical exact scalar subflow at P2,N32. Off-grid projection and extra disconnected graph invariants must be controlled before claiming a physical multi-species strong solver.

**NEXT_ACTION:** build and analyze a **single shared physical on-shell event list**, with all three species' physical mode capacities/energy resolution, generating the continuum formula and a discrete quantum forward/reverse collision operator. Prove finite-time strict Bose/Pauli realizability, exact two charges and energy, entropy nondecrease, correct μ stationarity and no spurious frozen components. Test weak N/Q convergence P2/8/16 and mesh/dt/IR+UV. If a particular positive finite-volume interpolation cannot satisfy this, report the obstruction rather than substituting arbitrary synthetic capacities. Only then proceed to finite SM/FRW, actual freeze-out/relic Ω, CLASS and observational comparisons.
