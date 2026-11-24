# First physically capacity-normalized **simultaneous** three-species finite-time Yukawa solver (2026-11-24 UTC)

## The actual scientific advance

Previous verified stages established the *continuum* joint χ1↔χ2+φ gain/loss vector field, exact Rφ and all three weak energy rates, quantum entropy production, a globally exact massive daughter inverse parent support, and conditional positive χ2 and φ subflows with the other two species frozen. None evolved all three quantum distributions **simultaneously**, with actual physical phase-space capacities and a closed number/energy ledger, even for one finite step. This stage implements a small first **fully coupled finite-time nonlinear moment solver** in a nonoverlapping energy-cell representation. The Yukawa kernel and occupations are still illustrative: this is a physically normalized **dark-sector-only numerical prototype**, NOT realistic SM/FRW freeze-out, an actual 86-eV relic or Bose condensate prediction.

Illustrative physical parameters: χ1 mass M=.1 GeV, χ2 mass m2=.02 GeV, real φ mass mφ=86e−9 GeV, Yukawa y=.001, g1=g2=4,gφ=1 with an isotropic inclusive rest-frame decay Γrest. Parent cutoff P=2/8/16 GeV, fermion energy cell grids N1=N2=4, scalar initial quadratic Nφ=8/12/16, corrected by inserting the **physical on-shell** scalar support edges kmin(P),kmax(P), physically inactive scalar intervals retained as actual independent modes.

## Physical per-species capacity and shared exact on-shell event geometry

Within every nonoverlapping massive energy cell [L,U] of species s, use the nonnegative quadratic Bernstein energy-hat basis

    B0=(1−x)^2, B1=2x(1−x), B2=x^2, x=(E−L)/(U−L).

Exactly ΣBa=1 and Σea Ba=E for the *physical* nodal energies (e0,e1,e2)=(L,(L+U)/2,U), for every E in the cell. Physical finite-volume moments are NOT arbitrary single-event counts:

    Msa(eta)=gs/(2π²)∫cell p² dp Ba(E(p)) fs[eta(E(p))],

with massive on-shell E=sqrt(p²+ms²). FD occupations are f=1/(expη+1), Bose b=1/(expη−1). A nonnegative Bernstein basis and strictly positive scalar η coefficients guarantee Bose η(E)>0 and b≥0; fermions are automatically Pauli-bounded. The physical DOS per energy is equivalent to gs E sqrt(E²−m²)/(2π²). Physical susceptibility per cell G_ab=gs/(2π²)∫p²dp BaBb f(1±f) is positive definite on nondegenerate cells. Physical entropy is integrated from exact quantum FD/BE entropy densities with the SAME DOS. In a fixed cell, entropy variables η are dual to these physical moments.

For each on-shell parent momentum p and boosted scalar rest angle/cosine, calculate actual daughter k, E2=E1−k. Split the radial parent integral **at exact scalar energy-bin shell crossings**, integrate each positive partially intersected scalar-energy cell with log-energy Gauss nodes (no sparse-angle false zero). Assign one and the same positive event rate W_a=g1 Γrest M/E1 p²dp/(2π²) *dk/(2B). At each microevent evaluate the three species' OWN physical Bernstein hats at true E1,E2,k, and quantum event flux

    F=f1(E1)(1−f2(E2))(1+b(k)),
    R=(1−f1(E1))f2(E2)b(k), J=F−R.

Deposit a common signed event extent WJ into the three moment vectors C1=−A1ᵀ(WJ), C2=+A2ᵀ(WJ), Cφ=+Aφᵀ(WJ). Since every species reproduces constants and physical energy exactly, the two dark charges N1+N2 and N1+Nφ, and total dark energy Q1+Q2+Qφ are EXACT at the algebraic finite-event/source level. Microscopic quantum entropy rate H=ΣWJ(−η1+η2+ηφ)=ΣW(F−R)ln(F/R)≥0 at any admissible interior state. These identities are not supplied by a post-step projection: all species share **one** physical event list and on-shell energies. Conditional identical occupation degeneracy g2=4 is an input convention, not a validated full spin/SM ultraviolet model.

## Finite-time construction (conditional, not an unconditional entropy theorem)

At state eta_old, compute the *joint* moment source C_old from the shared finite event list. Form explicit moment Euler targets M_target=M(eta_old)+dt C_old in **all three species simultaneously**. For each nonoverlapping physical cell solve its 3×3 nonlinear moment-inversion equations

    M(eta_new)−M(eta_old)=dt C_old

using the exact physical Bose/Pauli integral moment map, physical positive-definite susceptibility and a damped Newton method with strict ηφ coefficients >0. The direct difference M(eta_new)−M(eta_old) is integrated pointwise to avoid cancellation, not inferred by subtracting large thermal moments. An explicit physical quantum entropy integral is recomputed after the step. The adaptive driver **rejects** any trial with failed local positive inversion, significant physical charge/energy error, or numerically decreasing quantum entropy; it halves dt and retries from the same old state. The small cell-source roundoff rule (skip target |dt Ccell| less than max[1e−23,32 eps×|Moldcell|]) avoids spurious failures at exact chemical equilibrium; it introduces a documented floating-point resolution threshold, so no claim of *machine-exact* conservation for arbitrarily small physical moments is made.

**Important:** this is a *moment-explicit* finite step. The continuum H theorem does NOT guarantee global entropy increase under arbitrary explicit dt. The finite-step entropy and admissibility checks are A POSTERIORI, not a proof of unconditional stability, convergence of all modes or existence for every state. The method is a constructive positive physical-capacity **prototype**, not a fully independently converged cosmological strong simulation.

## Independent benchmarks and negative numerical evidence

**Weak moment normalization against independent analytical thermal reference** is verified at all parent cutoffs P=2/8/16 and three actual resolution triples (parent Gauss per smooth segment, scalar cell log-angle Gauss, original scalar N): (4,8,8),(6,12,12),(8,16,16). The worst relative error among Rφ,Qφ,Q1,Q2 for each cutoff is:

| P GeV | 4×8×8 | 6×12×12 | 8×16×16 |
|---:|---:|---:|---:|
| 2 | 5.05e−10 | 3.32e−13 | 3.81e−13 |
| 8 | 6.44e−8 | 1.54e−10 | 7.82e−14 |
| 16 | 1.51e−7 | 3.21e−10 | 1.22e−12 |

These results independently validate the shared discrete physical weak source against earlier analytic log/dilog daughter reference and its **three physical energy ledgers**. Scalar Nφ16 P2 has ~15,360 common positive events, 78 physical Bernstein moments. A comprehensive event graph rank audit on nine P×grid cases yields rank(nullspace)=**three expected dark invariants plus exactly the physically disconnected scalar hat DOFs**, and NO additional numerical spurious invariants under the stated singular-value cutoff. In particular P2,Nφ16 has nullity=12 = three physical charges/energy + **nine frozen φ DOFs** outside the hard on-shell support; other tested cases typically have nullity=9 = three invariants+six disconnected DOFs. The earlier informal summary that all grids have only six disconnected modes was too broad: both 6 and 9 occur, depending on the physical support edges/grid.

**Finite-time joint physical step:** for P=2,8,16, Nφ16, physical φ entropy-variable positive throughout, after a joint 1e5 GeV⁻¹ explicit-moment step:

| P GeV | ΔS_dark (natural number density units) | Charge1 and charge2 defect, worst absolute | Total energy defect |
|---:|---:|---:|---:|
| 2 | +1.86712643124e−7 | 1.32e−15 | 8.33e−17 |
| 8 | +9.20863988174e−7 | 2.36e−16 | −3.89e−16 |
| 16 | +9.41746922223e−7 | 9.24e−15 | −3.89e−16 |

All six thermal reverse and explicitly perturbed **nonthermal** initial states at these cutoffs pass accepted finite-step physical entropy nondecrease, both dark charges, total physical energy, Bose/Pauli bounds, and positive microscopic H_old/H_new at dt=1e5 GeV⁻¹. Note the numerical *absolute* conservation errors depend on the much larger integrated thermal numbers/energies: quote the measured values, not an unqualified exact machine-zero statement.

**Equilibrium stationarity:** at common T=1 GeV, μ1=−.15, μ2=.05, μφ=−.2, all event F/R match within floating precision and fixed state remains unchanged after a nominal dt=1e5. An initial attempt to invert 1e−22 unresolvable moment differences produced spurious Newton failure; the per-cell roundoff threshold above cured this, but it is a **numerical finite-precision policy** and NOT a proof of arbitrary-precision thermodynamic equilibrium.

**Negative large-step case and recoverability:** P2,Nφ16 dt=3e5 GeV⁻¹ triggers a failed *local positive Bose Newton line search*. No mathematical impossibility of a positive solution is established from this solver failure. The adaptive driver rejects this trial and safely integrates the same total t=3e5 in two accepted dt=1.5e5 steps with positive ηφ and entropy increases.

**First-order time convergence:** at P2,Nφ16, same final t=3e5, dtmax 18750/9375/4687.5 GeV⁻¹ yields actual ΔNφ ≈−3.63805967484e−6/−3.63315687463e−6/−3.63073792422e−6, with successive difference ratio **2.02683**, consistent with first-order explicit time accuracy. This numerical convergence is only in dt at one illustrative grid/cutoff, not proof of q-grid convergence, physical late freeze-out or nonthermal continuum limit.

Five focused local pytest tests passed in 10.56s. They cover: three parent cutoffs×three physical grids, *exact physical cell DOS number capacities* versus integrated Gauss weights to 2e−13, full event graph rank/three invariants, six nonthermal+thermal finite time steps with entropy/charge/energy/positivity, common-T nonzero μ equilibrium, adaptive large-step recovery, first-order dt convergence.

## Reproducibility and strict remaining gates

- Source: [first coupled Bernstein-P2 physical solver](phi_bose_coupled_bernstein_p2_finite.py)
- Tests: [focused 5-test suite](../tests/test_phi_bose_coupled_bernstein_p2_finite.py)
- Run: OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_phi_bose_coupled_bernstein_p2_finite.py
- Run: OPENBLAS_NUM_THREADS=1 python -m kinetics.phi_bose_coupled_bernstein_p2_finite
- New public Scientific verification CI must be checked independently (previous branch head completed/success, newest source/tests were in progress at this report).

**Scientific status:** A first **genuine fully coupled three-population finite-time physical quantum collision prototype** is implemented with exact same-event conservation and physically normalized capacities, positive Bernstein-P2 occupations, accepted entropy-increasing finite steps and first-order dt convergence. It remains finite-cutoff, illustrative fixed Yukawa vacuum width with conditional g1/g2 channel counting; excludes SM and FRW, thermal masses, realistic 2→2, physical zero-mode channels, nonthermal late-time relic normalization, observational constraints and independent cross-code full time simulation. It is NOT *CLASS-ready*, proof of BEC, proof that mφ=86eV yields DM, or completion of general finite-time numerical stability.

**NEXT_ACTION:** formulate and test a truly implicit all-species **backward Euler entropy-dual moment scheme** using the same physical event measure/capacities, with exact moment residual and convergence; by concavity of the Bose/FD entropy on physical moment manifold, a valid implicit root gives ΔS≥dt H_new≥0. Verify nonlinear existence/positivity and sparse Jacobian, timestep convergence, mesh/UV/IR and graph rank. Compare to explicit moment prototype; do not promise unconditional existence without a proof.
