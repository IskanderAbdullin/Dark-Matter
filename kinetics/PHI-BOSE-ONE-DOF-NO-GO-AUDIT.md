# One-DOF physical Bose cell closure: monotone positive number, but exact energy and thermal equilibrium are lost (2026-11-23)

## Scientific gate: why the tempting scalar constant-η repair is insufficient

The [prior physical DG-P1 counterexample](PHI-BOSE-DG-TANGENT-COUNTEREXAMPLE-AUDIT.md) exhibited a negative first inverse-fugacity tangent **even though the local Bose moment sources were negative**. Test a genuinely single-parameter, nonoverlapping scalar energy cell rather than declaring this particular mathematical fix sufficient. One scalar η per cell is indeed monotone with respect to its own particle number, but preserving energy and quantum detailed balance is a separate obligation.

For a physical massive scalar, define within one nonzero-width cell [E₋,E₊]:

- W = g/(2π²) ∫ p² dp, Q = g/(2π²) ∫ p²dp E(p), Ē=Q/W ∈ (E₋,E₊), where p²=E²−m².
- b(η)=1/(exp η−1), M(η)=W b, U(η)=Q b=Ē M.
- J=dM/dη=−W b(1+b)<0, so a cell with a strictly negative number source has η̇>0 **at the instantaneous tangent level**. For prescribed finite sink C, the physical solution still requires M+dt C≥0; this derivative fact is not unconditional positive time solvability.
- A true on-shell microevent carries energy Eφ∈(E₋,E₊). Depositing its scalar particle transfer ΔN into this one variable stores Ē ΔN, whereas the physical shell transfers Eφ ΔN. Any event with Eφ≠Ē produces an exact energy defect ΔU=(Ē−Eφ)ΔN; cancellations are accidental. Thus a one-DOF **cellwise constant Bose η** cannot have exact particle number and exact physical energy **for all** continuous shell events in a finite-width cell.
- At genuine common T and chemical equilibrium, quantum detailed balance requires ηφ(Eφ)=(Eφ−μφ)/T at every shell point (since E1=E2+Eφ and μ1=μ2+μφ). A constant ηcell cannot match all Eφ within a nonzero-width cell, even if fitted to reproduce the exact integrated thermal Bose particle number. It therefore generates a **false common-temperature source and entropy production**. Its formal discrete H≥0 is not evidence that the physical equilibrium is stationary.

The two independent mathematical failures (energy mismatch and equilibrium mismatch) mean the one-parameter proposal is **rejected as an exact closure**, not that all entropy-constrained multi-moment discretizations are impossible. An affine η(E)=a+bE in at least two coordinates can represent the thermal equilibrium and physical energy, but the previous DG-P1 results prove that its unconstrained inverse-moment tangent can leave η>0 near the infrared boundary. A future scheme must handle the **actual realizability cone**, not simply floor η.

## Reproducible physical-Yukawa negative control

The new [self-contained executable](phi_bose_one_dof_no_go.py) independently integrates the tree χ1↔χ2+φ event weights with positive boosted log-Eφ angle quadrature. Parameters: m1=.1 GeV, m2=.02 GeV, mφ=86e−9 GeV, y=.001, g1=4, T1=T2=Tφ=1 GeV, μ=0, parent pmax=2 GeV (also 8,16), nr24 and nangle48, scalar physical energy cell edges quadratic and exact p²dp scalar number fit by 64-point Gaussian radial integration. Fermions use exact FD at every event; the control uses exact BE(Eφ) instead of the cell projection.

| Cells | Thermal scalar energy bias | False φ source | False entropy rate | Energy defect / sum of |event energy flux| |
|---:|---:|---:|---:|---:|
| 8 | +1.5585% | +9.5375×10⁻¹³ | +2.2763×10⁻¹³ | +18.08% |
| 16 | +0.3971% | +3.2252×10⁻¹³ | +5.7558×10⁻¹⁴ | +8.72% |
| 32 | +0.09978% | +5.5582×10⁻¹⁴ | +1.4809×10⁻¹⁴ | +4.31% |
| 64 | +0.02498% | +5.4293×10⁻¹⁵ | +3.7285×10⁻¹⁵ | +2.14% |
| 128 | +0.00625% | +4.0361×10⁻¹⁵ | +8.6469×10⁻¹⁶ | +1.05% |

All equal-T *continuous* weak rates in the same test are ≈4.8×10⁻²⁷ (roundoff); the finite projected cell source is purely artificial. At N32 with independently sampled nr48 and parent pmax=2/8/16, the false H is respectively 1.45×10⁻¹⁴ / 1.45×10⁻¹³ / 2.88×10⁻¹³; this is not cured by changing the parent cutoff. The absolute numerical values are conditional on this artificial test and rate normalization, not a physical prediction. The 3 focused regression tests passed **locally in an independent clean Python sandbox** in 0.10 s; repository CI after commit must be checked separately.

Machine-readable compact data: [JSON](results/phi-bose-one-dof-no-go-2026-11-23.json). Reproduction: `OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_phi_bose_one_dof_no_go.py` and `python -m kinetics.phi_bose_one_dof_no_go`.

## Decision and subsequent constructive gate

Reject one-η **constant-in-energy** cell closure for exact continuous Yukawa kinetic physics, even though its scalar sink derivative is inward. Next design a genuinely realizable energy-resolving collision representation (e.g. positive on-shell reaction extents at *actual* event energies, with shared shell nodes or an explicit convex constrained two-moment inversion), independently check exact energy, both stoichiometric charges, Bose/Pauli positivity, detailed balance and entropy, and audit **collision graph spurious invariants** before time integration. Do not smuggle energy into neighboring Fermi variables or reset the thermal distribution; that would hide the same defect.

**Scope:** no physical eV abundance, true condensate, full FRW/SM/Z′ closure, frozen f(q), CLASS or observational likelihood. The strong negative already found for N32 DG-P1 with mapped low-Eφ daughters remains intact.