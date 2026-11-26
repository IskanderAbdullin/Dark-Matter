# Five-field physical Yukawa implicit Newton: float64 precision floor and cancellation-free FD/BE affinity (2026-11-24 UTC)

## Independently diagnosed negative test and resolution

Full 4-reaction 5-population χ1, antiχ1, χ2, antiχ2, real φ implicit physical-capacity strong code failed both genuine finite-time and independent KL-root tests (Scientific run #45405888389, **1221 passed, 2 failed**). Even 13 dt halvings did not satisfy an impossible componentwise relative residual threshold. A new diagnostic from run #45421473785 identified the specific **physical scalar high-energy moment DOF 74**, not a physical Bose positivity failure:

| Physical quantity, dt=10000 GeV^-1 | Observed |
|---|---:|
| Physical moment M_old | 0.010562502499539666 GeV³ |
| Entropy dual η | 4.086731509081753 |
| Original source C_old | 5.3317122527325266e−20 GeV⁴ |
| Absolute implicit residual | −9.497121774688923e−19 GeV³ |
| Original source-scaled relative error | 0.00178125175 |
| Newton η correction demanded | +1.1567543973059635e−16 |
| Smallest representable float64 η increment | **8.881784197001252e−16** |

The demanded update was **smaller than one floating-point ULP**. The physical moment residual is ~9e−17 relative to the moment itself. Further Newton iterations or smaller time steps cannot satisfy an unresolvable flux-relative target; that is a finite-precision numerical obstruction, NOT evidence of absent positive quantum kinetic roots, nonphysical negative Bose occupations or BEC.

Retain the original per-moment desired source scale s=max(dt|C_old|,1e−14|M_old|,1e−20), but define documented *absolute moment machine tolerance* ε=16 eps_float64 |M_old|. Converge only when max_active |R|/max(s,ε/rtol) < rtol, with rtol=2e−7. Report **both** the precision-aware scaled residual and the strict raw_scaled_residual, physical max_abs_moment_residual and precision_limited_dofs. This is explicit controlled loss of unattainable local relative accuracy, not a claim that the true source has 2e−7 relative accuracy in the tails. Actual finite-time **net Dirac U1 and total on-shell energy** acceptance, physical entropy/KL residual work and Bose eta>0 tests remain independent and enforced.

Focused, actually executed GitHub CI confirms:
- #45421790808: repaired true finite-time five-field and independent quantum KL tests success.
- #45422074906: numerical ULP regression, near-equilibrium affinity, five-field conserved charges and KL tests success.
- #45422175702: P2, P8 and P16 hard parent cutoffs with fully coupled physical quantum implicit steps, positive Bose, U1/energy and independent KL success.
The full Scientific verification 2-version suite can lag; check newest run status before claiming it is green.

## Numerically stable quantum detailed balance at every actual on-shell event

For both conjugate 1↔2 and the two allowed 2↔2 pair Yukawa channels, physical quantum FD/BE weights F,R satisfy exactly ln(F/R)=A=S·η, the on-shell chemical affinity. Direct floating subtraction F−R loses digits when near thermal equilibrium. The source now evaluates the same signed microevent extent with:
- A≥0: J=F[−expm1(−A)]
- A<0: J=R expm1(A).

This enforces J·A≥0 for every positive physical event (subject to ordinary floating rounding), hence the represented microscopic H=ΣW J A≥0, without altering the true common event E1/E2/φ four-leg geometry. Analytic sparse source derivatives still derive independently from F,R; the source and implicit Jacobian retain finite-difference verification. Tests check almost-balanced ±1e−12 scalar chemical potentials and exact common-U1 equilibrium.

**Scope:** reliable *conditional* vacuum-mass, finite-parent-cutoff 5-population dark-sector discrete strong roots at tested parameters. Neither positivity nor nonlinear root existence is unconditional; no actual Γthermal, scalar self-interaction phase, finite SM heat, FRW frozen f(q), absolute Ω, physical BEC or CLASS-ready claim. Additionally, nearcritical massive Bose scalar physical susceptibility needs a special IR quadrature, separately examined in PHI-BOSE-CRITICAL-IR-SUSCEPTIBILITY-AND-FIXED-GAUSS-GATE-2026-11-24.md.

Reproduce from kinetics/phi_bose_five_dirac_full_yukawa_implicit.py and tests/test_phi_bose_five_dirac_full_yukawa_implicit.py, tests/test_phi_bose_five_dirac_kl_certificate.py using OPENBLAS_NUM_THREADS=1 python -m pytest -q on the two tests. NEXT_ACTION: stringent 3-grid/quadrature/dt and nearcritical positive physical moment convergence, then realistic SM+FRW.
