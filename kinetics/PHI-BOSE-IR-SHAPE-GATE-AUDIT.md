# Independent Bose IR mode-curvature and mother-cutoff weak Yukawa audit (2026-11-22)

## Question, method, and model boundaries

The new physical p²dp entropy-dual Bose/Fermi moment representation and implicit χ1↔χ2+φ solver avoid the false 86-eV rest-node condensate and preserve quantum N/E/H. But its physically normalized synthetic time result still differs ~9.7% between N8 and N32 energy-mode grids. To disentangle possible causes, compare three independent gates: **Bose log-fugacity interpolation**, **parent p/c quadrature**, and **parent kinetic cutoff**. Exact detailed balance, positive entropy and energy conservation alone do not measure any of these errors.

Assumed illustrative masses mχ1=0.1 GeV, mχ2=0.02 GeV, mφ=86e−9 GeV, hypothetical y=0.001 and gχ1=4. Physical tree-level vacuum Γrest≈2.7501974166e−9 GeV. All events use exact boosted E1=E2+Eφ and full finite Bose/Pauli quantum source F=f1(1−f2)(1+bφ), R=(1−f1)f2 bφ, signed J=F−R and H=∫W J ln(F/R)≥0. The mother weight includes g1 p²dp/(2π²) × Γrest m1/E1 × dc/2. There is no finite SM e±γν/Z′ plasma, real FRW/H evolution, experimental coupling fit or physical condensate; 86 eV is a **historical illustrative mass**, not a predicted mass or observed DM signal.

The independent direct continuous-occupation source is the existing moving_yukawa_bose_pauli_weak.weak_moments. A second route evaluates η(E) at energy-hat nodes and reconstructs the same physical mother/daughter quantum events using the physical p-space Galerkin source. At parent pmax=2 GeV, nr24 and nangle48, refine Bose energy-hat intervals N16,32,64,128,256,512. Two nonthermal physical Bose inverse-fugacity shapes are tested against the affine thermal control:

    thermal:     ηφ(E) = E
    broad:       ηφ(E) = E + 0.30 exp[−((E−.10)/.06)²]
    IR-curved:   ηφ(E) = E + 0.40 exp[−((E−.04)/.02)²]

All are positive Bose inverse fugacities on the physical E≥mφ domain and all are **prescribed synthetic distributions**, not freeze-out solutions.

## Numerical findings: conditional Bose spectral alias

For the thermal affine ηφ(E)=E, independent weak Rφ and H agree with all N16–512 Galerkin hat calculations to ≤8×10⁻¹⁵ relative, validating exact affine reproduction; **this does not certify nonthermal spectra**.

For the sharply IR-curved shape, the original analytic-occupation weak reference at parent pmax2 is Rφ=−5.9510246418e−12, H=+1.0395604380e−12. The measured errors of Galerkin mode interpolation at identical parent quadrature are:

| E intervals | Relative Rφ error | Relative H error |
|---:|---:|---:|
| 16 | +68.1464% | +35.9721% |
| 32 | +29.5671% | +4.1296% |
| 64 | −11.3510% | −10.3124% |
| 128 | −6.2019% | −5.9463% |
| 256 | −1.9875% | −1.9049% |
| 512 | −0.5735% | −0.5569% |

The signed interpolation error changes sign between N32 and N64. All net rates remain negative and H positive; sign changes in *relative error* do not imply physical failure of the dark matter hypothesis. At N512 the shape has not reached a subper-mille source accuracy. An early two-point extrapolation would be misleading.

For the broader nonthermal η the reference is Rφ=−3.774816970e−12, H=+8.295360996e−13. At N16 the **number-source error is only +0.316% while H error is −16.83%**. Thus a seemingly accurate net production source does not certify the more sensitive entropy functional. At N64/128/256/512 number errors are −4.797/−1.396/−0.354/−0.0942%, entropy errors −4.168/−1.192/−0.309/−0.0835%. The last ratio is consistent with the expected O(h²) interpolation regime on this fixed physical p/c rule.

## Independent incoming quadrature and mother-tail evidence

For the IR-curved shape at pmax2, independently refine mother p/c nr24,nangle48→nr48,nangle96→nr72,nangle144: Rφ=−5.951024642e−12→−5.951959775e−12→−5.951964666e−12. The coarsest-to-finest shift ~0.016% is far below the N16 interpolation error; no statement about other cutoffs or more general distributions follows.

At fixed analytic-occupation IR-curved η and nr72,nangle144, vary only physical mother pmax, leaving **the Bose interpolation mesh entirely out of the test**:

| Parent pmax GeV | Rφ | Microscopic H |
|---:|---:|---:|
| 2 | −5.951964666e−12 | +1.039830393e−12 |
| 4 | −1.393433793e−11 | +3.555904102e−12 |
| 6 | −1.656245507e−11 | +4.937547463e−12 |
| 8 | −1.715801539e−11 | +5.381165733e−12 |
| 12 | −1.729468032e−11 | +5.517769135e−12 |
| 16 | −1.729875243e−11 | +5.523641876e−12 |

**Critical negative:** pmax2 contains only about 34.4% of the magnitude of the pmax16 net φ destruction under the same synthetic initial spectra. The truncated mother tail therefore causes ≈65.6% rate omission, a much bigger uncertainty than the fixed-pmax2 3-grid parent quadrature numerical difference. Parent pmax12→16 moves R by ~0.0235%; this is only a conditional cutoff test, not a fully validated Boltzmann solution.

## Reproduction, tests, decision

- Source: [phi_bose_ir_shape_gate.py](phi_bose_ir_shape_gate.py), using existing [physical moment source](phi_bose_moment_galerkin.py), [implicit solver tree weight](phi_bose_implicit_moment_step.py), and [independent moving weak](moving_yukawa_bose_pauli_weak.py).
- Four focused tests: [test_phi_bose_ir_shape_gate.py](../tests/test_phi_bose_ir_shape_gate.py), passing locally. They check thermal affine exactness, nonlinear curved IR alias, refinement of broader η and major independent parent tail loss.
- Full compact numerical data: [phi-bose-ir-shape-gate-2026-11-22.json](results/phi-bose-ir-shape-gate-2026-11-22.json), and [dated stage](../archive/stages/2026/11/2026-11-22-phi-bose-ir-shape-gate.md).

Reproduction commands (sandbox or fresh repository):
    OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_phi_bose_ir_shape_gate.py
    OPENBLAS_NUM_THREADS=1 python -m kinetics.phi_bose_ir_shape_gate

**Scientific interpretation:** the new physical-moment entropy scheme fixes false scalar number and guards finite-step entropy, but mode shapes and mother cutoff still require independent convergence. An exactly stationary thermal Bose spectrum says nothing about the interpolation of curved nonthermal η; a green CI or conserved H does not certify Ω or Bose condensation.

**Next high-priority gate:** at physically sufficient pmax≥12 independently converge parent p/c, construct IR-adaptive Bose energy or p-grid with ≥3 resolution/time comparisons in the fully nonlinear physical-moment implicit solver, isolate φ number/energy/entropy errors, then only couple with the existing χ1/χ2 shared-lepton strong Pauli, finite SM ν/e±γ, vector mediator and true FRW. Zero occupation b=0 and true δ³(p) condensate remain unresolved. No physical Γelastic, frozen f(q), ΩDM, BEC or CLASS-ready conclusion.
