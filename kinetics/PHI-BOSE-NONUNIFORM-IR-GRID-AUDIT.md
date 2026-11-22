# Nonuniform inverse-fugacity Bose IR energy hats: a major source improvement and a new near-critical moment-time obstruction (2026-11-22)

## Falsifiable question

The physically normalized moving χ1↔χ2+real φ Yukawa 1↔2 source is mathematically microreversible and the Bose p²dp entropy-dual moment representation eliminated a false first-energy-cell condensate. The [independent IR shape audit](PHI-BOSE-IR-SHAPE-GATE-AUDIT.md) nonetheless found that **uniform E hats grossly alias a legal curved Bose inverse fugacity**: a synthetic IR bump of width 0.02 GeV gives scalar source error +68% at N16 and −11.35% at N64, despite perfect affine thermal equilibrium. The previous implicit scalar time code also needs a feasible inverse-fugacity representation near ηφ(E=mφ)=mφ/T≈8.6×10⁻⁸ at mφ=86eV,T=1GeV.

Before testing, compare four **predeclared nonuniform energy-hat families** at identical domain, parent p/c quadrature and χ1/χ2 fermion grids: uniform E, quadratic concentration E=m+(Emax−m)x², cubic E=m+(Emax−m)x³ and two-zone with half of the bins in [m,0.25GeV]. Compute exact parent-on-shell φ net number source and true micro-entropy versus the same independent continuous-occupation weak reference for two prescribed nonthermal φ shapes at N16/32/64/128/256. Thermal affine η(E)=E must remain exactly reconstructed on ALL meshes. Also attempt a genuine **finite-step implicit** moment evolution with the physical tree Yukawa event weights and strict ηφ>0, recording solver failure rather than labelling it physical BEC.

## Mathematical conservation and physical scope

For any strictly ordered E_i beginning at the physical scalar mass, hats h_i(E) are nonnegative, form a partition of unity and satisfy ΣE_i h_i(E)=E. On nonuniform grids these identities remain **exact**: no change to the microscopic Fermi/Bose forward-minus-reverse source, daughter boost, or physical radial p²dp moment integration is needed. Moments M_i=∫p²dp h_i(E) b[Σh_j(E)η_j] and event deposit with the SAME hats still preserve quantum numbers and energy, with entropy gradient η_i and exact H=ΣW(F−R)ln(F/R)≥0. A thermal η=E/T has no interpolation error irrespective of node placement. The p-space Gaussian integral remains in physical momentum, not a forced energy-node Bose rest sample.

Assumed *illustrative* masses and couplings: mχ1=0.1GeV, mχ2=0.02GeV, mφ=86e−9GeV, y=.001, gχ1=4, Tχ1=.9GeV,Tχ2=Tφ=1GeV, physical parent pmax2GeV with nr24,nangle48. All nonthermal distributions and IR bump locations are chosen test fixtures, **not** a realistic FRW/SM evolution. The two bosonic inverse fugacity fixtures are ηbroad(E)=E+.30exp[−((E−.10)/.06)²] and ηIR(E)=E+.40exp[−((E−.04)/.02)²], both strictly positive on E≥mφ. The independent baseline is the already-existing moving Yukawa continuous-η source, not the same Galerkin code being evaluated.

## Independent quantified grid improvement

For the sharp IR η, percent error of net scalar rate Rφ relative to the same physically weighted independent parent weak integral:

| Scalar hat intervals | Uniform E | Quadratic E(x) | Cubic E(x) | Two-zone |
|---:|---:|---:|---:|---:|
| 16 | +68.146% | −5.144% | −6.827% | −11.704% |
| 32 | +29.567% | −3.357% | −2.507% | −5.243% |
| 64 | **−11.351%** | −1.004% | **−0.627%** | −1.594% |
| 128 | −6.202% | −0.274% | **−0.163%** | −0.454% |
| 256 | −1.988% | −0.076% | **−0.0430%** | −0.104% |

At N64 the cubic IR mesh reduces absolute source bias from 11.35% to 0.627%, an approximately **18-fold reduction with the same total modes**, and at N256 from 1.99% to 0.043% (~46×). For this specific fixture the cubic family is an effective local numerical accelerator, not a universal physically optimal mesh; this must be tested against the **dynamically produced** unknown φ spectrum. Nonlinear broad η at N64 has uniform R error −4.797%, two-zone −0.317%, and microscopic H error −4.168% vs −0.263%. Some methods can be worse at N16 because moving bins changes representation elsewhere; this is also preserved.

For thermal affine ηφ(E)=E, **all four mesh families reproduce both Yukawa net weak source Rφ and H exactly to relative rounding ≤10⁻¹²** at N32. The positive physical Bose number/energy measure is unchanged. Every curved test has negative net φ rate and positive microscopic entropy with exactly paired N/E on the same event graph; the source is still *conditional* on finite parent cutoff pmax2, itself previously shown to omit ~65.6% of the pmax16 inverse net rate.

## New negative: nonuniform near-critical Bose inverse moment map is not automatically time-integrable

Under the *physical tree width* and moving-mother p/c parent source, start with common-temperature ηφ(E)=E/T (zero φ chemical potential, ηφ(mφ)≈8.6×10⁻⁸), but χ1 colder than the common χ2/φ bath, so the microscopic φ source is net negative **at every event**. For N16 on the quadratic φ grid, the first scalar hat source is **exactly zero** and all other φ hat sources are ≤0; therefore there is **no positive microscopic source in the first cell** from which to infer genuine condensate formation.

Nevertheless, the fully implicit physical moment Newton solver has a strict feasibility problem:

| φ grid, dt in illustrative GeV⁻¹ | Result | Physical entropy |
|---|---|---|
| N16 uniform, dt=10³ | converged; ΔNφ≈−1.16484×10⁻⁸ | increases |
| N16 quadratic, dt=1 | converged | increases |
| N16 quadratic, dt=20 | converged, min ηφ≈2.604×10⁻⁸ | increases |
| **N16 quadratic, dt=30** | **Newton cannot decrease residual while maintaining ηφ>0** | no state accepted |
| **N16 quadratic, dt=10³** | **same fail-closed rejection** | no state accepted |
| N32 cubic, dt=10³ | converged, min ηφ≈4.719×10⁻⁸ | increases |
| N64 cubic, dt=10³ | converged, min ηφ≈8.608×10⁻⁸ | increases |

The failure occurs with **monotonically negative microscopic scalar sources**, showing that a physically valid sign of production and an H theorem alone do NOT ensure the nonlinear finite-dimensional inverse map admits an admissible endpoint for arbitrary mesh/time. The reason is the coupled physical moment susceptibility between adjacent hats: a zero source for the first integrated hat can require a negative change in its rest-node η to compensate changes at its neighbor. At N16 quadratic a linearized physical moment sensitivity gives first-node dηφ/dt≈−3.0×10⁻⁹ in this conditional normalization, while ηφ(mφ) starts only 8.6×10⁻⁸. A large dt therefore reaches the ηφ=0 boundary even though there is no physical Bose delta-mode or first-cell production in the model. A smaller dt≤20 succeeds, but accepting that result does not solve the general critical-boundary problem.

**This is a numerical realisability/IR moment-manifold obstruction, NOT evidence of physical Bose condensation**, and must not be silently resolved by clipping ηφ or inventing δ³(p) bosons. Physical BEC requires a separately derived and normalized zero-mode kinetic equation, threshold and qmin→0 tests.

## Source, tests, and next decision

- [Executable nonuniform E meshes and on-shell rate/time controls](phi_bose_nonuniform_ir_grid.py).
- [Six focused regressions](../tests/test_phi_bose_nonuniform_ir_grid.py): thermal reproduction, quantitative narrow/broad IR improvements, explicit fail-closed Newton and proof of nonpositive microscopic φ deposition.
- [Machine-readable family scan and time failures](results/phi-bose-nonuniform-ir-grid-2026-11-22.json).
- [Permanent stage evidence](../archive/stages/2026/11/2026-11-22-phi-bose-nonuniform-ir-grid.md).

    OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_phi_bose_nonuniform_ir_grid.py
    OPENBLAS_NUM_THREADS=1 python -m kinetics.phi_bose_nonuniform_ir_grid

Six local tests pass. No green GH Scientific CI is asserted until the stage commit is independently checked.

**P0 NEXT:** derive a physically admissible **IR zero-flux and entropy-dual moment-space time strategy** to handle near-critical ηφ→0 without fake BEC; inspect the source moment image/tangent cone and distinguish a grid-artificial boundary from genuine condensate physics. Independently repeat physical scalar full weak + nonlinear time on ≥3 IR meshes, **parent pmax12/16** with resolved p/c and 3+ dt, and couple true χ1χ2–SM Pauli strong/e±γν/Z′/FRW only after matching physical charges/energy/H. No physical Γelastic, ΩDM, BEC, frozen fφ(q) or CLASS-ready claim.
