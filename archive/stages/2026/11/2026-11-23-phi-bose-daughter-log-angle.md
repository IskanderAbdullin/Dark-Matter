# 2026-11-23 — Boosted 86-eV real Bose daughter: exact log-energy angular quadrature and false weak-rate precision

- **ID:** DM-KINETICS-20261123-MOVING-DAUGHTER-LOG-E-ANGLE
- **Active launch start:** 2026-11-23T06:50:13Z UTC, measured sandbox clock; actual overall elapsed recorded separately at safe final checkpoint.
- **Previous research commit:** [54de58d9](https://github.com/IskanderAbdullin/Dark-Matter/commit/54de58d9b44a323d8559ef01608c28e665a27270), public research/kinetic-equilibrium-audit-20261121 draft PR #1. Scientific CI #45128030398 in_progress at last check, do not call it successful before confirmed.
- **Task context:** [AGENTS](../../../../AGENTS.md) → [GOAL](../../../../GOAL.md) → [PLAN](../../../../PLAN.md) → [kinetics physical spec](../../../../kinetics/CONTINUOUS-RESEARCH-SPEC.md) → [guardrails](../../../../kinetics/RESEARCH_GUARDRAILS.md) → [NEXT_ACTION](../../../../kinetics/research-control/NEXT_ACTION.json).

## Predeclared numerical falsification criterion

The boosted Yukawa χ1→χ2+φ daughter has Eφ=a−b cosθ*; thermal bφ=BE(Eφ/T) sharply peaks in the forward CM direction despite affine thermal ηφ(E). Independent mother-radial p convergence does not bound the CM-angle error. Derive and verify a positive log-Eφ change of variables against closed exact I0 and Rayleigh–Jeans I1/E integrals, then full on-shell quantum Pauli/Bose inverse+stimulation Rφ/H at 3+ angular grids and several parent cutoffs and scalar IR shapes. Preserve failures of uniform Gauss and underresolved sharp nonthermal shapes.

## Mathematical and physical method

For fixed mother boost, E_min=a−b>0, E_max=a+b, b=p p*/mχ1. Let y=ln(Eφ/E_min)∈[0,L] with L=ln(E_max/E_min); then c=1−(Eφ−E_min)/b and **dc=Eφ/b dy** (positive reversed orientation). Positive Gauss weights w_j L Eφ_j/(2b) ensure exact continuous transformation. Control integrals: I0=∫dc=2 and I_RJ=∫dc/Eφ=L/b. At mother p=0, return uniform c. No approximate BGK, no occupation/moment clipping, no physically invented zero mode.

The scalar mass is **mφ=86e−9 GeV=86 eV**, χ1=.1GeV, χ2=.02GeV, yYukawa=.001, g1=4 inherited Dirac+anti numeric convention; Tχ1=.9 GeV, Tχ2=Tφ=1 GeV; parent p²dp Γrest m1/E1, nr48 Gauss, nangle8..256. Bose IR optional ηphi(E)=E+.4exp[−((E−.04)/.02)²], positive. These are synthetic test occupations, NOT actual cosmological fφ(q).

## Reproducible source and measured observations

- [Independent angle mapper + physical Yukawa weak integral](../../../../kinetics/phi_bose_daughter_log_angle.py).
- [10 focused regressions](../../../../tests/test_phi_bose_daughter_log_angle.py).
- [Derivation and comparison](../../../../kinetics/PHI-BOSE-DAUGHTER-LOG-ANGLE-AUDIT.md).
- [Raw exact-angle/weak numeric series](../../../../kinetics/results/phi-bose-daughter-log-angle-2026-11-23.json).

    OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_phi_bose_daughter_log_angle.py
    OPENBLAS_NUM_THREADS=1 python -m kinetics.phi_bose_daughter_log_angle

Ten local tests passed in 0.66s on clean sandbox, independent real GH CI pending. Analytic RJ control at eight nodes: at p1/2/8/16GeV uniform-c errors −14.7/−27.5/−46.5/−52.9%, mapped |relative error|≤3.1e−13. Exact ∫dc=2 mapped at p16,n8 to ~1e−7 relative.

**Critical underresolved prior rate:** at pmax2GeV and exact physical y, thermal Rφ old ordinary nangle12 **−1.164929972×10−11**; mapped nangle12 **−1.295453962×10−11**, converged mapped nangle64 **−1.295453955×10−11**. Earlier weak rate magnitude underestimated by **~10.08% despite radial nrad8..32 stability**. At pmax8/16, old nangle12 undercounts mapped nangle64 by **~17.35% / 17.40%**, respectively. This is a conditional numerical correction, not discovery of new cosmic physics. Mapped thermal n12 vs n64 at pmax2 differs only ~5e−9 relative.

**Retained negative:** IR-sharp nonthermal bφ curvature at pmax16 mapped nangle8 Rφ≈−1.87121441e−11, nangle96≈−1.72985519e−11 (≈8.17% relative), although nangle48≈−1.72986696e−11. A logarithmic mapping does not guarantee all nonthermal shapes at eight nodes. Parent pmax2→16 still misses ~65.6% of the IR-sharp inverse scalar rate even with nangle96. No full parent p/nc/q/UV/time budget or true cosmological Γelastic.

## Conclusion / NEXT_ACTION

Mathematically exact y=ln(Eφ/E_min) with positive Jacobian isolates a substantial historical angular alias and improves smooth near-thermal on-shell scalar weak integration. Keep original low-angle results as explicitly low-resolution algorithmic comparisons, do not silently overwrite them. Next combine the **mapped physical event weights** with both continuous-overlap and new discontinuous-cell physical Bose moments, compare at matched mode cost and 3+ scalar IR/time/pmax/angle meshes including nonthermal legal ηφ, investigate interface jumps and baryon/lepton/SM coupled model only after numerical uncertainties resolved. This does not establish Bose condensation, 86-eV ΩDM, frozen f(q), real Γelastic/H or CLASS-ready.
