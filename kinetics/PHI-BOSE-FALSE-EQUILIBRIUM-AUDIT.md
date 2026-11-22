# Exact chemical detailed balance is compatible with a false 86-eV boson abundance (2026-11-22)

## Scientific hypothesis and independent numerical criterion

The preceding [first-cell Bose audit](PHI-IR-FIRST-CELL-ARTIFACT-AUDIT.md) identified a false scalar count from assigning the **point** BE occupation at rest energy E=mφ to the entire first positive-volume energy cell. The current four-population χ1,χ2,real φ,one lepton nonlinear Pauli/Bose solver has careful microscopic detailed balance and a quantum H theorem. **Does this green equilibrium test certify the absolute scalar number?** It need not: detailed balance constrains *ratios of forward/inverse event probabilities*, not quadrature accuracy of an integrated density.

Before running: set identical Tχ1=Tχ2=Tφ=Tlepton=1 in a closed synthetic network, zero chemical potentials (μ1=μ2+μφ=0), elastic G=0 to isolate the moving real-φ Yukawa source, physical m1=.1, m2=.02, two scalar masses mφ=.001 and mφ=86e−9 GeV. Three uniform-energy meshes N16,64,256 with common Eupper12.1. Require nearly zero Rφ, entropy, charge and energy residuals, while **independently integrate** the physical scalar Bose p²dp thermal number over exactly the same energy interval. All mass, temperature and cutoff units are illustrative, not a cosmological history.

## Formulas, computation and independent controls

For each on-shell moving χ1↔χ2+φ event, F=fχ1(1−fχ2)(1+bφ), R=(1−fχ1)fχ2 bφ. At common T with μ1=μ2+μφ, energy conservation E1=E2+Eφ implies F/R=1 exactly. The code interpolates **fermion log-odds and the boson inverse fugacity ηφ(E)=ln[(1+b)/b]** with the same positive barycentric weights as event deposition, so it reproduces the identity on arbitrary finite E grids. The scalar number is separately computed as nφ,node=Σ_i w_i bφ(E_i), w_i=gφ/(6π²)[p(U_i)³−p(L_i)³].

The genuinely physical *excited, noncondensate* thermal number is

    nφ,true = gφ/(2π²) ∫_0^pmax dp p² / [exp(sqrt(p²+mφ²)/T)−1].

The preceding proof predicts nφ,node/nφ,true grows spuriously as mφ→0 at fixed cell width. The two calculations are logically independent: one verifies the local collision fixed point and the other the density measure. Their agreement is **not required** by microscopic balance.

New [independently reproducible driver](phi_bose_false_equilibrium_gate.py), [two tests](../tests/test_phi_bose_false_equilibrium_gate.py), [six-case raw JSON](results/phi-bose-false-equilibrium-2026-11-22.json), [permanent stage](../archive/stages/2026/11/2026-11-22-phi-bose-false-equilibrium.md). Commands:

    OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_phi_bose_false_equilibrium_gate.py
    OPENBLAS_NUM_THREADS=1 python -m kinetics.phi_bose_false_equilibrium_gate

**Two targeted local tests passed**. Real GitHub Scientific CI for the new commit is to be checked separately. No random seed, no observational external file, no actual FRW temperature assignment or particle-relic calculation.

## Actual results and negative result

In **all six** coupled-four-population stationary runs, Yukawa net scalar production |Rφ| is below 5×10⁻²⁶, absolute microscopic collision entropy below 4×10⁻⁴¹, total energy residual below 3×10⁻⁴¹ and combined charge residual below 2×10⁻⁴¹, while the **absolute scalar number is badly wrong**. These zero residuals are not post-hoc clipping.

| mφ | energy intervals | scalar density from discrete ledger | true excited thermal scalar density | false density ratio |
|---|---:|---:|---:|---:|
| .001 GeV | 16 | 1.041367 | .12174501 | 8.554× |
| .001 GeV | 64 | .1366430 | .12174501 | 1.1224× |
| .001 GeV | 256 | .1220201 | .12174501 | 1.00226× |
| **86 eV** | 16 | **10616.0262** | .12174521 | **87199×** |
| **86 eV** | 64 | **165.995785** | .12174521 | **1363.47×** |
| **86 eV** | 256 | **2.713570** | .12174521 | **22.289×** |

At the 86-eV/N64 fixed point the event source gives |Rφ|=1.7892×10⁻²⁶, H≈1.2563×10⁻⁴¹, energy drift≈7.17×10⁻⁴²; the scalar ledger reports Nφ=165.9958 while the true thermal gas has Nφ≈.1217452 (all density quantities in synthetic GeV³ units). **A perfectly stationary, microreversible and entropy-respecting source can therefore overstate the scalar number by 1363×.**

This is a **numerical false positive for cosmological abundance or condensate inference**, not a physical exclusion of the real scalar or a BEC phase transition. A true condensate is a separate momentum-space delta component, absent from both tested distributions. Thermal detailed balance alone is insufficient to determine correct φ number, physical Γelastic, ΩDM or frozen f(q).

## Next scientific gate

The earlier [first-cell integration analysis](PHI-IR-FIRST-CELL-ARTIFACT-AUDIT.md) showed replacing only the number contribution of the first cell by its exact p²dp integral nearly repairs the *initial thermal number*, but that is NOT a valid drop-in for the dynamic nodal Bose occupation: its quantum inverse fugacity would become inconsistent with on-shell F/R, and energy/entropy require independent verification. Next construct a **moment-consistent φ IR basis** (e.g. momentum-grid physically averaged occupation plus properly matched event projection) with independently verified Bose H theorem, energy, separate zero mode if physically required, and 3+ IR resolutions. Preserve mφ=86e−9 GeV and T=1 GeV as an *illustrative diagnostic*, not a validated cosmological benchmark. Finite SM bath, Γel/H, cosmological Ω, BEC and CLASS-ready all remain unproven.
