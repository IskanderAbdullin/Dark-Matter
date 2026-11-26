# A false 86-eV Bose first-bin “condensate” from energy-node collocation (2026-11-22)

## Falsifiable physics question

The current coupled χ1↔χ2+real φ kinetic pilot stores a boson occupation bφ(E_i) at uniform **energy nodes**, multiplies it by an exactly integrated radial phase-space cell capacity, and obtains the total boson number as Σ w_i b_i. In particular its initial thermal Bose function bφ(E)=1/[exp((E−μφ)/T)−1] is evaluated at **E0=mφ**. That node has nonzero occupancy (formally divergent for mφ→0) but represents a whole energy cell with strictly positive measure. Is the resulting Nφ first-bin spike a physical condensate or a collision-free quadrature artifact?

A true condensate, if present, requires an independently normalized zero-momentum measure proportional to a delta function and a separate kinetic equation; the continuous thermal excited-particle phase-space density itself can have b(p=0)→∞ and **still yield a finite total number**. See [Blaizot et al., Kinetic approach to a relativistic Bose–Einstein condensate](https://arxiv.org/abs/1510.04552), which explicitly separates the condensate and non-condensate distributions. We do not copy numerical code from that article.

Predeclared gate: compare the **existing first-node collocation convention** with independent high-accuracy physical Bose radial p²dp integrals inside that exact same first cell and over the full finite domain, across at least five positive masses and three energy-grid resolutions. Verify massless equilibrium against the exact zeta(3) integral. Do not infer a condensate from either a large node occupation or a number that varies with mφ/grid.

## Mathematical derivation

For a real scalar with degeneracy gφ and μφ=0 at temperature T in natural units, E(p)=sqrt(p²+m²) and excited (noncondensate) number density is

    n_exc(m,T;pmax)=gφ/(2π²) ∫_0^pmax p² dp/[exp(E(p)/T)−1].

With first energy-node E0=m and grid spacing h, the physical first cell spans E∈[m,m+h/2], p∈[0,sqrt((m+h/2)²−m²)]. The current code's exact radial cell capacity is

    w0 = gφ/(6π²) [((m+h/2)²−m²)^(3/2)].

But the **collocated first-cell number** is n0,node=w0/[exp(m/T)−1], whereas the correct cell integral is n0,true=gφ/(2π²)∫_0^p_hi p² dp/[exp(sqrt(p²+m²)/T)−1].

For m≪h≪T,
    w0≈gφ h³/(48π²), b(E0)≈T/m,
    n0,node≈gφ Th³/(48π² m),
    n0,true≈gφ Th²/(16π²), and therefore
    n0,node/n0,true ≈ h/(3m) → ∞ as m→0 at fixed h.

The actual continuous thermal massless Bose number is **finite**, n_exc(0,T;∞)=gφ ζ(3)T³/π². It is a separate mathematical object from a delta-function condensate. The pathology is caused by treating a divergent **point occupation at p=0** as typical of a whole positive-volume cell; the integrated radial measure itself is correct. A massless node has b(E0)=∞ and must NOT be silently regularized and called physical occupancy.

## Independent numerical experiment, with exactly declared illustrative units

Use gφ=1, T=1 GeV and energy endpoint Emax=12.1 GeV in a standalone thermal diagnostic. This is an **illustrative choice, not an FRW temperature history**. Masses mφ∈{.1,.01,.001,.0001,8.6×10⁻⁸} GeV; the last equals **86 eV**. N=16,64,256 uniform E intervals. Independent adaptive SciPy quadrature integrates p²dp/expm1(sqrt(p²+m²)/T) to εabs=2×10⁻¹³, εrel=2×10⁻¹¹. The original project's physical w_i from density_cells is reused as the collocated measure, while the Bose cell integral is calculated separately.

| mφ in GeV | energy cells N | first-cell n_node/n_true | total n_node/n_true |
|---:|---:|---:|---:|
| 0.1 | 64 | 1.621 | 1.00499 |
| 0.01 | 64 | 7.319 | 1.01669 |
| 0.001 | 16 | **287.19** | **8.554** |
| 0.001 | 64 | **65.75** | 1.12237 |
| 0.001 | 256 | 16.63 | 1.00226 |
| 0.0001 | 64 | **651.14** | **2.17694** |
| 86 eV = 8.6×10⁻⁸ | 16 | **3,331,266** | **87,199** |
| 86 eV = 8.6×10⁻⁸ | 64 | **756,352** | **1,363.47** |
| 86 eV = 8.6×10⁻⁸ | 256 | **184,651** | **22.289** |

For the illustrative 86-eV case with 64 intervals, **the first cell alone** contains n_node≈165.873845 in GeV³ units from node collocation, while the exact physical Bose number in that same cell is only ≈0.000219308 (ratio 756,352). The *entire* physically integrated thermal number is only ≈0.1217: the collocated first bin alone exceeds the entire real gas by three orders of magnitude. No exotic dynamics or true BEC was introduced; this is a pure discretization false positive.

As an independent massless limit check, adaptive p-space quadrature over 0≤p≤50T produces 0.12179382823357311 T³ for g=1, while ζ(3)/π²=0.12179382823357308. This confirms the true thermal phase-space number stays finite as mφ→0, unlike the first-node method.

## Partial numerical remedy and why it is NOT yet a replacement collision scheme

One can correct the **initial thermal zeroth-number moment only** by replacing the first occupied cell's number w0 b(E0) by its exact integral, or equivalently using the *number-weighted first-cell average* bbar0=n0,true/w0. Doing **only that** to the N=64 thermal example yields total-number error ≈+0.33995% for the 86-eV example, instead of +136,247% (total ratio 1363.47). At N=16,64,256 with 86 eV, the corrected total ratios are 1.027080,1.003399,1.000320. This implicates the first cell as the dominant source of false total number and shows a practical *diagnostic* improvement.

**However:** bbar0 is a cell average, not b(E0). In the nonlinear on-shell Pauli/Bose Yukawa operator the existing interpolated Bose inverse fugacity ηφ(E)=ln[(1+bφ(E))/bφ(E)] is treated as a nodal energy value; inserting a thermal cell-average into that interpolation is NOT guaranteed to preserve chemical detailed balance, physical Bose phase-space transition rates, quantum H or cell-level energy moments. We have **not** made this unverified modification to the production solver. A physically justified IR-resolved collision discretization (p-grid/thermal cell measures, separate zero mode where needed and exactly reversible event interpolation) must be independently derived and validated before evolving an 86-eV scalar.

## Reproduce, limitations, and next gate

- [Source with exact Bose p²dp integrals and original cell-volume comparison](phi_ir_first_cell_artifact.py).
- [Five focused tests](../tests/test_phi_ir_first_cell_artifact.py): passed locally in 0.60s, no Monte Carlo or external data.
- [Complete 15-case machine table](results/phi-ir-first-cell-artifact-2026-11-22.json), [dated permanent stage](../archive/stages/2026/11/2026-11-22-phi-ir-first-cell-artifact.md).

Commands:
    OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_phi_ir_first_cell_artifact.py
    OPENBLAS_NUM_THREADS=1 python -m kinetics.phi_ir_first_cell_artifact --full

**Scientific status:** A reproducible negative finding about the **existing finite-cell energy-node Bose approximation**, not a no-go theorem for scalar dark matter or for physical BEC. The ratio's mass dependence is an **algorithmic pathology** that could produce false cosmological densities and false “condensate” interpretations if transplanted to the 86-eV regime. No physical Γ_el, late Ωφ, genuine zero mode, finite-SM plasma, H/FRW, frozen fφ(q) or CLASS-ready result. Next stage must replace/validate the IR measure and maintain exact quantum 1↔2 gains, energy/stoichiometry and entropy under a full coupled 4-population evolution; never rescue it by a silent numerical floor.
