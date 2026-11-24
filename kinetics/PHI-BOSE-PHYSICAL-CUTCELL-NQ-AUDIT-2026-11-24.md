# Physical Yukawa cut-cell weak N and all three Q: converged shell-crossing audit (2026-11-24)

## Objective and scope

The previous positive invariant decay measure and exact unweighted daughter-cell geometry did not establish a physical daughter-bin-resolved Bose/Pauli **weak** rate, let alone a nonlinear **strong** kinetic solver. The old 12×12 parent-angle quadrature gave a false zero in the scalar first energy cell. This stage computes the conditional on-shell **number and all three physical energy sources**, with nonzero chemical potentials and reverse reaction, and independently checks an exact daughter primitive.

Masses/coupling (GeV): m1=.1, m2=.02, mφ=86e−9, Yukawa y=.001, g1=4. Parent pmax=2,8,16 GeV. Quadratic φ-bin energies E_j=mφ+[hypot(pmax,m1)+.2−mφ](j/N)². All populations are prescribed thermal FD/BE, with parent T1=.9, daughters T2=Tφ=1 and μ=0 for the IR benchmark. Separate chemical tests have T1=T0=1, μ2=.05, μφ=−.2 and μ1=−.09 (forward), −.21 (reverse) or −.15 (equilibrium).

## Two independent routes and a constructive numerical change

**Positive cell-shell route:** intersect each scalar bin [L,U] with actual massive boosted shell [kmin(p),kmax(p)], then integrate k in the clipped interval using positive Gauss weights on the exact log(k) mapping. Evaluate explicit Fermi/Bose forward and reverse occupations at true E1, E2=E1−k and k. Multiply by common positive parent measure [g1 p² dp/(2π²)] Γrest m1/E1 × dk/(2 p kstar/m1). Deposit each microevent into Rφ, Qφ=kJ, Q1=−E1 J, Q2=(E1−k)J. Pointwise Q1+Q2+Qφ=0 and N1+N2, N1+Nφ conserve. The physical scalar V is integrated exactly and energy capacity W by independent positive radial quadrature; these are not substituted for event energies.

**Essential independent parent partition:** for each interior daughter bin boundary below rest-frame Estar, bracket the true root kmin(p)=Eedge, and above Estar bracket kmax(p)=Eedge. The resulting parent segments remove discontinuous cut opening/closing from the middle of a Gauss interval. Coarse-to-fine tests use (Gauss parent order **per smooth segment**, log-angle nodes, φ bins) = (6,8,8), (10,16,16), (16,32,32). This is not cosmetic oversampling.

**Independent analytical reference:** integrate daughter FD×BE number via the exact endpoint logarithm and physical boson energy via a real-dilogarithm primitive; integrate over parent p via separate adaptive scipy.integrate.quad. This does NOT reuse the cut-cell event quadrature or its weights. Further independent direct daughter integration checks both primitives for p=.1,1,2,8,16 GeV.

## Results: all three physical energy moments, 3 cutoffs, 3 meshes

The table reports maximum relative error across Rφ, Qφ, Q1 and Q2 against the **independent analytical** benchmark on the medium mesh (10,16,16), plus the finest pmax16 (16,32,32). The coarse (6,8,8) mesh has maximum discrepancy below 1.44×10^−7 for all three non-equilibrium regimes and all three cutoffs.

| Thermal trial state | P=2 medium | P=8 medium | P=16 medium | P=16 fine |
|---|---:|---:|---:|---:|
| T1=.9, all μ=0 | 3.81e−13 | 7.77e−14 | 7.66e−14 | 7.64e−14 |
| common T, Δμ=+.06 | 2.96e−12 | 9.51e−13 | 9.38e−13 | 9.36e−13 |
| common T, Δμ=−.06 | 3.06e−12 | 9.93e−13 | 9.81e−13 | 9.78e−13 |

For the zero-μ T1=.9 benchmark at P=2 medium:
- Rφ = −1.295453955491667e−11 GeV⁴
- Qφ = −2.670878687857877e−12 GeV⁵
- Q1 = +1.682403591788566e−11 GeV⁵
- Q2 = −1.415315723002778e−11 GeV⁵

Q1+Q2+Qφ closes at floating-point precision; tested non-equilibrium physical entropy rate is positive. At P=8 and P=16, Rφ is respectively −3.212015333089339e−11 and −3.232104320686507e−11: P=2 is not parent-UV converged. For common-T chemical-forward P=16: Rφ=+3.474664979601323e−12, Qφ=+2.448885323000642e−12. Reverse P=16: Rφ=−3.332572619668482e−12, Qφ=−2.330252497766270e−12. Exact chemical equilibrium μ1=μ2+μφ gives all |Rφ| below 5e−27 GeV⁴ and all |Qi| below 8e−27 GeV⁵, a floating-point null, not a physical source.

## NEW negative control: cell DOS mean cannot replace collision-weighted energy

At P=2,Nφ=16 the physical first cell is **actually populated** in the signed inverse process:
- Rφ,first = −2.762149360447387e−12 GeV⁴; Qφ,first = −1.227799761951597e−14 GeV⁵.
- The old sparse 12×12 event graph assigned exactly zero to both.
- New partitioned integration obtains these within 2e−9 relative already at just 6 parent Gauss points per smooth segment ×8 log-angle points. By contrast a single **unsplit** 48-point parent integral with 16 log-angle points gives first-cell Rφ=−2.762068702645739e−12, a false 2.92e−5 relative cell error, even when integrated totals look converged.
- Exact physical DOS cell-average Wφ/Vφ = 0.006452696395 GeV. But the actual collision-weighted first-cell scalar energy Qφ/Rφ≈0.004445096 GeV.

Depositing the signed collision number at constant W/V rather than the event k would predict Qφ,first≈−1.7823311221e−14 instead of the true −1.2277997620e−14 GeV⁵: **45.16% spurious energy-flux magnitude bias**, defect −5.5453136e−15 GeV⁵ in this benchmark. This is a specific, independent **negative result** for naive one-DOF constant-energy cell kinetics, not a no-go theorem for richer energy-resolving discretizations.

At the same P=2, the first N64 φ cell is strictly kinematically forbidden and has exactly zero weak N/Q; it must not be conflated with the false N16 quadrature zero or with a condensate mode.

## Reproduction and limitations

Sources:
- [Standalone implementation](phi_bose_physical_cutcell_weak.py)
- [Five focused regression tests](../tests/test_phi_bose_physical_cutcell_weak.py)
- Test command: OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_phi_bose_physical_cutcell_weak.py
- Full audit command: OPENBLAS_NUM_THREADS=1 python -m kinetics.phi_bose_physical_cutcell_weak

Five focused tests **passed locally** (4.79 seconds, Python/SciPy). Test suite covers the 3×3×3 non-equilibrium N/Q matrix, independent direct daughter integrals, Pauli/Bose equilibrium/reverse sign and microscopic entropy, 12×12 false zero and unsplit false precision, and true N64 kinematic zero. Public Scientific CI must be independently verified for the new commit.

**Scientific boundary:** a positive *weak* scalar-bin rate with exact N/Q ledger is NOT a physically normalized, entropy-stable, finite-capacity **strong** nonlinear χ1/χ2/φ kinetic evolution. Daughter capacity V is physical, but no positive on-shell shared-state transition matrix, nonlinear time discretization, nullspace or continuum nonthermal convergence has been established. The finite P≤16 thermal parent shell cannot access the true 86-eV zero-momentum region; do not infer a condensate, actual relic Ω, frozen fφ(q), CLASS or galaxy profiles. Also no physical SM/plasma Γ(T), finite-temperature masses, FRW pressure work or observational likelihood has been calculated here.

**NEXT_ACTION:** construct a shared physical-capacity *energy-resolving* nonlinear collision operator, matching this N/Q cell reference and chemical equilibrium at fixed occupations, preserving 2 charges, exact energy and entropy under **finite-time bounded** steps. Audit graph rank/spurious invariants, cell energy realizability and qmin/UV/time convergence; retain failure as negative evidence rather than injecting a Bose floor or free thermal bath.
