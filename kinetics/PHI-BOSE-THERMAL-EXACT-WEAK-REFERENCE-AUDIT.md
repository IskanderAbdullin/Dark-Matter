# Exact angular Bose/FD thermal source primitive: independent physical Yukawa weak benchmark (2026-11-23)

## The bottleneck and new analytical reduction

After proving the hard-cutoff scalar low-energy support floor and a physically constrained toy exact-shell graph, a positive normalized on-shell **weak** collision rate is required as an independent benchmark before assigning any dynamical finite-volume phase capacities. This stage does **not** claim a physically accurate strong collision/time scheme. For m1=.1,m2=.02,mφ=86×10⁻⁹ GeV, g1=4, Γrest from tree scalar Yukawa y=.001, take physical isotropic FD mother T1=.9 GeV, FD daughter χ2 and BE scalar both T0=1 GeV, chemical potentials zero.

For any physical two-body shell E1=E2+Eφ, the **ratio** of exact quantum forward/reverse microevent factors is

`F/R = exp[A], A=−E1/T1+E2/T0+Eφ/T0 = E1(1/T0−1/T1).`

Consequently A depends only on the **parent** energy and is negative when T1<T0: the net reaction destroys φ and entropy production per microevent `(F−R) A ≥0`. The same common-T chemical equilibrium T1=T0 has exact F=R for every angle (not a post hoc equilibrium reset).

With daughter scalar energy k∈[Emin,Emax] from the Lorentz-boosted shell and T2=Tφ=T0, the inverse occupation integral is **analytically elementary**:

`I(E1;Emin,Emax)=∫ dk f_FD(E1−k;T0) b_BE(k;T0)`

`= T0 f_FD(E1;T0) [ L(Emax)−L(Emin) ],`

`L(k)=ln(expm1(k/T0))−ln(exp(k/T0)+exp(E1/T0)).`

This follows by setting t=exp(k/T0), A0=exp(E1/T0), so the occupation product times dk is `T0 dt/[(t−1)(t+A0)]`. The logs are evaluated as `ln(expm1 z)=z+ln(1−exp(−z))` and `logaddexp(z,E1/T0)` for stability; the [exact kinematic E_min](PHI-BOSE-PARENT-CUTOFF-KINEMATIC-IR-AUDIT.md) avoids catastrophic subtractive cancellation.

The **entire thermal net number source** is then only one parent-momentum integral:

`Rφ=∫₀^pmax dp [g1 p²/(2π²) Γrest m1/E1] [(1−f1) expm1(A)/(2b)] I(E1;Emin,Emax)`

with b=p p*/m1, and a finite p→0 limit. The 1D integral and its differential are independently coded using scipy adaptive quadrature, **not** by calling the existing mapped event source. The same script independently builds a positive Gauss parent×boosted log-Eφ angular event measure, deposits the physical on-shell parent/daughter energies and microscopic entropy and checks all three energy moments.

## Independent numerical results and cutoff bias

| pmax (GeV) | Exact angle + 1D adaptive parent Rφ | Independently sampled positive event Rφ (nr96,angle48) | Event Qφ | Mean daughter Eφ |
|---:|---:|---:|---:|---:|
| 2 | −1.295453955495×10⁻¹¹ | −1.295453955492×10⁻¹¹ | −2.670878687858×10⁻¹² | 0.2061732 GeV |
| 8 | −3.212015339411×10⁻¹¹ | −3.212015333089×10⁻¹¹ | −1.394186043243×10⁻¹¹ | 0.4340534 GeV |
| 16 | −3.232104327009×10⁻¹¹ | −3.232104320686×10⁻¹¹ | −1.438236797986×10⁻¹¹ | 0.4449846 GeV |

Relative source errors vs **independent analytical-angular** reference at nr96/angle48: about −2.6×10⁻¹² / −2.0×10⁻⁹ / −2.0×10⁻⁹. The larger pmax8/16 difference is consistent with finite parent integral/sampling tolerances, not a new physical effect. Three (nr,angle) resolutions (24,12), (48,24), (96,48) all show positive event weights, Qφ negative, two daughter energies negative, parent energy positive, H>0, and total **on-shell** Q1+Q2+Qφ <4×10⁻²⁶ in absolute units. These exact weak ledgers DO NOT certify a finite capacity strong source.

Importantly pmax2 omits about 59.9% of the **absolute** pmax16 rate (1.295/3.232≈.401), whereas pmax8 vs16 differs by <0.63%. The signed φ energy moment and mean Eφ are particularly cutoff sensitive: pmax2 mean≈.206 GeV vs pmax16≈.445 GeV. Thus reproducing just one total number rate at one parent cutoff cannot be mistaken for physical spectral convergence. This reference complements the [boosted log-angle audit](PHI-BOSE-DAUGHTER-LOG-ANGLE-AUDIT.md) and gives a lower-dimensional independent gold-standard test.

## Reproducibility and limitations

- [Standalone original code](phi_bose_thermal_exact_weak_reference.py), [three focused regression tests](../tests/test_phi_bose_thermal_exact_weak_reference.py), [machine data](results/phi-bose-thermal-exact-weak-reference-2026-11-23.json), [dated stage](../archive/stages/2026/11/2026-11-23-phi-bose-thermal-exact-weak-reference.md).
- Commands: `OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_phi_bose_thermal_exact_weak_reference.py` and `OPENBLAS_NUM_THREADS=1 python -m kinetics.phi_bose_thermal_exact_weak_reference`.
- The energy-integral analytic primitive was independently checked against direct adaptive daughter integration for p=.1/1/2/8/16; initially excessive independent `epsabs=1e−9` produced an **artificial** ~10⁻⁷ relative discrepancy at a tiny integral and was corrected to `epsabs=1e−15`. All **three local focused tests passed** (~1.32 s) after choosing an appropriate independent numerical tolerance.

This is a kinetic benchmark for a *prescribed* χ1/χ2/φ thermal population. No physical strong solver, gauge-complete portal, in-medium rates, finite SM, FRW, zero-mode/BEC or 86-eV relic/CLASS/observational likelihood. The next task is to turn this exact physically normalized **weak event measure** into a **shared positive energy-resolving finite phase-space capacity strong operator**, auditing unavoidable lattice resonance/disconnected graph issues and time+IR convergence.


## Continued independent physical Qφ energy reference: a real-dilogarithm primitive

Beyond the previously exact number primitive, let `t=exp(k/T0)`, `A0=exp(E1/T0)`. The positive angular Bose×FD energy moment is
`I1=∫_{Emin}^{Emax} dk k fFD(E1−k;T0)bBE(k;T0)=T0² fFD(E1;T0)[P(Emax)−P(Emin)]`,
where `P(k)=−Li₂(1−t)−ln(t)ln(1+t/A0)−Li₂(−t/A0)` using real dilogarithms. This follows directly by integrating `ln(t)/[(t−1)(t+A0)]` and is evaluated stably with `scipy.special.spence(z)=Li₂(1−z)`.

Consequently the separately integrated **physical scalar energy weak moment** is
`Qφ=∫₀^pmax dp [g1 p²/(2π²) Γrest m1/E1][(1−f1)expm1(E1(1/T0−1/T1))/(2b)] I1`.
For parent pmax=2,8,16 GeV the exact 1D Qφ values are respectively **−2.670878687857882e−12, −1.3941860432430647e−11, −1.4382367979860587e−11**. The completely independent 2D positive on-shell weak event quadrature at nrad96/angular48 agrees to relative **9.1e−15, 7.1e−14 and 1.93e−13**. Direct adaptive daughter-energy integrals at parent p=.1/1/2/8/16 agree with the Li₂ primitive to ~1e−14.

Two additional focused regressions have passed locally; **5/5 tests in 0.58 s** after extension. Previous source/JSON/stage are preserved and extended. This establishes **independent number and energy weak benchmarks**, not a coupled physical strong finite-mode time evolution, zero-mode, freezeout, Ω or CLASS.
