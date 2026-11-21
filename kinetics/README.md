# First genuine Bose-stimulated gain/loss calculation for dark matter production

**Automated safety of scientific interpretation:** [machine-readable adverse-result summary](results/spurious-invariants-2026-11-21.json) and [CLASS readiness gate](cosmology_readiness_gate.py) explicitly forbid marking the current adaptive f(q) spectra as cosmological predictions. A graph rank of three alone does not establish spectral/IR/FRW/microphysical convergence; missing evidence fails closed.


**NEW REJECTION GATE FOR ADAPTIVE GRIDS:** [the collision-invariant rank audit](ADAPTIVE-INVARIANT-AUDIT.md) proves that apparently conservative adaptive 1↔2 grids can have additional, nonphysical conserved quantities. The 27-node sparse mesh has 14 invariants rather than 3; the 56-node three-block mesh has 4. [Unique-mother adaptive ODE](adaptive_resonant_lattice.py) confirms an extra charge is conserved by the full nonlinear solver. Both illustrative spectra are rejected as cosmological predictions despite correct charge/energy/H checks.


**Dark-matter abundance requirement:** a separate [conditional normalization to present Ωh²](../cosmology/FRW-RELIC-ABUNDANCE-BRIDGE.md) finds nφ/T³≈0.12268 for a 86-eV stable scalar at assumed g*Sprod=55 and Ωh²=.12. At initial mother M/T=.05 this corresponds to ≥33.6% closed-bath depletion. Thus matching abundance alone cannot justify the fixed thermal mother assumption or the cold primordial PSD.


**Closed-reservoir budget:** [necessary number-conservation bound](FRW-RESERVOIR-BUDGET.md) compares the fixed-bath FRW output to the thermal mother inventory at initial x=.05. At κ=.2 production exceeds the full original mother number, so external parent replenishment is unavoidable. At κ=.1, a closed reinterpretation would use ~66% of initial mothers, clearly outside passive fixed-source dynamics.


**CLOSED-LATTICE PROGRESS:** [three-species resonant lattice](CLOSED-THREE-SPECIES.md) now evolves chi1 alongside chi2 and phi, with two independent number invariants, total discrete energy conservation and H-theorem tests. Its uniform grid is too coarse to assert physical cold Bose-burst production; the [next detailed specification](CLOSED-FRW-NEXT-SPEC.md) addresses adaptive discretization and closed FRW.


**COSMOLOGICAL PILOT:** [FRW fixed-thermal-mother evolution](FRW-FIXED-BATH-PILOT.md) integrates the on-shell 2D source in comoving q=p/T, radiation H proportional to T² and a dimensionless coupling κ. Weak κ yields mean frozen q~2.5 versus instantaneous flat-bath q~0.003. The mother reservoir, temperature history and energy budget are not self-consistent; this is not CLASS input.


**NEW ON-SHELL 2D SOURCE:** [two-dimensional on-shell collision operator](ONSHELL-COLLISION-2D.md) removes parent/angular remapping entirely and integrates exact allowed daughter-energy cell areas. The independent continuum vacuum source (0.0164959084596) matches the historic quadrature to ~1.7e−10, with positivity, detailed-balance and H-theorem tests. Finite binning and fixed mother bath remain; this does NOT prove full physical production or an observational DM model.


**NEW ANGULAR-QUADRATURE RESULT:** the previous Gauss–Legendre decay-angle nodes miss the soft daughter-momentum endpoint, making the cold-burst mean momentum angle-grid dependent. [Angular IR audit](ANGULAR-ENDPOINT-IR-AUDIT.md) and [log-daughter-energy integration](LOGQ-ANGULAR-QUADRATURE.md) give the analytic explanation and a substantially better convergent *surrogate* benchmark. These do not solve the physical parent-energy/inverse-kinetic gap.


**ERRATUM 21 Nov 2026:** The separate-bin Bose/Pauli reaction rates break exact *physical* parent-event energy conservation despite the exact geometric barycentric projection. Original event-level detailed balance does not establish discrete-ODE equilibrium. See [the energy defect audit](ENERGY-PROJECTION-DEFECT.md), [the equilibrium audit](EQUILIBRIUM-AND-IR-AUDIT.md), and [the entropy-balanced finite-grid surrogate](ENTROPY-BALANCED-LATTICE.md). The original solver is preserved for reproducibility; its output is NOT a physical CLASS spectrum.


**Status: controlled benchmark, 2026-11-21.** Unlike the preceding ad hoc f_alpha cosmology scans, this test explicitly solves a discretized momentum-dependent quantum-statistical collision source from a *specific decay* process studied by Wen Yin, JHEP 05 (2023) 180 https://doi.org/10.1007/JHEP05(2023)180, eqs. (2.2)–(2.4), (2.9)–(2.11), section 3. However it is NOT yet the full paper's numerical calculation: important simplifying physical assumptions are stated below.

## Reaction and exact checks

Assume \(\chi_1\leftrightarrow\chi_2+\phi\), where \(\chi_1\) is a *thermally maintained* Dirac fermion, \(\chi_2\) another fermion, \(\phi\) a real boson. Choose \(M_1/T=0.1\), \(M_2=0\), neglect \(m_\phi\) during production, and measure time in units \(\Gamma_{\rm rest}^{-1}\). Ignore cosmic expansion during this diagnostic. Fix mother occupation \(f_1(E)=1/(e^{E/T}+1)\) but evolve both daughter distributions \(f_\phi(p,t)\) and \(f_2(p,t)\) from zero.

The exact **quantum** gain–loss combination is
\[
S=f_1(1-f_2)(1+f_\phi)-(1-f_1)f_2f_\phi.
\]
It includes Bose stimulated emission, fermionic Pauli blocking, and inverse decay. The 2-body kinematics are
\[
E_\phi=(E_1+p_1\cos\theta)/2,\quad E_2=E_1-E_\phi.
\]
The source integrates the mother thermal phase-space density and time-dilated decay rate \(M_1\Gamma_{\rm rest}/E_1\) over momentum and isotropic rest-frame decay angle. Two *identical weighted reaction events* are deposited into the boson and fermion daughter momentum bins, preserving their number difference exactly. Paired low/high momentum-cell collision events distribute the reaction into four weighted subevents. Crucially, gain/loss is evaluated using the *occupation of each receiving daughter bin*, rather than an interpolated common occupation: this preserves both a nonnegative Bose population and the fermionic Pauli upper bound at the ODE level while conserving the daughter number difference exactly. The grid still causes residual shape dependence.

An exact event-level detailed-balance test at \(f_1,f_2\) Fermi–Dirac and \(f_\phi\) Bose–Einstein at the SAME T vanishes to roundoff, since \(E_1=E_2+E_\phi\). Event energy conservation is similarly machine-exact (bin-centroid energies are not). Numerically the conservation of \(n_\phi-n_{\chi_2}=0\) remains below \(10^{-8}\) relative during both stimulated and control simulations.

## What is measured

The public [reproducible script](one_to_two_fixed_bath.py) evolves up to \(t\Gamma_{\rm rest}=0.03\). The **control** artificially omits the factor \(1+f_\phi\) while keeping the rest of the gain/loss, so it is explicitly NOT an independent physically consistent collision model.

The benchmark uses three independent grid configurations (100/140/72, 150/180/96, 180/240/120) for boson energy bins / parent momentum nodes / decay angle quadrature nodes. The output includes time series of exact daughter number conservation, finite fermion occupancies, boson number, low-momentum fraction, mean momentum, squared-momentum moment, and gross bosonic enhancement.

**Crucial warning:** preliminary independent local tests find the *integrated* boson number agrees at the ~1% level between grids, and the strong Bose-enhanced production and very cold number fraction are qualitatively robust. But **the detailed momentum-distribution shape / mean momentum and peak occupation are not grid converged**; different grids exhibit a 20%-level difference in mean momentum and much larger peak changes. Therefore this distribution must NOT yet be supplied to CLASS or interpreted as a prediction for 86-eV DM. The code reports this limitation explicitly and does not falsely mark the shape as converged.

## What remains to turn this into a physical production model

1. Replace fixed bath by explicitly budgeted mother-state kinetics or a physically justified fast bath thermalization model.
2. Implement the expanding FRW Universe through comoving momentum and temperature-dependent mother width, while conservatively tracking entropy and mother/daughter energy.
3. Include nonzero daughter mass \(M_2\) and coupling-dependent \(\eta=1-M_2^2/M_1^2\), thermal mother mass, and inverse-driven hot-tail cooling.
4. Improve full momentum-grid convergence of \(f_\phi(q)\) (not merely its integrated abundance) and physical IR regulation/condensate handling at immense Bose occupation.
5. Only **after** that use the resulting frozen \(f(q)\) in CLASS and compare to the externally computed halo-like dark matter distribution, including galactic thermalization times. The existing cold-parameter α scans and mixed-hot-tail CLASS scenarios are explicitly phenomenological comparison benchmarks, NOT derived kinetic outputs.

Related primary papers: Yin (2023) https://doi.org/10.1007/JHEP05(2023)180; Sakurai & Yin (2025) https://doi.org/10.1007/JHEP03(2025)202.


**Improvement (energy-conserving projection):** Instead of interpolation in log momentum, the discrete 1↔2 event is deposited barycentrically between the *exact momentum-shell mean energies* q̄ᵢ=(3/4)(q_hi⁴−q_lo⁴)/(q_hi³−q_lo³). This gives E_phi,event and E_chi2,event exactly after deposition, for every quadrature event inside the grid, while still matching parent energy E1 event by event. A numerical unit check enforces a maximum daughter projection error below 10⁻¹² in dimensionless energy. The gain/loss reaction is evaluated at all four receiving bin pairs and stiff BDF time stepping enforces positive occupations and Pauli ceiling. These are numerical **conservation** properties for the flat, fixed-mother-bath model, not the cosmological energy budget; the external mother χ1 reservoir is not evolved. The detailed f(p) remains grid-sensitive until an independent momentum-resolution study passes.


### Verified run: correct collision factor and conservation

The current [verified quantum-stimulated production benchmark](VERIFIED-QUANTUM-BURST.md) passed GitHub Actions #44892708035. Three grid resolutions conserve daughter number and enforce the fermionic Pauli ceiling, but the nonlinear operator does NOT exactly conserve physical parent-event energy, and yield a >760-fold number enhancement with over 99.8% of produced bosons having p/T<0.1 in this fixed-T fixed-source laboratory. Its integrated number shows ~1.06% grid spread; the mean boson momentum ~21.8% grid spread. Therefore the detailed spectrum **is not yet ready as CLASS input**, despite the robust qualitative burst. This separation of mathematical code validity and physical predictivity is intentional.
