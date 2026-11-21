# First genuine Bose-stimulated gain/loss calculation for dark matter production

**QUANTUM REMAP ENTROPY:** [discrete Bose/Pauli entropy check](QUANTUM-REGRID-ENTROPY.md) verifies that positive q³ convex occupation averaging cannot decrease each species' quantum coarse-grained entropy under explicit vacuum UV extension; Pauli 0≤f≤1 and arbitrary Bose f≥0 preserved. This is numerical Jensen mixing, not a physical collision entropy theorem or validated f(q).


**ALL-THREE-SPECIES GRID EXTENSION:** [DAUGHTER-UV-REGRID.md](DAUGHTER-UV-REGRID.md) performs q³-accurate finite-volume transport of arbitrary positive Bose occupations (even f≫1) and Pauli-bounded fermions when daughter's upper energy domain grows to repair spurious reaction invariants. Tests exact individual numbers, both charges, parent phase volumes and a short closed new-grid collision continuation. UV-vacuum assumption and artificial energy/shape changes are reported explicitly.


**DYNAMIC UV-GRAPH REPAIR:** [automatic scan across instantaneous M1/T](DYNAMIC-RANK-UV-REPAIR.md) identifies how increasing daughter qmax reconnects the actual on-shell 1↔2 graph and removes spurious invariants beyond x=.025. A separate parent FD-number/energy/pressure check highlights that rank=3 alone is insufficient. This is a mesh feasibility tool, NOT evolved or CLASS-ready Bose dark-matter PSD.


**SEPARATE PHYSICAL WORK FROM NUMERICAL ERROR:** [FRW-ENERGY-WORK-AUDIT.md](FRW-ENERGY-WORK-AUDIT.md) validates free massive-mother x-derivative against independent q×x integrals and adds actual physical-energy endpoint differences/closure to both short FRW split schemes. Collision-shell and mother-remap defects are displayed rather than hidden as a falsely conserved FRW energy.


**TIME-INTEGRATION ALTERNATIVE:** [symmetric collision–Liouville–collision split](STRANG-FRW-COMPARISON.md) compares step-refinement of the early moving-on-shell 3-species FRW pilot against the original Lie operator split. It supports regridding when mother energy-mode counts change and reports both collision-half energy errors and Liouville/remap errors. This does not yet certify time convergence or full FRW physical energy closure.


**FRW TOPOLOGY CHANGES:** [cross changing mother-mode sets](FRW-MODE-CHANGES.md) by exact q³-overlap: shared parent energy states can change from 385 to 383 with two charges intact, while all new on-shell graphs are independently rank audited. At larger x the daughter grid develops fake invariants and is rejected. Explicit time-grid negative result: 2→4 split steps changes nphi by 25.66%, so no convergence claim.


**COUPLED EARLY-FRW PILOT (STRUCTURAL):** [NARROW-CLOSED-FRW-SPLIT.md](NARROW-CLOSED-FRW-SPLIT.md) combines instantaneously rebuilt on-shell collision channels, all three dynamic species and mother Liouville q³ remapping in a tiny RD x-window. It tests two full comoving number charges and separates collision/discrete energy, true free redshift pressure work, and remap-energy bias. It refuses parent-state changes rather than silently fabricating conserved physics.


**TWO-MOMENT FEASIBILITY:** [Pauli-limited shell energy bounds](PAULI-ENERGY-BOUNDS.md) compute exact finite-grid minimal and maximal energy for fixed mother number; an optional bounded transfer can enforce number+energy if feasible, but alters phase-space shape and must not be silently used in the actual FRW collision solver.


**CONSERVATIVE LIOUVILLE REMAP:** [parent_q_remap.py](parent_q_remap.py) moves freely streaming χ1 occupations between fixed-energy parent meshes at different x=M1/T by conservative q³ overlaps. Pauli and comoving number preserved; finite-volume energy drift and numerical mixing explicitly reported; no hidden claim of full coupled FRW conservation or CLASS readiness.


**TIME-DEPENDENT FRW ON-SHELL GRAPH REQUIRED:** [fixed-energy resonance drift no-go](FRW-RESONANCE-DRIFT.md) derives an exact nonzero physical energy defect when the same massive parent comoving momentum/old daughter-energy triplet is reused after x=M1/T changes. The original Schur-BDF flat-space on-shell network cannot be converted to FRW by merely changing Γ/H. Need dynamic parent/comoving remeshing with full pressure-work and quantum invariants.


**ENTROPY-AWARE LIOUVILLE:** [variable gstarS collisionless χ1 benchmark](FRW-ENTROPY-RELEASE.md) derives d ln a/d ln x=1−(1/3)d ln gS/d ln x and x d(ρ/T⁴)/dx=(1+α)(ρ/T⁴)+(α−3)(P/T⁴). It verifies true comoving n a³ constancy while n/T³ varies with SM entropy dilution. gstarS step is a test function, not the actual SM thermal history or closed 1↔2 model.


**FINITE-VOLUME FRW MOTHER:** [constant-q-bin collisionless Liouville scheme](FRW-MOTHER-FINITE-VOLUME.md) preserves initial FD number at ~1e−12 and exactly obeys x d(ρ/T⁴)/dx=(ρ−3P)/T⁴ for every grid, but independently benchmarked pressure still errs +1.33% at 100 bins (x=5). Conservation is insufficient to certify phase-space resolution. Γ=0 only; no 1↔2 collisions yet.


**COLLISIONLESS MASSIVE FRW LIMIT:** [exact Liouville parent background](FRW-COLLISIONLESS-MOTHER.md) distinguishes free f1(q)=constant from external continuous FD rethermalization as M1/T rises; independent number/energy/pressure quadratures obey x d(ρ/T⁴)/dx=(ρ−3P)/T⁴ and dρ/dln a=−3(ρ+P). Required independent benchmark before coupling the closed 1↔2 collision solver to expanding FRW, not a dark matter PSD.


**POSITIVE RESONANT RATE CALIBRATION:** [cut-cell FD averaged K_ij](CALIBRATED-CUT-CELL-RATES.md) retains exactly the same physical center-on-shell 1↔2 channels and unique mother states. Reference-thermal empty-daughter source matches the included continuum cut-cell integral; forward/reverse rates share one positive constant K, so conservation/H theorem remain. Dynamic nonequilibrium and zero-mode/freezeout accuracy remain unproven.


**CUT-CELL DIAGNOSTIC:** [exact threshold geometry + Fermi-Dirac source averages](CUT-CELL-SOURCE-DIAGNOSTIC.md) independently integrate the original on-shell domain over every cell (with 16→32 Gauss convergence), separating finite-cell FD center approximation from physical allowed regions incorrectly masked because the *central* reaction is off-shell. This is not a shortcut to including those off-shell center reactions in a conserved discrete kinetic network.


**UV SOURCE ERROR DECOMPOSED:** [exact rectangular continuum source versus finite-cell bias](UV-SOURCE-DECOMPOSITION.md) analytically integrates the on-shell 1↔2 FD kernel with both daughter energies truncated at grid upper EDGE Q. The finite-cell error plateaus near −0.0103% across Q~5–10T while missing physical high-q phase space separately falls exponentially. Independent double integration tests validate the analytic expression. No nonlinear frozen f(q) claim.


**SOURCE QUADRATURE AUDIT:** [vacuum source and nonlinear cutoff sensitivity](VACUUM-SOURCE-CONVERGENCE.md) compares exact continuum on-shell Rvac, finite-bin ΣK_ij fFD(ei+ej), and the actual closed nonlinear ODE initial derivative. It shows that even a near-correct initial vacuum rate does not ensure nonlinear production convergence; more bins can make center/cut-cell source errors worse.


**PHYSICAL ON-SHELL IR THEOREM:** [vacuum small-q source bound](VACUUM-IR-TAIL-BOUND.md) integrates the exact Fermi–Dirac on-shell daughter-energy source and proves a sharp exponential-integral upper bound for empty daughters. For M1/T=.1, initial birth below q/T=.0005 is 6.1e−7 of total, but this does NOT constrain Bose-stimulated later f(q). Independent numeric tests enforce the bound.


**SPECTRAL SHAPE DIAGNOSTICS:** [finite-bin CDF/number quantile audit](SHAPE-DIAGNOSTICS.md) reconstructs number-weighted momentum q01/10/50/90/99 from the unchanged Schur-BDF short-time daughter spectra. It checks exact q²dq shell-volume CDF and separates bin-dependent Bose heights from integrated observables. No physical CLASS prediction; must still verify IR/UV/time/FRW and condensate dynamics.


**QUANTIFIED SHORT-TIME CUTOFF ERROR:** [full dynamic 1↔2 cutoff/time audit](CUTOFF-AND-TIME-ERROR.md) and [reproducible script](cutoff_and_time_convergence.py) compare independent parent FD-number/energy/pressure error, BDF time error and actual closed boson source for Rmax 5/7/8.5/10T. Near Rmax=8.5→10 the source number changes by ~0.44%, while the initial mother number error is tiny; full shape/freeze-out still not demonstrated.


**SCHUR-BDF PROGRESS:** [exact daughter-block Schur complement](SCHUR-BDF-RESULTS.md) accelerates unchanged closed nonlinear quantum kinetics: 299 daughters + 10,199 unique mother states. A four-grid full dynamical study (Rmax 5000/7000/8500/10000) finds significant cutoff dependence of tGamma=.001 Bose abundance despite independently well-converged mother FD moments. No physical primordial spectrum/CLASS approved.


**SPARSE KINETIC DYNAMICS:** [analytic sparse Jacobian + closed three-species evolution](SPARSE-CLOSED-KINETICS.md) verifies the same nonlinear RHS against the original adaptive operator and independently checks Jacobian finite differences. A 5,497-state numerical test with 33,191 active 1↔2 pairs conserves physical charges/energy and dynamically depletes χ1. This is short fixed-T kinetics, NOT a cosmological freezeout or CLASS spectrum.


**NEW INVARIANT-CLEAN GRID:** [resonant graph repair / initial parent FD moments](RESONANT-GRAPH-REPAIR.md) solves the *specific* spurious-invariant problem for tested bridged fine-core/coarse-comb meshes. Three physical graph invariants (not 4/14), no isolated daughters, and parent number/energy/pressure errors below 4e-5 by 269 daughter modes. Runtime cost is thousands of mother-energy states; this is NOT an evolved PSD or a physical DM spectrum.


**Closed-form spurious conserved charge:** the [three-band modular-invariant witness](modular_invariant_witness.py) proves exactly that a residue-weighted pseudo-energy is conserved by every on-shell reaction and by the full nonlinear 3-species ODE. This falsifies the apparent uniqueness of the kinetic equilibrium even though energy, two number charges and H theorem hold. See [detailed derivation](ADAPTIVE-INVARIANT-AUDIT.md).


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
