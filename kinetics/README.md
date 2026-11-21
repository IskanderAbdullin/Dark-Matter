# First genuine Bose-stimulated gain/loss calculation for dark matter production

**SM BATH CONSERVATIVE RE-THERMALIZATION BUILDING BLOCK:** [SM-NUMBER-CONSERVING-BATH.md](SM-NUMBER-CONSERVING-BATH.md) implements a finite-q fermionic BGK elastic relaxation with chemical µ determined from actual conserved mother number, exact positive/Pauli-bounded time substep, decreasing fermion KL distance, and explicit opposite SM energy transfer. No hidden thermal FD reset; still phenomenological γ_el, SM bath temperature/gstar and 2→2 microphysics NOT simulated. [Literature requirements](LITERATURE-PHYSICAL-SCOPE-2026-11-21.md) distinguish thermal-plasma Yin (2023) from our previously confirmed collisionless-mother RD-BEST scenario.


**ACTUAL BEST RD-FRW 3-CUTOFF NO-BEC TEST:** [BEST-FRW-QMIN-SCAN.md](BEST-FRW-QMIN-SCAN.md) / [Actions #44924098843](https://github.com/IskanderAbdullin/Dark-Matter/actions/runs/44924098843) / [measurements](results/best-rd-frw-qmin-real-2026-11-21.json). At end tΓ=4e−5, N0 scales as qmin³ across qmin=1e−4,5e−5,2e−5 with N0/qmin³ scatter 0.00384%; total phi comoving number steady to 0.000133% across the cutoff sweep. This supports NO macroscopic early Bose-Einstein condensate in the particular RD model; further time/Nq/BEST integrator convergence and physical SM rethermalization remain necessary.


**SECOND-ORDER EXTERNAL BEST FRW IR VALIDATED:** [BEST-SECOND-ORDER-FRW-PILOT.md](BEST-SECOND-ORDER-FRW-PILOT.md) / [Actions #44924314432](https://github.com/IskanderAbdullin/Dark-Matter/actions/runs/44924314432). Native exprb2 (two quantum rate passes/time step) has fφ(qmin)=1.019622 at 32 steps versus separate physically redshifted early analytic 1.019684, only −0.00606% relative. 4→8→16→32 difference decreases ~5× per dt/2, unlike first-order exprb needing 128 steps with +0.49% error. [Permanent evidence](results/best-frw-exprb2-real-2026-11-21.json). The zero-mode n0 remains ~1% distinct from regular q<qmin analytic; full qmin/SM/frozen dark matter still unproven.


**Actual FRW 64→128-stepping and analytic redshift crosscheck:** [BEST-RADIATION-FRW-PILOT.md](BEST-RADIATION-FRW-PILOT.md), [RD-COLLISIONLESS-BOSE-REFERENCE.md](RD-COLLISIONLESS-BOSE-REFERENCE.md) and [permanent external data](results/best-rd-time-64-128-real-2026-11-21.json). fphi(qmin) 1.029882→1.024675, observed 0.508% temporal difference; 1st-order extrapolated 1.019468 vs independent analytic 1.019684 (0.021% apart). The exact early FRW kernel uses mother's initial FD frozen in qcom, not repeated thermal reservoir; no fully numerical or physical DM convergence yet.


**INDEPENDENT ACTUAL FRW SOFT SOURCE BENCHMARK:** [RD-COLLISIONLESS-BOSE-REFERENCE.md](RD-COLLISIONLESS-BOSE-REFERENCE.md) derives the exact massive two-body Bose gain A(q,t) for initially thermal massive mother **held fixed in comoving momentum**, not rethermalized during prescribed radiation expansion. Integrating time-dependent A yields fphi(qmin=1e−4,tΓ=4e−5)≈1.019684 and sub-qmin number≈1.709258e−14, close to asymptotic extrapolation of independent real BEST 4/8/16/32 data. The initial A(qmin) drops 27047→10929 due to redshift; fixed-a Bose exponential 1.9502 would be wrong in FRW. No closed/SM-realistic/frozen 86-eV/CLASS claim.


**Independent BEST ZERO MODE qmin³ NEGATIVE TEST:** [actual 3-cutoff experiment](BEST-QMIN-ZERO-MODE-SCAN.md), [pinned real data](results/best-qmin-real-2026-11-21.json), shows N0/qmin³ invariant to 0.019% for qmin=1e−4, 5e−5, 2e−5, whereas number in all φ remains ~6.6055e−7. Thus early zero-mode is smooth massive Bose IR phase volume, NOT a condensate. **BEST RD-FRW time-step survey** [4/8/16/32 steps](BEST-RADIATION-FRW-PILOT.md) resolves ~first-order φ (Δφ 0.259% for 16→32) versus only ~second-order pressure-work quadrature; [measurements](results/best-rd-time-refinement-2026-11-21.json). Next stage higher temporal resolution plus q grid and SM/zero-mode physical closure before 86eV PSD/CLASS.


**ACTUAL CLOSED FRW BEST MPI BASELINE:** [BEST-RADIATION-FRW-PILOT.md](BEST-RADIATION-FRW-PILOT.md) reports true third-party BESThep a=sqrt(t/t0), positive massive Bose + Pauli inverse collisions and dynamic mother, zero-mode; comoving charges conserved and physical d(a⁴ρ)/dln a=a⁴(ρ−3P) independently integrated. At tGamma=4e−5, a=1.1832: 4 vs8 steps phi 7.0678e−7 vs7.1429e−7 (**1.05% difference**); IR f and N0 **8.42% difference**, physical work residual decreases ~4×. [Permanent observed data](results/best-rd-frw-real-2026-11-21.json). Continue 16/32 steps; no SM gstarS / late-time 86-eV relic.


**PHYSICAL FRW EXPANSION BENCHMARK, NOT JUST H=0:** [BEST-RADIATION-FRW-PILOT.md](BEST-RADIATION-FRW-PILOT.md) sets BEST a(t)=sqrt(t/t0) for the entire dynamic χ1,χ2,massive φ/zero mode collision network. Independently computes P from q⁴/E physical quadratures, proper Ecom=a⁴ρ work law, 2 comoving charges, and all kinetic spectra. This is a fixed RD background without SM entropy exchange or backreaction and not yet a converged dark-matter freezeout.


**REAL CLOSED NONLINEAR TIME EVOLUTION:** [BEST-CLOSED-MASSIVE-DYNAMIC.md](BEST-CLOSED-MASSIVE-DYNAMIC.md) and [persisted run](results/best-closed-massive-dynamic-real-2026-11-21.json) now record actual pinned third-party BEST χ1/χ2/massive φ evolution at fixed a=1 on 4 versus 8 time steps and exprb_full; corrected finite initial daughter seed 1e−12 exp(−q). qmin occupancy fφ≈1.9506 and zero-mode N0≈3.238e−14; [analytical frozen-bath Bose exponential](BEST-EARLY-BOSE-EXPONENTIAL.md) predicts 1.95022 and N0≈3.2705e−14, with measured ~0.02% and ~1% differences, not a physical condensate or FRW dark-matter spectrum.


**DYNAMIC FULL 3-SPECIES BEST PILOT:** [BEST-CLOSED-MASSIVE-DYNAMIC.md](BEST-CLOSED-MASSIVE-DYNAMIC.md) is a pinned BEST MPI-protocol for actual evolution of massive phi + massless chi2 + massive chi1 with quantum inverse channels, separate zero mode n0, total physical energy and two conserved charges after each step, no external thermal bath or FRW. Includes fail-closed validation and explicit time-step/spectral convergence TODO; no macroscopic BEC inference from a single low-q occupation.


**REAL MASSIVE φ INDEPENDENT BEST CROSS-CHECK:** [BEST-MASSIVE-DAUGHTER-REFERENCE.md](BEST-MASSIVE-DAUGHTER-REFERENCE.md) verifies 2 actual phi masses × 2 q grids using pinned BEST v1.3.5 MPI in [Actions #44917857964](https://github.com/IskanderAbdullin/Dark-Matter/actions/runs/44917857964); 96-q physical φ number error 6–7e−6 relative, physical energy sum 5.6e−8 or 3.47e−9. [Measurements](results/best-massive-phi-real-2026-11-21.json). [IR physical cell audit](MASSIVE-BOSE-IR-CELL.md) shows huge initial q0 occupancy does not imply macroscopic condensed particle number; qmin occupation differs from q=0 by ~1.08% at mφ/M1=.01.


**OPTIONAL REAL BEST MASSIVE BOSON:** [BEST-MASSIVE-DAUGHTER-REFERENCE.md](BEST-MASSIVE-DAUGHTER-REFERENCE.md) specifies an independently pinned BEST source run with nonzero mφ/T and exact physical width-adjusted spin-summed amplitude. Checks mode-resolved φ production, true finite mφ deep-IR occupation and 3-species energy. The qmin species point is NOT a condensate n0, and zero-mode dynamics or physical 86-eV ratio is still unresolved.


**FINITE MASS BOSE IR SOURCE:** [MASSIVE-DAUGHTER-SOURCE.md](MASSIVE-DAUGHTER-SOURCE.md) derives the exact thermally averaged momentum-differential source for χ1→φ(mφ)+χ2(mχ), with continuous fixed-width number rate and finite φ occupation growth at q→0 whenever mφ>0. The strictly massless approximation has a nonuniform q0 limit: initial fφ(0) source can be enormous but the integrated low-q number is phase-volume suppressed; physical BEC is NOT proven.


**INDEPENDENT REAL FINITE-Bose/Pauli GAIN–LOSS:** [Actions #44916922599](https://github.com/IskanderAbdullin/Dark-Matter/actions/runs/44916922599) actually executes pinned BEST v1.3.5 with FD mother/FD chi2/BE phi at µ1=-.15,+.05 affinity, and µ1=-.25,-.05 affinity; both decay and inverse rate signs match exact 2D analytic source. At 96-q/24000 samples, best phi net 0.000856491254 vs exact 0.000856502725, but finite-grid physical energy −2.34e−7 and entropy relative bias −0.59% remain. [Full report](BEST-CHEMICAL-AFFINITY.md), [long-lived measured evidence](results/best-affinity-real-2026-11-21.json). No FRW or zero mode yet.


**AWAY-FROM-EQUILIBRIUM BEST PROTOCOL:** [best_chemical_affinity_reference.py](best_chemical_affinity_reference.py), [physics and instructions](BEST-CHEMICAL-AFFINITY.md) compare real pinned external 1↔2 Bose enhancement + fermionic Pauli + inverse decay against independently derived thermal-affinity exact R, 3 separate physical energy fluxes and entropy sign. Analytic reference verified by original 2D gain/loss quadrature; external non-equilibrium BEST run is not yet executed at this stage.


**EXTERNAL BEST 2×2 GRID RESULTS:** [genuine independent number/energy/thermal detailed-balance measurements](BEST-FIRST-REAL-RESULTS.md) in public [successful GitHub Actions run](https://github.com/IskanderAbdullin/Dark-Matter/actions/runs/44916394684). With grid n=56 and neval=12000, phi number −6.30e−6 relative and parent energy +6.26e−6 relative to exact source, but at n=28 energy still biases +0.004. Equilibrium µ1=µ2+µphi shows net/gross χ1=2.36e−5, φ=2.34e−6, χ2=3.40e−6. [Long-lived machine-readable evidence](results/best-independent-2026-11-21.json). No independent-seed, late-time FRW or zero-mode validation.


**ACTUAL EXTERNAL BEST SCIENTIFIC EXPERIMENT:** [GitHub Actions #44915936135](https://github.com/IskanderAbdullin/Dark-Matter/actions/runs/44915936135) succeeded running independent BEST v1.3.5 MPI source. First 28-node result reproduces Bose newborn number rate within 0.0291%, but physical mother energy flux differs +0.3964% and the 3-species physical energy rate sum is nonzero. [Run log, table, expanding to independent 2-grid × 2-neval and equilibrium checks](BEST-GITHUB-ACTIONS.md). Unlike earlier protocol-only tests, this is genuine external execution but still NOT numerical convergence or physical dark matter PSD.


**ANALYTIC THERMAL NONLINEAR CROSS-CHECK:** [THERMAL-CHEMICAL-AFFINITY.md](THERMAL-CHEMICAL-AFFINITY.md) derives a closed one-dimensional formula for finite-Bose/Pauli 1↔2 net source on exact two-daughter mass shell. Affinity Δµ=µ1−µchi2−µphi fixes gain/loss sign; qchi integral is analytic with a stable A=B limit. Direct independent quadratures and vacuum-limit tests validate it; not a coupled FRW or zero-mode computation.


**EXTERNAL BEST QUANTUM DETAILED BALANCE:** [BEST-THERMAL-BALANCE.md](BEST-THERMAL-BALANCE.md) and [optional runner](best_thermal_balance_reference.py) define equilibrium mother FD, daughter FD and Bose distributions with µphi=-.3, µchi2=.1, µ1=-.2. On-shell quantum gain/loss cancel exactly (also tested under chemical-potential violation), while real BEST mode-resolved gain/loss integration is pending independently installed MPI/Vegas.


**EXACT INITIAL DIFFERENTIAL 1↔2 SOURCE:** [BEST-VACUUM-SPECTRUM.md](BEST-VACUUM-SPECTRUM.md) and [best_vacuum_spectrum.py](best_vacuum_spectrum.py) analytically compute dR/dq=g1*MΓ/(2π²)*log(1+exp[-q−M²/(4q)]) with exact maximum q=M/(2T)=.05 at m=.1. Independent integration returns original Rvac=.0164959084596. The external BEST runner records a pointwise reference and robust q-shape diagnostic, not an observed cosmological freeze-out.


**REPRODUCIBLE PROGRESS MANIFEST:** [status document](RESEARCH-STATUS-2026-11-21.md) and [JSON audit trail](results/research-status-2026-11-21.json) enumerate this continuation's 14 Git commits, analytic vacuum number+energy benchmarks, fail-closed independent BEST protocol, 25.66% Lie x-step nonconvergence and missing full FRW energy closure. They explicitly prevent claiming verified physical DM or a CLASS-ready PSD.


**BEST ENERGY GATE:** [external validation protocol](BEST-SOURCE-COMPARISON.md) now fails closed unless the independently executed BEST output includes rates of physical energy of χ1/φ/χ2, agrees with ΓM1nFD and has zero total energy flux within tolerance. CI fixtures are artificial and external BEST has NOT yet run.


**BEST NUMBER+ENERGY CROSSCHECK:** [external BEST source protocol](BEST-INDEPENDENT-VALIDATION.md) now includes independent exact energy loss (M1 Gamma nFD), which at mass/T=.1 equals 0.0364681857069, in addition to 0.0164959084596 number rate; running BEST can reveal incorrect spin/dof factors or physical energy bookkeeping. No BEST external numerical result has yet been produced.


**BEST INDEPENDENT SOURCE VALIDATOR:** [BEST-SOURCE-COMPARISON.md](BEST-SOURCE-COMPARISON.md) documents pinned external Git SHA check before executing BEST and a [multi-run comparator](best_source_comparison.py) that fails closed on missing data, wrong spin normalization, incomplete 3-species rate stoichiometry, or insufficient q/neval variation. Synthetic test fixtures are not independent physics evidence; external BEST has not run.


**EXTERNAL INDEPENDENT VALIDATION PATH:** [BEST v1.3.5 protocol](BEST-INDEPENDENT-VALIDATION.md) inspects an MIT-licensed open-source independent collision+FRW Boltzmann solver with full quantum statistics and zero-mode examples. New [optional benchmark](best_external_vacuum_reference.py) prepares a mass=.1, gχ1=4, empty daughter 1↔2 rate test against exact analytic Rvac=.0164959084596. Ordinary scientific CI only tests PROTOCOL; no actual external Monte Carlo comparison has run.


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
