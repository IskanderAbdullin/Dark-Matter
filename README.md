# Bosonic dark-matter halos with a dominant noncondensed component — research program

**Independent BEST MASSIVE 86-eV-like scalar validation prepared:** [pinned BEST conditional mφ>0 on-shell 1→2 benchmark](kinetics/BEST-MASSIVE-DAUGHTER-REFERENCE.md) uses physical two-body matrix element M²=8π g1 M1²Γ/p* for fixed rest width, rather than invalid massless amplitude. Measures finite-mass differential φ source and analytically defined q0 occupation seed against independent formula, with three species number and physical energy. This reference is protocol-only until actual independent BEST run; no condensate/FRW/CLASS claim.


**Finite 86-eV boson mass and deep-IR nonuniformity:** [exact massive two-body χ1→φ+χ2 vacuum source](kinetics/MASSIVE-DAUGHTER-SOURCE.md) gives physical two-body daughter energy spectrum and its nonzero finite-mφ q→0 Bose occupation source, contrasting the strictly massless q=0 limit. It analytically preserves number and physical energy with fixed rest width for arbitrary daughter masses below threshold. This is NOT condensate number, a physical chosen mφ/Mχ1 model or cosmological freezeout.


**Genuine independent NON-EQUILIBRIUM Bose/Pauli 1↔2 verified:** [BEST GitHub Actions #44916922599](https://github.com/IskanderAbdullin/Dark-Matter/actions/runs/44916922599) directly ran externally pinned BEST at finite BE+FD occupation with nonzero chemical affinity and inverse decay. [Actual observed vs independent exact 2D rate+energy and entropy numbers](kinetics/BEST-CHEMICAL-AFFINITY.md), [permanent 3-case numerical evidence](kinetics/results/best-affinity-real-2026-11-21.json). At 96 q nodes/neval=24000, φ number rate differs −0.001339% and differential shape RMS/peak 9.60e−6; physical total energy residual −2.34e−7 and entropy still errs ~0.59%. Reverse-chemical driving changes signs as predicted. This validates a snapshot quantum source, NOT nonlinear cosmological time evolution, BEC or CLASS.


**FIRST independent non-equilibrium quantum source protocol:** [BEST finite-Bose, Pauli and reverse-decay check](kinetics/BEST-CHEMICAL-AFFINITY.md) can run real pinned BEST at common T but nonzero chemical affinity Δµ=.05 and compare all three number+energy rates and differential phi source against [a separate exact analytic/full-2D thermal quantum collision benchmark](kinetics/THERMAL-CHEMICAL-AFFINITY.md). Targets nphi dot=0.00085650272535, physical energies χ1 −0.00181019828605, φ +0.00067304518789, χ2 +0.00113715309817. Tests at this stage are protocol-only, not an executed non-equilibrium BEST run or a FRW relic.


**ACTUAL independent 2×2 BEST physics result:** [real MPI/Vegas 28-versus-56 q-grid survey, energy and inverse-decay detailed balance](kinetics/BEST-FIRST-REAL-RESULTS.md), [permanent numerical evidence JSON](kinetics/results/best-independent-2026-11-21.json), [successful Actions run #44916394684](https://github.com/IskanderAbdullin/Dark-Matter/actions/runs/44916394684). On 56 q nodes, neval=12000, independent φ rate differs −0.0006303% and parent energy loss +0.0006258% from our analytic 1↔2 source; on 28 q nodes energy error is still ~+0.40%. Nonempty Bose+Pauli equilibrium sources cancel to net/gross ~2e−5 or better. This validates initial reference only, NOT a nonlinear frozen DM or CLASS result.


**FIRST ACTUAL INDEPENDENT BEST RUN PASSED:** [public GitHub Actions #44915936135](https://github.com/IskanderAbdullin/Dark-Matter/actions/runs/44915936135) executed pinned external BEST v1.3.5 with two MPI ranks, separate from our model. At M/T=.1, 28 q nodes and neval=3000, φ number rate 0.0164911055 vs exact 0.0164959085 (−0.0291%); mother energy loss 0.0366127458 vs exact 0.0364681857 (+0.3964%) with nonzero total energy imbalance. [Exact raw run and next 2×2 independent convergence workflow](kinetics/BEST-GITHUB-ACTIONS.md). This single coarse run is NOT a validated nonlinear or cosmological freeze-out; expanded 2×2 plus inverse-decay BEST checks are queued in code.


**Exact finite-occupation 1↔2 source derived:** [thermal Bose/Pauli chemical-affinity benchmark](kinetics/THERMAL-CHEMICAL-AFFINITY.md) analytically integrates the full nonlinear forward-minus-inverse kernel over one on-shell daughter momentum at finite BE and FD occupations, proving its sign is sign(mu1−mu2−mu_phi) for common T. Independent direct double quadrature gives Rphi=-0.0008263650 / 0 / +0.0008565027 for mu1=-.25 / -.2 / -.15; its dilute daughter limit reproduces Rvac=.01649590846. This is instantaneous kinetics only, not frozen dark matter.


**Independent quantum reverse-channel benchmark:** [exact thermal Bose–Fermi detailed balance with nonzero daughter occupations](kinetics/BEST-THERMAL-BALANCE.md) now tests the full on-shell A=f1(1−f2)(1+b) vs B=(1−f1)f2b identity at µ1=µphi+µchi2 and T=1, and deliberately breaks chemical equilibrium to verify gain/loss signs. An optional pinned BEST v1.3.5 runner tests mode-resolved physical gain/loss, all three number and energy rates. External BEST has NOT yet been executed; this is an analytical and CI protocol.


**Independent pointwise Bose source benchmark:** [exact initial dRvac/dqφ and peak at q=M1/(2T)](kinetics/BEST-VACUUM-SPECTRUM.md) extends the external BEST comparison beyond integrated rates. At M1/T=.1 the empty-daughter thermal source peaks at q=.05 and is exponentially suppressed at deep q, while the later stimulated Bose population may be much colder. The optional BEST runner now returns its grid-level differential source and robust relative shape errors; external BEST has still NOT run.


**Research snapshot / audit trail (2026-11-21):** [14 consecutive commits, negative FRW step-convergence finding, known kinetic/energy gaps, precise next tasks and independent BEST status](kinetics/RESEARCH-STATUS-2026-11-21.md); [machine-readable evidence manifest](kinetics/results/research-status-2026-11-21.json). New independent BEST v1.3.5 is identified, but has NOT been executed; the full 86-eV physical DM model, CLASS PSD and observational likelihood remain **unvalidated**.


**Independent BEST checks strengthened:** [multi-run validator](kinetics/BEST-SOURCE-COMPARISON.md) now also requires all three energy source derivatives and their total conservation, plus analytic mother energy-loss magnitude ΓM1nFD=0.0364681857069. Missing/wrong energy is a failed preliminary comparison even when all three particle-number rates match. Still no external BEST run; this tests protocol/math only.


**Independent on-shell energy target:** the optional [BEST external rate benchmark](kinetics/BEST-INDEPENDENT-VALIDATION.md) now compares not only newborn daughter number Rvac=0.0164959084596 but also the exact mother energy loss rate ΓM1 nthermal=0.0364681857069 at M/T=.1,g1=4,T=1. All three species' BEST energy moments and their sum are part of the external JSON, pending a real run.


**External BEST comparison now fail-closed:** [pinned SHA verification and multi-run source consistency checker](kinetics/BEST-SOURCE-COMPARISON.md) accept only supplied outputs created under matching BEST v1.3.5 Git SHA and compare thermal vacuum 1↔2 rates for mother and both daughters, with ≥2 grids/≥2 Vegas neval settings for a preliminary agreement label. No actual BEST Monte Carlo results have yet been supplied; current external status remains UNVERIFIED, and no nonlinear FRW or CLASS validation is claimed.


**Independent solver identified (not yet executed):** [BEST v1.3.5 momentum-resolved quantum Boltzmann framework](kinetics/BEST-INDEPENDENT-VALIDATION.md), source best-hep/best pinned c70875be. It supports on-shell 1↔2, evolving mother/daughters, Bose/Pauli, FRW, baths and zero-mode examples. A safe optional [standalone vacuum-source benchmark protocol](kinetics/best_external_vacuum_reference.py) specifies the independent analytic rate 0.0164959084596 at M/T=.1 and required stoichiometric tests. External BEST has NOT been run or verified; no claimed independent agreement.


**Quantum entropy audit on adaptive UV remaps:** [Bose/Fermi concavity under q³ averaging](kinetics/QUANTUM-REGRID-ENTROPY.md) shows nonnegative coarse-grained entropy change for all three species even with huge Bose f, alongside exact individual numbers and charge conservation. This is Jensen numerical mixing, NOT physical entropy generation or proof of physical FRW energy conservation; high-q domain still assumes unrepresented vacuum.


**Three-species UV regrid + collision continuation:** [Bose/Pauli-safe daughter and mother momentum-space remapping](kinetics/DAUGHTER-UV-REGRID.md) now crosses x=.02→.03 with daughter Rmax 300→500 to keep an invariant-clean reaction graph. It exactly preserves all three species counts and two independent charges under q³ overlap while explicitly assuming unrepresented high-q cells are vacuum. Parent phase-volume normalization is synced to the new closed ODE, and a short new-grid collision continuation is tested. Physical momentum shape/energy remains unvalidated; not a CLASS relic.


**Instantaneous FRW graph connectivity repair:** [dynamic upper-cutoff rank scan](kinetics/DYNAMIC-RANK-UV-REPAIR.md) finds how expanding the daughter momentum domain prevents fake conserved collision charges as M1/T rises: on the tested fine-core/comb meshes, x=.025/.03/.04/.05/.08/.1 first pass three-invariant gate at Rmax 400/500/1000/1500/3500/5000 respectively. The minimum graph-valid mesh is NOT necessarily thermal-moment or nonlinear-source converged; physical daughter remapping and full FRW energy closure remain unsolved.


**FRW physical-energy audit:** [independent massive-mother pressure-work integration](kinetics/FRW-ENERGY-WORK-AUDIT.md) verifies the exact collisionless identity d(rho/T^4)/dx=(rho-3P)/(x T^4) via separate q/x quadratures. Both left-endpoint and symmetric split-step solvers now explicitly close a telescoping physical-energy ledger: legitimate pressure work, physical collision-shell bias and numerical mother-remap drift remain separate. Zero ledger residual does NOT imply the full physical energy equation is solved or a CLASS-ready relic exists.


**Symmetric early-FRW split diagnostic:** [Strang-type collision–Liouville–collision](kinetics/STRANG-FRW-COMPARISON.md) tests the same closed dynamically rebuilt 1↔2 networks as the left-endpoint Lie scheme and measures, rather than assumes, step refinement. Both half-collision discrete-energy errors, mother physical shell-energy defects, free FRW work, projection energy drift and two global comoving charges are tracked separately. Threshold nonsmoothness prevents an automatic second-order claim; no CLASS-ready relic is inferred.


**FRW parent-mode threshold crossing:** [dynamic q-conservative regridding across a changing parent-energy basis](kinetics/FRW-MODE-CHANGES.md) now allows mother resonance modes to disappear as GRID VARIABLES (385→383) without deleting physical mother particles, provided each freshly rebuilt on-shell graph still has precisely three invariants. The numerical x-stepping is NOT converged: 2 vs 4 split steps at κ=500 changes produced Bose number by ~25.66%. No full energy-conserving FRW relic or CLASS PSD is claimed.


**First explicitly coupled early-FRW structural pilot:** [narrow closed collision+Liouville splitting](kinetics/NARROW-CLOSED-FRW-SPLIT.md) now rebuilds on-shell χ1↔χ2+φ reaction graphs as x=M1/T changes across .0200→.0202, evolves all species with sparse BDF, and conservatively remaps the massive mother on fixed q=p/T domains. Both independent comoving charges and collision-step discrete energy are checked, while physical shell-energy and remap errors are explicitly reported. This is NOT a self-consistent FRW relic calculation or validated spectrum.


**Two-moment Pauli feasibility:** [exact attainable massive-parent energy interval at fixed number](kinetics/PAULI-ENERGY-BOUNDS.md) derives E_min/E_max by shell filling under fermion occupancy bounds. An optional artificial N+E correcting remap is diagnostic only: a changed momentum distribution can preserve both integrated moments but is NOT physical collisionless Liouville or guaranteed entropy/detailed balance.


**Liouville momentum remap:** [number-exact Pauli-safe O(N) parent regridding](kinetics/PARENT-Q-REMAP.md) now transfers massive χ1 occupations between moving energy-resonant meshes using q³ cell overlaps. The scheme preserves comoving particle number and 0≤f≤1, but introduces a separately measured numerical energy/entropy error. This is a collisionless component, not a coupled energy-conserving FRW dark-matter solver.


**A fundamental coupled-FRW obstacle is now proven:** [fixed resonance energies cannot stay on shell as a massive mother redshifts](kinetics/FRW-RESONANCE-DRIFT.md). For fixed daughter comoving energies e_i+e_j and the same freely streaming parent momentum, the latter's energy changes as sqrt(E_parent(x0)^2−x0²+x²), creating nonzero energy mismatch for every original channel whenever x!=x0. A future closed FRW 1↔2 solver must dynamically construct/reassign physically resonant channels and preserve full energy, positivity and charge. This is an exact kinematic obstruction, not a solved cosmological kinetic model.


**Entropy-release FRW invariant solved:** [collisionless massive-parent background with variable SM gstarS](kinetics/FRW-ENTROPY-RELEASE.md) proves fixed q=p/T is invalid during entropy dumps; the conserved free variable is q_c=p/[T(gstarS/gstarS_dec)^(1/3)]. An independent synthetic 55→3.909 entropy-step test conserves n*a³ and exactly verifies energy-pressure continuity. This is no-collision physics, not a complete Bose production model or observed dark-matter relic.


**Finite-volume FRW Liouville milestone:** [collisionless massive χ1 fixed-comoving-q bin model](kinetics/FRW-MOTHER-FINITE-VOLUME.md) now phase-volume-averages the initial FD state, preserving comoving number and physical FRW pressure-work continuity to roundoff. But at x=M/T=5, just 100 bins give 0.275% energy and 1.33% pressure error even though the initial number is accurate to ~1e-12. This is a strict Γ=0 prerequisite, NOT the coupled Bose-production or CLASS spectrum.


**FRW Liouville benchmark for MASSIVE dynamically conserved mother:** [collisionless FRW χ1 energy/number/pressure check](kinetics/FRW-COLLISIONLESS-MOTHER.md) proves that after decoupling fχ1(q=p/T) stays fixed while E/T=sqrt(q²+x²) changes. Comoving parent n/T³=.3652061 remains constant, whereas artificially reapplying instantaneous FD at x=5 leaves only 7.35% of that inventory. Exact dρ/dln a=−3(ρ+P) is tested. This is a no-collision benchmark, NOT coupled closed Bose FRW production.


**New positive quadrature variant:** [thermal cut-cell calibration of on-shell quantum reaction rates](kinetics/CALIBRATED-CUT-CELL-RATES.md) integrates the exact empty-daughter FD source within retained center-resonant cells and sets a positive constant rate that preserves the graph's two charges, energy, detailed balance, Pauli bounds and H theorem. It does NOT recover omitted off-shell-centroid regions or establish nonlinear/FRW primordial PSD accuracy.


**On-shell cut-cell error isolated:** [threshold-cut source diagnostic](kinetics/CUT-CELL-SOURCE-DIAGNOSTIC.md) now integrates the Fermi-Dirac gain exactly in one daughter energy inside every geometrically admissible cell and compares summed high-order Gauss quadrature against the independent rectangular continuum source. It separates centroid FD sampling bias from physically allowed cells discarded because their centroid is off-shell. This is a SOURCE-ONLY test, not an alteration or validation of the nonlinear Bose kinetic operator.


**Physical source-error decomposition solved analytically:** [finite daughter-domain cutoff vs cell discretization](kinetics/UV-SOURCE-DECOMPOSITION.md) uses the exact continuum qphi,qchi<Q on-shell FD integral. On the repaired 3-invariant mesh, Rmax=8500 incurs −0.0487% true phase-space truncation and −0.0103% finite-cell error; at Rmax=10000 both are around −0.011%. These are EMPTY-daughter source diagnostics, not bounds on the nonlinear Bose occupation, dark-matter PSD or cosmological relic.


**Initial source vs nonlinear rate calibration:** [independent 1D vs finite-cell on-shell source audit](kinetics/VACUUM-SOURCE-CONVERGENCE.md) now verifies the complete CLOSED ODE's t=0 φ derivative equals its summed thermal reaction rates and compares both to exact continuum 0.0164959084596. Changing UV cutoff 8.5→10T alters this initial vacuum rate ~0.038% but the nonlinear tΓ=.001 φ number ~0.44%: this is a numerical sensitivity diagnostic, NOT the Bose enhancement factor.


**Analytic infrared control added:** [rigorous bound for the EMPTY-daughter on-shell vacuum source](kinetics/VACUUM-IR-TAIL-BOUND.md) yields a fractional birth rate below q/T=.0005 of ~6.1e-7 (bound 4.1e-6) for M1/T=.1. It is explicitly NOT a bound on the later nonlinear stimulated Bose population or zero-mode condensate; the factor 1+f can radically change production.


**Beyond a single mean momentum:** [closed Bose-bin CDF and q01/q10/q50/q90/q99 diagnostics](kinetics/SHAPE-DIAGNOSTICS.md) now use the exact q²dq number measure of each kinetic shell, with independent analytic tests on constant-occupation spectra and normalized CDF comparisons across high-energy cutoffs. Snapshot analysis remains at short tΓ=.001, NOT cosmological freezeout, zero-mode BEC or CLASS input.


**Independent dynamic-convergence audit:** [source evolution versus static FD moment errors](kinetics/CUTOFF-AND-TIME-ERROR.md) shows the closed Bose number at tΓ=.001 shifts ~0.44% between high-energy cutoffs 8.5 and 10T even though initial mother FD number is already accurate to ~5.6e−6 fraction. Tightening BDF rtol reduces number changes to ~0.014%, so cutoff error currently dominates numerical time error. This is NOT a converged primordial PSD or CLASS result.


**Major numerical advance:** [exact Schur-complement implicit BDF](kinetics/SCHUR-BDF-RESULTS.md) now evolves the invariant-clean closed 1↔2 network with >10,000 mother energy modes by eliminating their diagonal Jacobian block *only inside Newton linear solves*, not from physics. The short-time M1/T=.1 production number changes by 0.44% from upper mesh cutoff 8.5 to 10T despite initial mother thermal-density error below 1e-5, so the frozen spectrum remains unvalidated.


**Sparse closed evolution now runs on the repaired 3-invariant graph:** [sparse analytic-Jacobian solver](kinetics/SPARSE-CLOSED-KINETICS.md) evolves 199 daughter + 5,099 shared mother-energy states at M1/T=0.1. A short tΓ=.001 BDF run yields nφ/T³≈0.00859347, fully conservative at ~4.5e-16 for charge/energy, without a fake reaction-graph invariant. The initial mother energy grid is still biased (about −1.07%), and no FRW freezeout or CLASS PSD is claimed.


**Adaptive graph obstacle partially solved:** [resonant graph repair](kinetics/RESONANT-GRAPH-REPAIR.md) bridges the infrared core (100 consecutive nodes) to a coarse energy comb. Tested grids with 199–299 daughter modes have exactly **three** active collision invariants, zero isolated daughters, and independently converging initial thermal mother number/energy/pressure (at 269 nodes: errors 5.6e-6, 3.8e-5, 3.8e-5). Important: these networks require 5–10k shared parent-energy states. No dynamic freeze-out or CLASS validation is claimed.


**Constructive proof of grid aliasing:** for the three-band adaptive mesh, a modulo-2500 weighting of daughter and mother energy-cell indices yields an additional *exactly conserved* pseudo-charge unrelated to physical number/energy. It is verified on the complete nonlinear collision ODE, not merely found by SVD. See [analytic formula](kinetics/ADAPTIVE-INVARIANT-AUDIT.md#8-явная-аналитическая-формула-для-лишнего-инварианта-трёхдиапазонной-сетки) and [regression program](kinetics/modular_invariant_witness.py). This mesh is rejected for any physical CLASS initial spectrum.


**Reproducible negative-result artifact:** [JSON spurious-invariant benchmark](kinetics/results/spurious-invariants-2026-11-21.json) and [fail-closed CLASS readiness gate](kinetics/cosmology_readiness_gate.py) now block both unvalidated adaptive spectra, even when the nonlinear ODE conserves energy. All missing physical/IR/FRW evidence defaults to a failed release criterion.


**Adaptive-kinetics warning (2026-11-21):** A graph-rank audit of exact-resonance 1↔2 grids found additional numerical conservation laws beyond the two physical charges and energy. A sparse 27-bin cold/hot grid has 11 extra invariants, and a 56-bin three-band grid has one, despite conservation of the intended moments and the discrete H-theorem. A concrete extra invariant is conserved by the full nonlinear adaptive ODE. These grids are explicitly rejected as physical CLASS inputs. See [rank/ODE audit](kinetics/ADAPTIVE-INVARIANT-AUDIT.md), [code](kinetics/collision_invariant_rank.py) and [adaptive proof-of-concept](kinetics/adaptive_resonant_lattice.py).


**Conditional 86-eV relic-density link:** [freeze-out number to Ωh² normalization](cosmology/FRW-RELIC-ABUNDANCE-BRIDGE.md) uses the existing checked thermal scalar benchmark and standard assumed entropy dilution, finding required nφ/T³≈0.12268 for Ωh²=.12 at mass 86 eV, g*Sprod=55. A closed initial-mother population would have to convert ≥33.6%, so source depletion matters even near the target. This is conditional amplitude normalization, not a kinetic, structure, or observational validation.


**New mother-budget falsifier:** [closed-system particle-number audit of the fixed-bath FRW pilot](kinetics/FRW-RESERVOIR-BUDGET.md). A closed χ1↔χ2+φ model beginning at M1/T=.05 with empty daughters has initial nχ1/T³≈0.365206. The high-coupling fixed-bath κ=.2 run produces nφ/T³≈0.36948, **exceeding the initial mother inventory**; it therefore necessarily relies on external replenishment. This is not an inconsistency of the explicitly open pilot, but prohibits presenting it as a closed relic-production calculation.


**Closed three-species benchmark now implemented:** [resonant-energy on-shell decay lattice](kinetics/CLOSED-THREE-SPECIES.md) evolves the mother and both daughters without an external thermal reservoir. Its full nonlinear ODE conserves two independent number charges and total discrete energy, obeys Pauli/Bose bounds and a quantum H theorem. The deliberately coarse uniform lattice is a structural test, NOT an IR-resolved cosmological spectrum. [Next closed FRW specification](kinetics/CLOSED-FRW-NEXT-SPEC.md).


**FRW next-step pilot:** [radiation-dominated 1↔2 kinetic evolution in comoving q=p/T](kinetics/FRW-FIXED-BATH-PILOT.md) explicitly shows that a cold burst at fixed T does not automatically imply a cold frozen primordial population. The mother is still externally thermal, the coupling κ=Γ/H(T=M1) is illustrative, and the energy/entropy budget is not closed.


**New on-shell collision construction:** the 1↔2 measure transforms exactly to a two-dimensional daughter-energy domain (4q_phi q_{chi_2}ge M_1^2). An analytic 1D vacuum-source integral matches the original parent-angle quadrature at relative 1.7e−10; a finite-grid physically on-shell collision network now has exact daughter invariant, discrete balance and entropy tests. [Full derivation and numerical audit](kinetics/ONSHELL-COLLISION-2D.md). Still a fixed thermal mother bath; not FRW or validated 86-eV DM.


**New numerical finding:** finite Gauss–Legendre angle sampling misses the strongly backward-collinear soft-boson endpoint in the 1↔2 source. A verified change of variables to log daughter energy dramatically reduces angle-resolution bias in the *finite-lattice surrogate*, but does not remove its independent physical parent-energy defect. See [angular endpoint audit](kinetics/ANGULAR-ENDPOINT-IR-AUDIT.md) and [new log-q quadrature](kinetics/LOGQ-ANGULAR-QUADRATURE.md).


**Research status update (2026-11-21):** A separate audit found that the historic 1<->2 kinetic benchmark's event-level detailed balance and geometric energy remap do **not** imply discretized equilibrium or physical nonlinear energy conservation. See [audits](kinetics/EQUILIBRIUM-AND-IR-AUDIT.md), [energy defect](kinetics/ENERGY-PROJECTION-DEFECT.md), and [an entropy-balanced lattice surrogate](kinetics/ENTROPY-BALANCED-LATTICE.md). The latter has exact lattice equilibrium and H-theorem but changes mother energy in the finite-grid statistical factor; it is NOT a validated cosmological production model.

**Date:** 2026-11-21 | **Status:** substantive preliminary research, five preprint-style research drafts; not submitted or peer reviewed.

## Provenance
Primary source: LaTeX archive \`rv2_JCAP_041P_0721(1).zip\` uploaded by the user, containing \`rv2_JCAP_041P_0721.tex\`, \`lit2.bib\`, and 28 EPS figures. Publication: I. G. Abdullin and V. A. Popov, [*Boson dark matter halos with a dominant noncondensed component*](https://doi.org/10.1088/1475-7516/2021/11/055), JCAP 11 (2021) 055, [arXiv:2107.06237](https://arxiv.org/abs/2107.06237).

Source claims: noncondensed Bose component with equilibrium Thomas–Fermi density; condensed particles as compact Bose stars; 22 dark-matter-dominated LSB rotation curves; 20 retained in a joint mass scan yield preferred \`m=86 eV\`; two outliers modeled with additional Bose-star component. These are *reported* in the source, not independently reproduced observational results.

## Files
- \`00-critique-and-roadmap.md\`: independent scientific audit, hazards, priorities, milestones.
- \`literature/100-screened.md\`: exact 100-item screening inventory, and review-level honesty.
- \`papers/01-critical-core-asymptotics.md\`: **analytical result**: nonanalytic linear central density at critical degeneracy.
- \`papers/02-finite-halo-truncation.md\`: phase-space truncation and boundary-value formulation to avoid infinite-mass isothermal tails.
- \`papers/03-cosmological-transfer.md\`: thermal/nonthermal relic audit and small-scale-structure falsifiability.
- \`papers/04-hierarchical-galaxy-test.md\`: rigorous baryonic and global-mass inference against SPARC.
- \`papers/05-condensate-population.md\`: Bose-star kinetic/population/lensing predictions and survival criteria.
- \`research-log.md\`: observed facts, transformations, calculations, unverified requirements.

## Ground rules
No fabricated catalog fits, posterior constraints, N-body simulations, particle-production model or observational detections. Paper 1 includes a self-contained analytical derivation; the other four are theory/protocol manuscripts with *testable* claims and clearly listed outstanding computations. The original particle-mass estimate is **not** presented as a cosmologically established mass measurement. All new equations are for clearly stated approximations.


## Public migration and automated scientific checks

The 141 research files in this public repository were transferred **byte-for-byte** from the dedicated research folder of the private Info repository at commit d3e341adfce7b6dc012a4b193e7b9d55a8369a2d. The original private repo and its other projects were left unchanged. Public-only workflow, provenance, and tests are maintained here.

[![Scientific verification](https://github.com/IskanderAbdullin/Dark-Matter/actions/workflows/scientific-ci.yml/badge.svg?branch=main)](https://github.com/IskanderAbdullin/Dark-Matter/actions/workflows/scientific-ci.yml)

See [automatic test coverage and limits](TESTING.md), [continuation V](continuation-05/README.md), [theoretical continuation IV](continuation-04/CONSOLIDATED-RESULTS-RU.md), and [SPARC/kinetic continuation II](continuation-02/CONSOLIDATED-REPORT-RU.md).

**Research status:** exploratory mathematical and numerical work. A green CI run verifies specified reproducibility checks, not particle discovery, galaxy formation or cosmological viability.


## New energy-space research and deeper CI

[Continuation VI — tidal energy density of states](continuation-06/README.md) derives the static orbit-energy measure for finite Bose halos and the universal shallow-cutoff scaling. A quadratic energy cutoff yields approximately 19 times more states near escape than linear lowering in the same internal normalization, though their dimensionless masses are nearly equal. **This is not a measured evaporation rate.**

[Weekly deep theoretical numerics](https://github.com/IskanderAbdullin/Dark-Matter/actions/workflows/deep-numerics.yml) reruns seven heavier mathematical checks every Sunday at 04:23 UTC or by manual dispatch. It keeps downloadable JSON outputs for 14 days. The lightweight push/PR suite now contains ten scientific checks, plus synthetic SPARC regressions.

The follow-up [open energy-flux scaling study](continuation-06/02-generalized-open-diffusion-boundary.md) examines which boundary smoothness a power-law diffusion coefficient selects under a prescribed nonzero escape flux. Its equations are conditional transport models, not galactic evaporation measurements.


## New cosmological transfer benchmark

[Relativistically decoupled one-real-scalar 86-eV Bose relic](cosmology/README.md) — independently checked entropy, abundance, Delta Neff and correct CLASS noncold phase-space normalization, with automated [CLASS transfer workflow](https://github.com/IskanderAbdullin/Dark-Matter/actions/workflows/thermal-bose-class.yml). This is a **conditional standard thermal relic model**, not a conclusion for all Bose-stimulated cold-production scenarios or a data exclusion.


### Verified CLASS benchmark and colder-spectrum research

The public [thermal 86-eV Bose CLASS analysis](cosmology/THERMAL-CLASS-VERIFIED.md) is a successful full linear Einstein–Boltzmann calculation (not just ballistic free-streaming), giving a conditional amplitude-half point k≈0.433/Mpc and correct analytic density. A [positive colder toy phase-space family](cosmology/TOY-COLDER-PSD.md) holds the present number density fixed while varying the momentum distribution; its CLASS transfer scans test sensitivity to nonthermal initial conditions without asserting a realizable production mechanism.


## Research status: kinetic production and primordial hot-tail constraints

The [consolidated research status and falsifiability roadmap](THEORETICAL-STATUS-2026-11-21.md) integrates the latest [verified ten-case CLASS cold+hot transfer study](cosmology/VERIFIED-HOT-TAIL-CLASS.md) and the first [verified Bose-stimulated 1↔2 quantum collision benchmark](kinetics/VERIFIED-QUANTUM-BURST.md). The former demonstrates potentially severe small-scale structure suppression from even a residual hot component within an **illustrative** initial momentum family. The latter shows a number-conserving, Pauli-respecting and discrete-energy-conserving very cold burst in a fixed thermal bath, but the momentum distribution is **not yet grid converged**, and no self-consistent cosmological production is claimed.
