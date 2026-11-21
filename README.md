# Bosonic dark-matter halos with a dominant noncondensed component — research program

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
