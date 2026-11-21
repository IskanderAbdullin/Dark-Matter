# Bosonic dark-matter halos with a dominant noncondensed component — research program
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
