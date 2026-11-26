# Exact-shell lattice resonance and graph connectivity: analytic negative control (2026-11-23)

## Question
Can an exact-energy uniform lattice, by itself, certify thermalisation in the proposed physical finite-volume Yukawa strong operator?

## Exact arithmetic
Let the three energy lattices have common spacing h and offsets a1,a2,aφ, with E1_i=a1+ih, E2_j=a2+jh, Eφ_k=aφ+kh. Exact on-shell events exist only if (a1-a2-aφ)/h is an integer. When it is an integer m, every allowed event satisfies i=j+k-m. If the offset mismatch is not an integer, **there are zero exact triplets** irrespective of positive capacities or time integrator. Thus a generic common-spacing lattice has a measure-zero resonance condition on its offsets. Even when resonant, a finite interval, physical mass thresholds, and angular kinematics can remove most triplets.

For distinct spacings h1,h2,hφ, exact events require a1-a2-aφ+j h2+k hφ-i h1=0. Incommensurate ratios can leave isolated or no triplets. This is an arithmetic obstruction, not a statement about the continuum collision kernel.

## Nullspace test before time evolution
For each allowed event a=(i,j,k), form stoichiometric column s_a=(-e1_i,+e2_j,+eφ_k). A linear conserved quantity c·n is exactly invariant iff cᵀs_a=0 for every event. With N species modes and rank r of S, the number of independent linear invariants is N-r. At minimum, two charges (N1+N2 and N1+Nφ) and total energy are independent in a nondegenerate physical discretisation, hence nullity >=3. If nullity >3, the finite event graph has additional frozen linear combinations and cannot uniquely establish the expected chemical/thermal equilibrium. Conversely nullity=3 is necessary but not sufficient for continuum accuracy, ergodicity, or the correct weak source.

**Minimal counterexample:** two disjoint exact triplets (1a,2a,φa) and (1b,2b,φb), each with positive event weights and exact shell energy. There are six modes and two independent columns, hence nullity=4. Both components separately conserve their own charges. Each event satisfies the H theorem and quantum bounds in the interior, but the system cannot exchange populations between components. A thermal-looking stationary state does not prove global thermalisation.

## Consequence for the next physical construction
Do not choose a lattice by forcing exact resonance and then interpret its graph as a converged physical strong collision operator. Instead derive the common invariant positive event measure and species phase capacities from continuum phase space, audit the graph rank/connectivity, and compare the summed weak number and three energy moments against the independent exact-angle and dilogarithm references at multiple UV/IR and radial/angular/daughter resolutions. A finite-bin energy-conserving interpolation must also demonstrate microreversibility and a discrete H theorem; conservation of mean energy alone is inadequate.

This note is an **analytic negative control**, not a numerical run or an implementation of the strong solver. No new claim about 86-eV relic abundance, BEC, or CLASS follows.

References within repository: kinetics/YUKAWA-STRONG-DISCRETIZATION-GATE.md; kinetics/PHI-BOSE-THERMAL-EXACT-WEAK-REFERENCE-AUDIT.md.
