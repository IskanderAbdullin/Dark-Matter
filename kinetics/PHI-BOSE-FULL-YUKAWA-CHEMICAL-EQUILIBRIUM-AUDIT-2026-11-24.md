# Full real-scalar Yukawa chemical equilibrium is **not** the earlier synthetic nonzero-μφ fixed point (2026-11-24 UTC)

## Scientific question

The physically normalized *restricted* 3-species χ1↔χ2+φ Yukawa numerical operator has common-T chemical equilibrium μ1=μ2+μφ. Its tests use e.g. μ1=−.15GeV, μ2=+.05GeV, μφ=−.2GeV. This is **exactly** a fixed point of the isolated particle-number-only decay network, and remains a valuable numerical conservation/detailed-balance regression. It is, however, NOT the chemical equilibrium of the natural **Dirac+antiDirac+real-scalar** microscopic theory once antiparticles and allowed pair reactions are accounted for.

The [exact 5-species stoichiometric audit](PHI-BOSE-YUKAWA-CONTINUOUS-CHARGE-VS-ACCIDENTAL-COLLISION-2026-11-24.md) has four reaction rows: χ1→χ2+φ, antiχ1→antiχ2+φ, χ1+antiχ1→2φ and χ2+antiχ2→2φ. Its integer rank is 4, nullity1; the only continuous vector chemical potential is the common Dirac charge μU1 times (1,−1,1,−1,0). A real scalar cannot carry this fermion U(1) charge, and its equilibrium chemical potential vanishes in the full pair-and-inverse network.

## Algebraic exact proof and independently tested physical occupations

For species (χ1,antiχ1,χ2,antiχ2,φ), every chemical equilibrium process must have stoichiometric affinity Ar=Σμ_initial−Σμ_final=−S_r·μ=0. With a vector Dirac U(1) and anti-particle potentials opposite, μ=(μ1,−μ1,μ2,−μ2,μφ). The restricted previously employed point

    μ1=−.15, μ2=+.05, μφ=−.2 GeV

gives **four affinities**

    (A_decay_particle,A_decay_antiparticle,A_pair_chi1,A_pair_chi2)
       = (0,+.4,+.4,+.4) GeV.

Only the *particle* decay is stationary; the anti-decay and both 2↔2 pair processes are forward-favoured (provided the actual kinematic phase spaces are populated). This is a sharply scoped **negative microphysical-consistency result** for treating the old μφ=−.2 numerical thermal fixed-point as a *universal equilibrium* of the Yukawa+anti network; it is not a bug in the restricted solver.

At a true common-T full equilibrium satisfying all four rows:

    μphi=0,    μchi1=μchi2=μU1,
    μanti_chi1=μanti_chi2=−μU1.

The pair detailed-balance relation is independently verified at actual on-shell energies E2+Eanti=k3+k4. Its quantum forward/reverse factors are

    F_pair=f2(E2) f_anti2(Eanti) [1+bφ(k3)][1+bφ(k4)],
    R_pair=[1−f2(E2)][1−f_anti2(Eanti)] bφ(k3)bφ(k4).

Using exact FD and BE odds,

    F_pair/R_pair = exp[(μ2+μanti2−2μφ)/T] = exp[−2μφ/T].

At earlier μφ=−.2,T=1 this is **exp(.4)=1.4918247**, so a physical positive pair amplitude drives a nonzero forward reaction. At μφ=0 the ratio is exactly one at ALL on-shell energies (not just integrated totals). The pair quantum entropy integrand (F−R)ln(F/R) is nonnegative in interior, tested for multiple on-shell energies and μphi.

The audit is exact **stoichiometry + quantum statistics**. It DOES NOT prove the allowed pair channel is cosmologically effective at all temperatures, or that a late-time non-equilibrium Bose zero mode is excluded: actual thermal pair rates/medium propagators, SM portal, expansion and possible conserved charges need to be derived.

## Reproduction and NEXT_ACTION

- [Pure chemical affinity and quantum gain/loss reference](phi_bose_pair_chemical_affinity.py)
- [Four focused tests](../tests/test_phi_bose_pair_chemical_affinity.py)
- [Independent finite 2→2 tree cross-section & forward MB rate](PHI-BOSE-YUKAWA-PAIR-ANNIHILATION-MB-HUBBLE-AUDIT-2026-11-24.md).
- Commands: OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_phi_bose_pair_chemical_affinity.py; python -m kinetics.phi_bose_pair_chemical_affinity

**NEXT_ACTION:** couple physical pair F/R with amplitudes and invariant 2↔2 phase-space to the existing positive dark-sector 1↔2 implicit solver, including explicitly separated particle/antiparticle occupancies and net vector U1. Quantify real pair chemical relaxation versus H including quantum blocking/Bose stimulation, medium thermal mass and finite SM+FRW. Retain the old nonzero μ fixed-point only as restricted-single-channel synthetic regression; no real 86-eV Ω/BEC or CLASS-ready inference.
