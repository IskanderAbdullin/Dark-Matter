# An arbitrary-spectrum vacuum bound on real-scalar quartic multiplication

For a real scalar of mass μ=86 eV, include the **tree diagrams with
two** λφ⁴/4! **vertices** for φφ→4φ. Do not assume a thermal spectrum,
a translationally coherent field, or a kinetic-equilibrium closure.
The result bounds the vacuum rate of this single interaction, or
its Bose-stimulated rate **only** with a uniform pointwise upper
bound F on all four final occupations over the whole period.
For a complementary thermal treatment of this channel see
[Arcadi, Lebedev, Pokorski and Toma](https://arxiv.org/abs/1906.07659).

## Positive phase-space bound

With all six legs labeled, the tree amplitude sums ten partitions
of the legs into distinct unordered triples. Cauchy bounds the
amplitude squared by ten times the sum of squared individual
diagrams (times λ⁴). There are four timelike partitions containing
the sum of three final momenta, and six mixed partitions in which
one incoming momentum is subtracted from a final pair.

For a timelike partition with triple invariant A≥9μ², factorize
the four-body phase space into two-body and three-body parts. Use
the upper measures Φ₂≤1/(8π), Φ₃(A)≤A/(256π³), and integrate A
with measure dA/(2π). The squared propagator is (A−μ²)⁻².
The result is at most

`I_time = [log((s−μ²)/(8μ²)) + 1/8 − μ²/(s−μ²)]/(4096π⁵)`.

For a mixed partition regard the two final pairs as clusters with
invariants A,B≥4μ². The transfer is spacelike: if it were timelike
and future pointing it would require a single incoming scalar to
decay to a heavier pair plus a positive-mass state; a past-pointing
transfer similarly contradicts the other vertex. If t₊,t₋ are
the exact endpoints, their product is
`(A−μ²)(B−μ²) + μ²(A−B)²/s`.
The exact angular mean of (μ²−t)⁻² is
1/[(μ²−t₊)(μ²−t₋)], thus no larger than
1/[(A−μ²)(B−μ²)]. Factor the 2+2 phase space, bound each
two-body measure by 1/(8π), and extend the cluster-mass domain
to a rectangle to get

`I_mixed = log²((s−μ²)/(3μ²))/(2048π⁵)`.

Consequently the physical center-of-mass rate is at most
`λ⁴(10/24)(4 I_time+6 I_mixed)/s` for s≥16μ², with 1/24
for the identical final scalars. In the bath the Møller rate
is smaller by s/(4E₁E₂)≤1. Maximizing the simple logarithmic
upper terms gives a uniform bound K at all incident momenta.
An independent four-body Sobol integration includes all diagrams
and their interference at six energies and lies below the analytic
bound. At high energy its scrambles differ due to the sharp angular
regions; it is a cross-check, not the proof or a precision rate there.

## Optimistic cosmological comparison

Allow every initial dark particle to be scalar immediately, and
suppress the 4→2 reverse reaction. If the initial total yield is Y₀,
then dY/dt ≤ K(1+F)⁴ s_visible Y² (two net particles per event times
the 1/2 identical incoming pair factor). With adiabatic visible
entropy and H at least the photon floor, integrate this inequality
to obtain Y≤Y₀/(1−τ) when τ<1. In the prescribed-source comparison,
place **all** the independently estimated source particles at
reheating, which maximizes the time available for this number
growth. The assumed cumulative source is not itself bounded by
this calculation if new interactions are added upstream.

Use the deliberately generous [project halo and curvature comparison](https://github.com/IskanderAbdullin/Dark-Matter/blob/90c141722376dab64e4cd10422d545474817787c/kinetics/PHI-BOSE-QUARTIC-CURVATURE-HALO-SCATTERING-TENSION-2026-11-24.md)
for the two reference λ values; these are conditional inputs, not
a new constraint from a halo likelihood.

| Input and quartic | Vacuum growth τ upper | Scalar Ωh² upper with all particles scalar |
| --- | ---: | ---: |
| 12 MeV archive, λ=7.652e-8 | 2.39e-8 | 5.250e-4 |
| 12 MeV archive, enlarged λ=1.085e-6 | 9.68e-4 | 5.255e-4 |
| All prescribed TR50 counts placed at reheat, λ=7.652e-8 | 1.34e-7 | 5.271e-4 |
| Same TR50 total, enlarged λ=1.085e-6 | 5.41e-3 | 5.299e-4 |
| All prescribed TR9.019 counts placed at reheat, enlarged λ | 7.13e-7 | 5.160e-7 |

With the smaller λ and whole prescribed TR50 count, reaching 99% of
0.120 requires a final-occupation cap greater than about 51 for this
inequality even to permit it. For the larger λ the corresponding
necessary threshold is about 2.68. These are **not** sufficient
conditions, nor estimates or upper bounds on actual occupations.
An occupation at an inaccessible momentum does not enhance the rate.

This is a strict conditional gate for a subset of the model. It does
not include an effective φ⁶ contact, Yukawa loops, a broken-symmetry
cubic, heavy-pair recycling, coherent or nonperturbative dynamics,
extra injection, or arbitrary Bose enhancement. It is not a
stand-alone exclusion of scalar dark matter.

```sh
OPENBLAS_NUM_THREADS=1 python -W error scalar_quartic_vacuum_multiplication_bound_reference.py EOS.late.dat results/conditional-relic-normalization-2026-11-26.json
```
