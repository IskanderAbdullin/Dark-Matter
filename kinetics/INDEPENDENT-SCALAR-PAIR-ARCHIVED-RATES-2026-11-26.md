# Stable fermion annihilation to a scalar pair on archived distributions

Evaluate the real off-diagonal scalar Yukawa channel χ₂χ̄₂ ↔ φφ at
fixed distributions. The inputs are represented two-moment cells,
with the exterior occupation set to vacuum. The code does not assume
a common temperature or evolve either distribution. Produced particles
outside the occupied cells are tallied, not reflected at the grid edge.

## Independent checks and conventions

The sum of t- and u-channel amplitudes can be written as a slash-vector
plus a scalar propagator term B. If the incoming fermion momentum in
the center of mass is along z, its spin-averaged square reduces to the
sum of two positive pieces, proportional to
E² a_perp² + (m a_z + p B)². This avoids cancellation near threshold.
An independent four-by-four Dirac trace agrees within 5.3e-13 at the
sampled energies and angles. For equal exchanged and incoming fermion
masses and a massless scalar, the low-velocity result approaches
3 y⁴ v_rel²/(128 π m²), as in [Biondini and Shtabovenko, Eq. 4.17](https://arxiv.org/html/2106.06472).
That normalization control is not used in place of the actual unequal
masses. With real FD/BE populations on a shared manufactured thermal
state, forward and reverse rates cancel to 5.3e-7 relative or better
in number and energy. The inverse computation includes four outgoing
spin states and the initial identical-scalar pair factor. The forward
calculation includes the final identical-scalar factor.

Final-state stimulation and Pauli blocking are integrated in the lab
frame with the angle split at both daughter shell boundaries. The net
shell output gives number and kinetic-energy rates for each species,
including both exterior bins. Summing restores conserved total number
and total energy; individual species numbers are not conserved.
The normalization uses the visible-radiation Hubble rate at the
reported temperature.

## Fixed-state results

| Initial archive | Vacuum mean σv (GeV⁻²) | Inverse events / stimulated forward | Net fermion loss / (nχ H) | Net scalar gain / (nφ H) | Mean forward scalar energy |
| --- | ---: | ---: | ---: | ---: | ---: |
| 50 MeV reheat, 12 MeV endpoint | 1.3784e-12 | 0.7589 | 0.3371 | 0.8768 | 43.7 MeV |
| 50 MeV reheat, 9.3 MeV endpoint | 1.3373e-12 | 0.3822 | 0.6508 | 1.6980 | 37.6 MeV |
| 9.019 MeV reheat, 3.157 MeV vacuum-source endpoint | 7.2580e-13 | 2.43e-11 | 0.0002591 | 57.80 | 25.3 MeV |

The high-reheat net is appreciably smaller than the gross forward
rate, **not** because of small stimulation, but because the existing
hard scalar population supports inverse reactions (76% and 38%
of the corresponding forward counts). The forward Bose increase is
only 0.00370% and 0.00252% in these snapshots. The 12 MeV populated
vacuum q32 control and its exact same-inner-state q16 crop give
0.333409 and 0.333446 for the normalized all-output stable net loss.
This crop check concerns the same fixed incident inner state, not
arbitrary nonthermal tails at another time.

At low reheating, even the small instantaneous stable number loss
would be a sizable source relative to the tiny preexisting scalar
count: about 58 times that count per Hubble time at this snapshot.
It produces hard scalars, with no represented direct injection below
0.25 MeV in these cases. The low-reheat forward count landing outside
the source-only grid is roughly 0.20%. These rates are not integrated
survival fractions or a new relic estimate. A coupled expanding run
with this operator is needed before revising the scalar distribution
or repeating a condensation claim. The separate conservative
[vacuum-only bound](INDEPENDENT-ARCHIVED-VACUUM-DEPLETION-2026-11-26.md)
continues to apply only under its stated assumptions.

## Reproduce

Prepare [the EOS inputs](INDEPENDENT-BENCHMARK-SETUP.md) and run from
`kinetics`:

```sh
export OPENBLAS_NUM_THREADS=1
python -W error archived_scalar_pair_reaction_reference.py EOS.dat results/{12-quantum-48-shape,early-9p3-quantum-shape,wide32-embedded-vacuum-shape}-2026-11-26.json --crop-last 16
python -W error archived_scalar_pair_reaction_reference.py EOS.late.dat results/low-reheat-matched-vacuum-shape-2026-11-26.json --eos-lower-knot-limit .001
```
