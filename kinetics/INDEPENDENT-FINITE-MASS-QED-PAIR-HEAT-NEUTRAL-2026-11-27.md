# Zero scalar number exchange does not imply zero photon heat exchange

On the [hard charged-pair crossing](INDEPENDENT-FINITE-MASS-QED-CROSSED-ANNIHILATION-SOURCE-2026-11-27.md),
hold the charged and scalar occupations, masses, and cold scalar
momentum fixed. Keep the same explicit anti-daughter proxy and
one charge assignment. Change only the temperature of the photon
distribution. This is a sensitivity comparison, not a candidate
temperature for the complete plasma.

For the signed scalar source divided by its held occupation, write
$J(k)=(1+1/n_s)A + n_\gamma[(1+1/n_s)A-B]$, where $A$ and $B$
contain the forward charged product and reverse Pauli product,
respectively. The temperature-independent part is positive. The
coefficient of the increasing photon Bose occupation is negative
at every allowed event for these archives. Hence the scalar number
moment and its positive-$k$ photon energy moment each decrease
strictly with the imposed temperature. Each changes sign between
vacuum and the archived temperature fixture, giving unique but
not necessarily equal zeros. At the scalar number zero,
weighting the same positive and negative contributions by
photon energy need not cancel them.

| Archive | Scalar number zero | Photon heat zero | Photon heat per held scalar occupation at the number zero |
|---|---:|---:|---:|
| fine | 0.49937298 | 0.52507558 | +0.15859248 |
| third | 0.50481099 | 0.52954867 | +0.15630340 |

The temperatures and last column are in archive units and omit
couplings. The last column uses the unweighted scalar occupation
as a normalization; it is not an energy fraction of a bath.
The number integral is still negative at each heat zero.
The largest two-mesh root drift is 3.58e−11. For the
two [hard energy comparators](INDEPENDENT-FINITE-MASS-QED-HARD-REAL-PHOTON-HEAT-2026-11-27.md),
the positive photon-energy weighted omitted tails at temperature
1 are bounded by 4.79e−22 (inverse charged pair, k>68)
and 3.39e−8 (parent plus scalar inverse, parent energy above 56)
before coupling factors. These follow by adding an energy factor to
the existing polynomial-times-Boltzmann bounds, not by treating the
number tail as a heat tail.

Neither moment zero ensures pointwise detailed balance. In contrast,
the separate [common-temperature real heat gate](INDEPENDENT-FINITE-MASS-QED-SIGNED-REAL-PHOTON-HEAT-2026-11-27.md)
uses consistent equilibrium weights for every species and each
real direction. No self-consistent transport solution or
thermal history is claimed.

Replay from kinetics:
~~~sh
OPENBLAS_NUM_THREADS=1 python -W error finite_mass_qed_pair_heat_neutral_reference.py
~~~
[Fine receipt](results/finite-mass-qed-pair-heat-neutral-fine-2026-11-27.json),
[third receipt](results/finite-mass-qed-pair-heat-neutral-third-2026-11-27.json).
