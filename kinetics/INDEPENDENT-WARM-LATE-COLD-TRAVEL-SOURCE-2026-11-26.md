# Select a travelled distance instead of a temperature label

The [photon-heating calculation](INDEPENDENT-WARM-LATE-COLD-SCALAR-BALLISTIC-2026-11-26.md)
showed why p/T<.1 alone is not a dynamical definition
of cold. Here specify how far a collisionless scalar
can travel from its creation to equality in the same
reference background. This is a distance diagnostic,
not a matter power scale or a halo prediction.

The comoving path increases with conserved momentum q
and decreases when birth is later. For a conservative
upper on the allowed injection momentum at any time
from 50 to .1 MeV, discard all travel before the final
EOS knot, a=2.35e−8, later than the latest considered
birth. The effective radiation contribution is constant
after this knot. On each scale-factor interval use the
larger right-end denominator in the positive path
integrand, and integrate the remaining speed factor
analytically with asinh. This gives a lower bound on
path length at fixed q and therefore an upper bound on
the q eligible for a given maximum distance. Refinement
and an independent full-time numerical root establish
the direction and size of this deliberately permissive
cut. The held archive count uses its actual numerical
cut, not the more permissive production selection.

| Maximum reference path | Archive p/T at this cut | Direct lepton upper / held count | Conditional lepton+pion upper / held count |
| --- | ---: | ---: | ---: |
| .0556 Mpc | .005862 | 9.0e−43 | 3.08e−36 |
| .1000 Mpc | .011247 | 3.61e−20 | 5.76e−16 |
| .3000 Mpc | .038602 | .006791 | .006792 |
| 1.0000 Mpc | .153327 | >1 (inconclusive) | not evaluated |

Fold the same direct vacuum lepton source with the
[off-pole and cold-angle bounds](INDEPENDENT-WARM-LATE-COLD-VISIBLE-OFFPOLE-2026-11-26.md).
The separate pion allocation includes the chosen KLOE
fixed-bin values plus twice their marginal errors, the
two independent [positive low-gap allowances](INDEPENDENT-WARM-LATE-COLD-PION-DISPERSIVE-RECOIL-2026-11-26.md),
and a conditional uniform R_pi≤1 above the measured
table. High-band results are upper bounds per unit
R_pi and can be rescaled without moving the other rows.
The most stringent path rows are dominated by deliberately
loose upper high-band terms, not by a detected pion
population there. At .3 Mpc, the lepton bound dominates
the resulting 0.68% ceiling.

Neither source bound includes hadronic channels other
than the conditional vacuum charged-pion model, an
initial heavy population, dark rescattering, or
changes to the expansion or perturbations. The archive
denominator selects an existing subset for comparison,
not a measured cosmological abundance. The .0556 Mpc
row is only a comparison to the earlier thermal path
diagnostic, not a Lyman-alpha exclusion.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
D=results
python -W error warm_late_cold_distance_source_reference.py EOS.dat \
 "$D/warm-multishell-corrected-q64-fine-stage-48-2026-11-26.json.gz" > paths.json
python -W error warm_late_cold_pion_distance_bound_reference.py EOS.dat paths.json \
 "$D/warm-late-cold-pion-dispersive-recoil-2026-11-26.json" KLOExSec.data
~~~

[Path and lepton limits](results/warm-late-cold-distance-source-2026-11-26.json);
[conditional pion transfer](results/warm-late-cold-pion-distance-2026-11-26.json).
