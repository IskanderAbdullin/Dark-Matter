# Stress-test a high inverse source when only light occupancy and energy are given

The [late direct-gain recalibration](INDEPENDENT-WARM-LATE-DIRECT-LIGHT-GAIN-HEADROOM-2026-11-26.md)
leaves an above-cut heavy allowance of about 8.96% of the total
energy envelope at the late cut 0.70779 GeV, under its stated
inherited-band assumptions. Can an inverse-decay source-to-decay
comparison show a small high band using only the total available
energy and the pointwise light cap, without relying on the
archived resolved momentum distribution?

For a deliberately nonthermal diagnostic, spend one quarter of
the available energy in frozen light modes at f=1/2, and one
quarter in a narrow scalar packet. Place the packet centers at
unequal daughter energies whose sum is the chosen parent energy.
Pay the exact energy cost of each packet in the massive radial
measure; leave the other half of the budget unused by the
prescribed daughters. Verify throughout each occupied rectangle
that the on-shell incoming angle is allowed. All resulting parent
energies are above the matched comoving-band cut.

The total positive inverse energy source in that band is the
exact two-dimensional integral in daughter energies, independently
rechecked by a positive Gauss rule. For **fixed** daughter packets
the heavy collision equation is affine in its occupation. The
sum of inverse blocking and direct depletion at a tagged parent
is bounded above by (Γ M/E_H)(1+⟨f_scalar⟩). Enclosing the
scalar angular average and taking the slow variable to be total
heavy energy in the band gives `dY/dt >= Q - lambda_upper*Y`.
With initially zero heavy modes, the frozen comparison floor
approaches `Q/lambda_upper`. It is only a source-to-loss
diagnostic with prescribed input, not self-consistent evolution.

| Synthetic parent energy center (GeV) | Frozen source / upper loss coefficient / available energy |
| ---: | ---: |
| .75 | 1.054e-5 |
| 10 | .001814 |
| 30 | .01543 |
| 100 | .10556 |
| 200 | .18628 |
| 300 | .21701 |

Thus this route cannot infer a universally tiny inverse reservoir
from only the light occupancy cap and total energy. Momentum
location matters as well as the density. The synthetic incident
packets at hundreds of GeV are far outside the represented range
and are not evidence those particles are ever produced. Energy
exchange would deplete the fixed daughters, change the source,
and require a self-consistent analysis. This test neither falsifies
a dynamically derived small heavy fraction nor covers extra
scattering channels or high-energy changes to the model.

## Exact contrast with represented late snapshots

In the corrected q32 endpoint the largest represented daughter
energy sum is .431 GeV, so the band is geometrically inaccessible.
In the independently advanced q64 endpoint it is open (maximum
sum .860 GeV). The reconstructed q64 bare inverse energy rate
with stationary blocking is only 4.130e-32 in units of retained
dark energy times H. The corresponding frozen stationary heavy
energy in this band is 2.639e-45 of retained dark energy. Two high
quadrature rules agree at displayed precision; this does not
assess a true heavy solution or untracked incident distributions.

~~~sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error \
  warm_high_inverse_pressure_test_reference.py results/warm-corrected-anchor-signed-floor-2026-11-26.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error \
  warm_corrected_archive_high_parent_recheck_reference.py EOS.dat \
  results/warm-corrected-anchor-signed-floor-2026-11-26.json \
  results/warm-corrected-q32-time48-step48-2026-11-26.json.gz \
  results/warm-multishell-corrected-q64-stage-24-2026-11-26.json.gz
~~~

[Geometry and frozen reservoir controls](results/warm-high-inverse-frozen-pressure-test-2026-11-26.json).
[Represented high-band snapshot check](results/warm-corrected-archive-high-parent-recheck-2026-11-26.json).
