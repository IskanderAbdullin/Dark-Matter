# Recompute the partial direct light gain at the corrected late anchor

Repeat the [direct pair subtotal](INDEPENDENT-WARM-LIGHT-DIRECT-PAIR-DOSE-2026-11-26.md)
at the exact q64 7.394049 MeV state in the
[matched corrected trajectory](INDEPENDENT-WARM-CORRECTED-TIME-AND-SUPPORT-2026-11-26.md).
Reuse the original conservative gross available-energy budget,
transported by the radiation entropy clock. Keep the heavy per-spin
cap at one-half as a hypothesis and drop all direct-pair losses.
For the b⁻³, b⁻⁴ and b⁻⁶ gain weights, independently check that the
local EOS interval floors improve over the original uniform floor.

The light occupancy of represented modes at the anchor is at most
4.1021e-5. On the refined clock the positive scalar and selected
visible pair contributions to the integrated occupation bound are
0.05018 and 0.01027. The heavy-to-light term can still contribute
7.667 if granted the full energy envelope. A hypothetical uniform
heavy energy fraction 0.10 reduces the whole direct subtotal to
0.41391. It reaches the remaining half-cap slack at a heavy energy
fraction of about 0.11772 (0.11770 on the uniform clock). This is
a sufficient threshold for this subtotal alone, not an inferred
heavy fraction or light-cap proof.

The older [finite comoving band envelope](INDEPENDENT-WARM-HEAVY-COMOVING-BAND-ENERGY-2026-11-26.md)
can be inherited only if its original premises held on the entire
interval before the later anchor. In particular, an initially empty
old band does not make the late band empty. Its comoving cut rescales
by the relative expansion 1.4128443. Keep the worst of that initial
allowance and the newly integrated source over the conservative loss
coefficient. The kept fraction must be subtracted from the partial
direct threshold before allocating any budget above the cut.

| Cut at the late anchor (GeV) | Below-cut fraction including inherited bound | Remaining fraction for above-cut modes, direct subtotal only |
| ---: | ---: | ---: |
| 0.42468 | 0.005254 | 0.11247 |
| 0.70779 | 0.028122 | 0.08960 |
| 0.99091 | 0.095733 | 0.02199 |

The remaining fraction is a required target, not a measurement
or constraint on the high band. Unknown elastic transfer into the
heavy band, omitted incoming modes, prompt decay/inverse gain and
light elastic redistribution are not bounded here. The earlier
CLI's direct-pair JSON was reproduced byte for byte after extracting
its reusable gain formula, and the default earlier radiation-clock
comparison still gives 1.121197 and 2.111328.

~~~sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error \
  warm_corrected_later_light_gain_headroom_reference.py EOS.dat \
  results/warm-multishell-corrected-q64-stage-24-2026-11-26.json.gz
~~~

[Local clock and transported band controls](results/warm-corrected-late-direct-light-gain-2026-11-26.json).
