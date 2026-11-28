# Warm-only subdivision isolates and reduces the geometric-node moment bias

The [matched numerical moment comparison](INDEPENDENT-EARLY-POLE-ENDPOINT-SCALAR-TRANSPORT-AND-NODAL-MOMENTS-2026-11-28.md) finds that the saved scalar histogram and its geometric-node surrogate have the same total count but different higher moments. Hold the evolved S4 x=.05 occupation histogram, scalar mass, and charge cloud fixed. Duplicate each warm parent's constant occupation onto all its children, preserving every native edge. Split only *complete* native cells whose lower edge is ≥.01. This leaves both the 4-way cold scalar refinement below 1e−5 and the two-subcell charged charts unchanged.

| Warm subcells per selected native cell | Total scalar bins | Nodal q² bias | Nodal pressure-numerator bias |
| ---: | ---: | ---: | ---: |
| 1 | 145 | −8.3474% | −3.7717% |
| 2 | 167 | −2.2149% | −.9770% |
| 4 | 211 | −.5624% | −.2472% |
| 8 | 299 | −.1411% | −.06265% |

The exact piecewise-constant moment totals do not change under these divisions (relative numerical difference <5e−16). The change in the nodal estimates is therefore a representation effect, not a change in the frozen distribution. It does not predict the error or reaction rate of a dynamically evolved warmer grid. Unsplit intermediate shells remain part of every estimate.

The opt-in startup/continuation driver accepts `--warm-floor .01 --warm-div N` when constructing a **new** chart, and checkpoints bind these options to prevent restarting on a different grid. The old default still subdivides only the cold shells. The optional `--fast-capacity` reuses the [exact capacity guard](INDEPENDENT-EARLY-POLE-ENDPOINT-CAPACITY-FAST-PATH-2026-11-28.md), with the original operation still the default and the short identical-input CLI tested byte-for-byte. The [frozen fixture](results/early-pole-endpoint-warm-nodal-fixture-2026-11-28.json), not a full warm evolution, is reproduced from `kinetics/` by:
```sh
OPENBLAS_NUM_THREADS=1 python -W error early_pole_endpoint_warm_nodal_fixture.py \
  results/early-pole-fine-charge2-scalar4-x05-state-2026-11-28.json.gz
```
