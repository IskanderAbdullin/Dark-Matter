# One bounded shape change can reverse both heavy stress comparisons

The moment-hull atomic witnesses establish separate pressure and pseudo-pressure ambiguities at heavier formal masses. This check uses no atoms and no per-observable optimization. For each original warm cell, define two bounded distributions in its normalized reference excess-energy coordinate. The first uses a uniform interval around the stored mean; the second puts a uniform interval at each end and assigns their two weights to preserve that exact mean. Both retain the same exact total count. Use those same two densities throughout the full probe-mass sweep.

With maximum interval width 0.20 in the normalized energy coordinate:

| Formal mass | G4 center vs G3 spread (pressure, pseudo-pressure) | G4 spread vs G3 center |
| ---: | ---: | ---: |
| Archived reference | −0.125385%, −0.125574% | −0.125395%, −0.125594% |
| 1 | −0.229313%, −0.585792% | −0.073621%, +0.201925% |
| 10 | −1.90882%, −8.38302% | +1.40518%, +7.91655% |
| 100 | −2.34200%, −11.9375% | +1.82033%, +12.1653% |

Thus at 10 and 100 the same number/energy-preserving change can reverse both comparisons simultaneously, even with finite-width distributions. This family does not reverse pressure at 1 (where the separately optimized atomic witnesses do). The held measurement-mass sign is unchanged in both families, consistent with the independent hull. All values reproduced by paired quadratures to better than 3e-12 relative. None is a claim about an evolving cosmological distribution, or that these arbitrary replacements are dynamical solutions.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
p=results/early-pole-fine-charge2-scalar4; d=2026-11-28
python -W error -m unittest test_early_pole_endpoint_warm_bounded_shape_witness
python -W error early_pole_endpoint_warm_bounded_shape_witness.py \
 --native "$p"-x05-state-$d.json.gz --gauss "$p"-gauss{3,4}warm-x05-state-$d.json.gz
```
