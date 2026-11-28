# Frozen factor changes across the fine cold-band interior

The [integrated fine d1 to d2 comparison](INDEPENDENT-EARLY-POLE-ENDPOINT-FINE-D2-COLD-BAND-2026-11-28.md) does not identify a mechanism. Here the same saved states are crossed at four held q values: 2.5, 3.5, 4.7 and 8.2 in units of 1e−6. Independently switch the fitted scalar pole mass, the transverse photon energy and residue as one factor, and the charged occupation callback and support. Compute all eight corners. Average their logarithmic changes over all six replacement orders; the minimum and maximum over orders are kept. Factor out the corresponding scalar Bose multiplier or absorption occupation exactly.

The mass changes from 8.78487705e−4 to 8.72936483e−4 (−0.632%). Its isolated average gross-forward factors are 0.1999, 0.2038, 0.2102 and 0.2332. Photon-pole factors are 0.9972–0.9975; switching only the continuum charged callback and support changes the gross-forward coefficient by roughly 1e−8 relative at the tested points. The population multiplier is 0.9665–1.0538. At q=3.5e−6, all factors give 0.20018 gross forward and 0.17371 inverse relative to d1. The tiny charged-callback factors are conditional on the two constructed closures and are not an error bound on real charge distributions.

This is numerical evidence that the gap, rather than a direct change in the charged callback, dominates among the tested **pointwise** replacements. It does not attribute the integrated comparison or prove a continuum limit. The two angular and radial quadratures agree on each of the 32 crossed cases within 8.47e−13 relative. A crossed case is not necessarily a self-consistent physical state.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error -m unittest test_early_pole_endpoint_cold_point_factorial
python -W error early_pole_endpoint_cold_point_factorial.py \
  results/early-pole-cell-endpoint-fine-x05-state-2026-11-28.json.gz \
  results/early-pole-cell-endpoint-fine-d2-x05-state-2026-11-28.json.gz
```
The [receipt](results/early-pole-endpoint-cold-point-factorial-2026-11-28.json) includes bare corner data, paired differences, and order sensitivities.
