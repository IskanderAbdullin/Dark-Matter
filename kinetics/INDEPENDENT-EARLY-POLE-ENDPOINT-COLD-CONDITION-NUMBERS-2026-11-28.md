# Percent accuracy of the cold coefficients needs much tighter pole control

Use the fine D2 charge and scalar occupations at the held transverse pole k=.0015. For each of q=3.5e−6 and 4.7e−6, perturb only the scalar pole mass or only the photon energy, independently holding the photon residue fixed. Use paired (pm h,pm h/2) log offsets with h=2.5e−5. The two steps and Richardson extrapolate agree within 4.35e−6 relative; independently paired radial/azimuth orders agree within 4.48e−13 on the point coefficients. The original larger h=1e−4 did not pass the slope gate and is not used.

| q | Log-mass slope of gross forward | Log-energy slope of gross forward | Log offsets for 1% changes (mass, photon energy) |
| --- | ---: | ---: | --- |
| 3.5e−6 | +266.6336 | −710.0078 | 3.73e−5, 1.40e−5 |
| 4.7e−6 | +263.1914 | −700.9376 | 3.78e−5, 1.42e−5 |

The inverse slopes are larger in magnitude: +284.8 to +290.9 for the gap and −774.0 to −757.9 for photon energy. The log offsets are local first-order scales, not certified accuracy requirements for the coupled model. Holding the original residue and evaluating the two endpoints of the previous charge-moment photon-energy interval [0.00230413885981,0.00230438226143] changes the forward coefficient at q=3.5e−6 by +2.597% and −4.822%. At q=4.7e−6 the changes are +2.563% and −4.763%. Thus a small frequency interval can still matter even though the actual photon shift from one fitted evolution to another was much smaller. The endpoint evaluation is a controlled sensitivity sample, **not** a rate enclosure that propagates occupation and residue errors or a claim about the physical photon distribution.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error -m unittest test_early_pole_endpoint_cold_condition_numbers
python -W error early_pole_endpoint_cold_condition_numbers.py \
  results/early-pole-cell-endpoint-fine-d2-x05-state-2026-11-28.json.gz \
  results/early-pole-endpoint-cold-threshold-enclosure-2026-11-28.json
```
