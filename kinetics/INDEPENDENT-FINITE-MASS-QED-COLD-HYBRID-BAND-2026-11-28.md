# The integrated cold discrepancy mostly follows the dispersion pair

The [point tests](INDEPENDENT-FINITE-MASS-QED-COLD-SENSITIVITY-2026-11-28.md) isolate scalar pole mass, photon dispersion and charged populations at fixed q. This calculation weights the full cold interval [1e−6,1.5e−5] at k=.0015. Use all four combinations of the charged histogram and the **coupled dispersion pair** (scalar pole mass and photon frequency/residue). For each geometry use both scalar histograms as weights on an identical outer scalar q mesh. The scalar grid is refined at the union of both sets of scalar edges and the exact thresholds of the chosen charged support. Only the all-fine and all-third arms are native archive combinations.

The two quadratures for each geometry check the full and centered scalar energy separately. Both native arms replay the archived cold number, inverse and energy coefficients to 1.11e−13 relative or better. The other six arms are counterfactuals, not solutions of another coupled plasma model.

For each positive extensive coefficient, changing the three factors from third to fine in different orders gives the same total change but different individual marginal changes. Average the log changes over all six orderings; the following are the exponentials of those averaged effects.

| Integrated quantity | Charged occupations | Dispersion pair | Scalar histogram |
| --- | ---: | ---: | ---: |
| Gross forward scalar-number coefficient | 1.00831 | 7.17675 | 1.04527 |
| Gross inverse scalar-number coefficient | 1.04862 | 8.46450 | 1.05006 |
| Conditional number-neutral photon n | 1.03997 | 1.17943 | 1.00458 |

The pointwise prioritization survives the full archived cold weighting: the dispersion pair dominates the measured difference in absolute coefficients. Much of that common effect cancels in the conditional neutral ratio, so coefficient convergence and neutrality convergence are distinct. The separate mass and photon effects are not disentangled by this integrated calculation. These results do not assign the actual n, cover larger q or other k, or quantify the full reaction rate.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error finite_mass_qed_cold_hybrid_band_receipt_pilot.py \
  results results/finite-mass-qed-fixed-pole-cold-scalar-band-2026-11-28.json
# To regenerate each original panel log:
# python -W error finite_mass_qed_cold_hybrid_band_pilot.py CHARGED DISPERSION MESH RAW_PATH
```
The [receipt](results/finite-mass-qed-cold-hybrid-band-2026-11-28.json) is reconstructed from the eight compressed panel logs for fine/third charged inputs, fine/third dispersions, and meshes 0/1.
