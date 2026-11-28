# Most resolved warm-to-cold feedback acts through the hard populations

Extend the [short shared-grid comparison](INDEPENDENT-EARLY-POLE-JOINED-EXTERIOR-PAIRED-COLD-2026-11-27.md) on the fine final n216 snapshot. Evaluate each cold source increment by replacing the baseline hard state, then the scalar pole, then the existing cold populations with their promoted-arm values. This telescopes to the direct evolved cold difference. The order defines the finite counterfactual pieces; they are not separate physical channels.

| Native hard state, n6, exposure 2e−6 | Direct cold feedback | Hard subtraction piece | Independent anchored pole piece |
| --- | ---: | ---: | ---: |
| q=.1–5 markers | −7.49e−19 | −7.37e−19 | −4.33e−21 |
| Same eleven shells, energy-matched q4 | −5.57e−19 | −7.06e−19 | −4.12e−21 |
| Broad marker support q>1.00001e−5, through 5 | −1.08e−17 | −9.37e−18 | −1.26e−18 |

For the last column, differentiate the common gap equation at fixed comoving momenta with respect to the explicitly evolved count difference. Multiply that mass shift by the centered slope of the evaluated cold packet rate. The two relative mass perturbations, 1e−4 and 5e−5, give pole pieces within 1.3% of one another. This independent small-difference estimate avoids subtracting two nearly identical independently solved pole values. The directly subtracted pole piece is particularly unreliable for the resolved narrow shell (+1.50e−19 with the opposite sign to the anchored slope). Replacing it *as a sensitivity diagnostic*, not a new trajectory, gives about −7.10e−19 total. These calculations support hard-population dominance but not a precise percent difference between narrow shell representations or an accurate direct pole signal there.

The native broad support retains negative feedback at n12 (−9.42e−18), although the n6→n12 movement is about 13%; neither is a fine-limit estimate. On the continuous hard reconstruction and narrow marker support, n6/n12 give −7.68e−19/−7.83e−19. The n12 result changes by only about 1.0e−23 when the optional per-step energy gate is tightened. To test the numerical resolution, a manufactured multiplier on reciprocal warm collision packets (both scalar and charged arms, never the initial occupation) yields cold response per unit multiplier of −7.87e−19/−7.91e−19 at x16/x64; x0 gives zero warm gain and a residual cold signal of only −2.80e−24. The broad n6 response scales to −1.07e−17 per unit multiplier at x16, within about 1.5% of x1. The multiplied runs are not physical alternatives.

From `kinetics/`, with `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1` and `a=results/early-pole-joined-fine-216-memory-2026-11-27.json.gz`:
```sh
p=early_pole_live_scalar_exterior_paired_difference_reference.py
python -W error "$p" "$a" --joined --duration 2e-6 --steps 6 --order 8 --attribute-feedback --stable-pole-tangent
python -W error "$p" "$a" --joined --duration 2e-6 --steps 6 --order 8 --attribute-feedback --stable-pole-tangent --warm-order 4
python -W error "$p" "$a" --joined --duration 2e-6 --steps 6 --order 8 --attribute-feedback --stable-pole-tangent --outer 1.00001e-5 5
python -W error "$p" "$a" --joined --duration 2e-6 --steps 6 --order 8 --attribute-feedback --stable-pole-tangent --outer 1.00001e-5 5 --warm-multiplier 16
python -W error "$p" "$a" --joined --duration 2e-6 --steps 12 --order 8 --outer 1.00001e-5 5
for n in 6 12; do
  python -W error "$p" "$a" --joined --duration 2e-6 --steps "$n" --order 8 --projection cell-conservative
done
for m in 0 16 64; do
  python -W error "$p" "$a" --joined --duration 2e-6 --steps 6 --order 8 --projection cell-conservative --warm-multiplier "$m"
done
python -W error "$p" "$a" --joined --duration 2e-6 --steps 12 --order 8 --projection cell-conservative --elastic-energy-budget 5e-15
```

Raw: `results/early-pole-paired-warm-feedback-attribution-2026-11-27.json.gz`. No inferred abundance or general error bound follows.
