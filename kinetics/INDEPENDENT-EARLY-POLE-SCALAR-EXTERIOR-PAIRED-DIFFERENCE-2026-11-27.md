# Resolve the small feedback of live exterior decay on the cold count

The two separate adaptive trajectories in the [live exterior test](INDEPENDENT-EARLY-POLE-LIVE-SCALAR-EXTERIOR-DECAY-2026-11-27.md) resolve the warm gain, but their much smaller difference in cold gain is below the numerical variation of either trajectory. Instead evolve the common cold-only base and a dedicated difference variable on one RK4 grid. Its derivative is the difference of the cold-only operators at the two reconstructed states, plus the additional reciprocal exterior operator. Track warm gains from zero separately. Both arms begin from the identical final accepted five-channel hard snapshot and unchanged scalar partition; the continuation itself still uses decay and inverse decay only.

| Fine reconstruction, q=.1–5 promotion | Cold-gain difference, 24 common steps | Cold-gain difference, 48 common steps | Warm gain (48 steps) |
| --- | ---: | ---: | ---: |
| Native steps | −9.806e−16 | −9.938e−16 | 1.978738e−11 |
| Count-conserving continuous | −1.007e−15 | −1.001e−15 | 1.987930e−11 |

Exposure is 7.2e−5. The differential correction is about −1.2e−10 relative to the cold gain in either setup, and is consistent with the separately measured initial directional derivative. The two common-grid resolutions differ by less than 1.4% on this tiny correction, even though the absolute cold integral error judged by comparing independently solved arms is larger than the correction. Both reported absolute energy residuals are below 3e−16. As an additional finite-support check, moving the wider q>1e−5 exterior in the native state gives cold feedback −2.980e−17 and warm gain +8.779e−13 over a much shorter 2e−6 comparison, again consistent with the initial linear probe. This does not determine feedback for an expanding or multi-channel history or the final warm population.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1
b=results/early-pole-joined-fine-216-memory-2026-11-27.json.gz
for p in steps cell-conservative; do
 for n in 24 48; do
  python -W error early_pole_live_scalar_exterior_paired_difference_reference.py "$b" \
   --projection "$p" --steps "$n"
 done
done
python -W error early_pole_live_scalar_exterior_paired_difference_reference.py "$b" \
 --duration 2e-6 --steps 12 --outer 1.00001e-5 5
```
Raw: `results/early-pole-live-scalar-exterior-paired-difference-2026-11-27.json`.
