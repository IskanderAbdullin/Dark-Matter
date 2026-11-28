# Joint coarse-clock perturbations introduce a small, resolved interaction at two cold points

On each of S2 and S4, compare all four combinations of 2048 and 1024 steps for x=.02–.04 and x=.04–.05. Three corners are already independently archived. For the fourth, start the final 1024-step integration from the actual earlier-coarse .04 state; the digest identifies that restart and the saved terminal state. All corners retain the same spatial chart. The two q-point coefficients are recalculated at their separately fitted optical gaps and poles with the same two quadratures.

For forward number, the results in fractional **log** units are:

| Scalar chart | q | Earlier clock alone | Final clock alone | Joint interaction |
| --- | ---: | ---: | ---: | ---: |
| S2 | 3.5e−6 | −0.00461850 | −0.00283432 | −1.79532e−5 |
| S2 | 4.7e−6 | −0.00453631 | −0.00278699 | −1.80986e−5 |
| S4 | 3.5e−6 | −0.00455615 | −0.00285174 | +1.06177e−7 |
| S4 | 4.7e−6 | −0.00447809 | −0.00280238 | +3.27578e−6 |

The combined log shift is the sum of these three columns. On S2 the interaction is about 0.24% of the combined point log shift; on S4 it is smaller and changes sign relative to S2. The gross inverse and affine columns are retained separately. The worst point quadrature discrepancy is 8.38e−13 relative, and the earlier-only point edge independently reproduces the archived comparison. This resolves a small non-additivity on the tested states and points; it is not a whole-band proof of temporal convergence.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
d=2026-11-28; p=results/early-pole-fine-charge2-scalar; s=results/early-pole-cold
for n in 2 4; do
  last="$p$n-coarse-x05-state-$d.json.gz"
  [ "$n" != 4 ] || last="results/early-pole-geometric-warm-w1-coarse-x05-state-$d.json.gz"
  b="$s-s$n-prior-clock"
  python -W error early_pole_endpoint_cold_two_segment_clock_rectangle.py --level "$n" \
    --states "$p$n-x05-state-$d.json.gz" "$last" "$b-x05-state-$d.json.gz" "$b-x04-state-$d.json.gz" \
    --history "$b-report-$d.json" --points "$s-prior-clock-points-$d.json" \
    --final-fixture "$s-scalar-band-clock-contrast-$d.json" --combined "/tmp/S$n-both-coarse.json"
done
```
