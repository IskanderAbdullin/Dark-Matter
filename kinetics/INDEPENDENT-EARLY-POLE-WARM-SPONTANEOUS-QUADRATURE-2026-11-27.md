# The warm source difference is in the spontaneous quadrature

At the same final fine n216 hard state and exactly the same original scalar pole, evaluate the eleven native scalar shells selected by marker q∈[.1,5]. Compare one marker per shell with Gauss nodes in q³ on those same shell volumes. This selection is by marker, not a clip at q=5. For each quadrature node the signed count source is (S+f_phi S-f_phi I), where (S=f_P(1-f_d)) and (I=f_d(1-f_P)) enter the same positive packet measure. Evaluating at (f_phi=0) and (1) in the frozen rate kernel separates the terms; it does not change the scalar state or the pole.

| Hard reconstruction, outer order 12 | Marker minus resolved net | Spontaneous part | Occupation-dependent part |
| --- | ---: | ---: | ---: |
| Native steps | 1.33848925e−8 | 1.33850378e−8 | −1.45299e−13 |
| Continuous | 1.34394670e−8 | 1.34396125e−8 | −1.45414e−13 |

The net marker-to-resolved change is approximately −4.87%; the signed occupation-dependent change cancels roughly 11 parts per million of that gap. On **identical resolved nodes**, changing from constant shell occupation to the positive energy-matched occupation changes the net rate by only about 1.22 ppm of the rate itself. The warm inverse rate is about 2.9 ppm of the full resolved net. Changing the inner radial order from 8 to 16 leaves the gap unchanged at the precision quoted. Outer orders 4, 8 and 12 move the overall resolved conclusion by less than .01% of the rate; this does not establish a fully converged continuum rate.

The [earlier exterior continuations](INDEPENDENT-EARLY-POLE-EXTERIOR-JOINED-AND-DISTRIBUTED-2026-11-27.md) tested distributed warm shells and live hard channels in separate runs. Combining both on the same native step snapshot (2e−6 exposure, six joined macro steps) gives warm gains (5.2291145752	imes10^{-13}) and (5.2287379651	imes10^{-13}) at shell orders 4 and 8. Those are 4.87–4.88% below the joined marker gain. Adding the harder channels increases each distributed-shell gain by about .0046% relative to its decay-only counterpart. Terminal energy closure is at most (4.45	imes10^{-16}) in magnitude across these combined arms. Differences between separately advanced cold populations at the (10^{-18}) level are not resolved by this test. The results are for the same finite hard-box model, not an inferred abundance or a measurement.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
a=results/early-pole-joined-fine-216-memory-2026-11-27.json.gz
for p in steps cell-conservative; do
  python -W error early_pole_warm_seed_measure_reference.py "$a" --projection "$p"
done
for n in 4 8; do
  python -W error early_pole_live_scalar_exterior_decay_reference.py "$a" --joined --warm-order "$n" --panels 6
done
```

Raw: `results/early-pole-warm-spontaneous-quadrature-2026-11-27.json.gz`.
