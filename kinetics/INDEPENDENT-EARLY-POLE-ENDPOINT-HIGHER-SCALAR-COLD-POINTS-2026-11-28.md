# Higher scalar resolutions still move individual cold coefficients

Hold the fine charge d2 chart and x=.05 terminal clock fixed while independently varying the scalar cold subcells S4, S8 and S16. At k=.0015 and two q samples, switch each fitted scalar gap, photon pole, and charged continuum population through all eight corners; account separately and exactly for the scalar Bose/occupation factor. Retain both high and low radial-angle orders. This extends the [S1/S2/S4 point check](INDEPENDENT-EARLY-POLE-ENDPOINT-SCALAR-S4-COLD-BAND-2026-11-28.md) without filling in any higher-order q-band panels by interpolation.

| Transition | Gross forward at q=3.5e−6 | Gross forward at q=4.7e−6 |
| --- | ---: | ---: |
| S4→S8 | −5.9618% | −25.3174% |
| S8→S16 | +3.6309% | −14.0928% |

The occupation-only factors are 0.940434/0.746867 and 1.036313/0.859075, respectively. This nearly exhausts the observed point changes; the refitted gap/photon/charged factors are all much smaller. The lower-q changes reverse sign, and the upper-q difference remains sizable. None of these four differential comparisons establishes monotone or spatial convergence of an integrated band. All crossed numerical mesh differences are at most 4.99e−13.

Independently compare actual S4 endpoints with 1024 and 2048 final steps from the same archived .04 restart. The gross forward coarse/fine differences are −0.28477% and −0.27985%, respectively, and the inverse differences are −0.30979% and −0.30265%. These are an order of magnitude or more smaller than the last scalar-mesh movements at these points. The fine-time values reproduce the independently calculated S4 factorial weights to 1.22e−13. This rules out the last S4 time step as a sufficient explanation of the larger observed point shifts; it does not bound accumulated time error or errors on higher charts. Neither the pointwise nor the clock comparison specifies a physical net photon reaction.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
p=results/early-pole-fine-charge2-scalar; d=2026-11-28
python -W error early_pole_endpoint_scalar_band_point_factorial.py --higher \
  "$p"{4,8,16}-x05-state-$d.json.gz
python -W error early_pole_endpoint_s4_cold_clock_points.py \
  results/early-pole-geometric-warm-w1-coarse-x05-state-$d.json.gz "$p"4-x05-state-$d.json.gz \
  results/early-pole-geometric-warm-w1-pressure-clock-$d.json \
  results/early-pole-endpoint-scalar-band-point-factorial-$d.json
```
