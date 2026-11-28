# A warm tail evolved with its comoving energy exchange

From the same accepted final fine n216 snapshot, compare two short joined continuations with all dielectric tables held fixed: a fixed log scale and a prescribed path log(a)=t. Both have six macro steps and t=2e−6, and both use the [same comoving reciprocal decay packets](INDEPENDENT-EARLY-POLE-WARM-EXTERIOR-REDSHIFT-SNAPSHOT-2026-11-27.md), hard photon pairs, Compton, flavor conversion and the two screened elastic channels. Use the fixed radiation-clock normalization of the reference, not a solved Friedmann clock. Each path has a control with the same selected warm distribution frozen and a second arm where it evolves and exchanges number and energy with the hard states. The cold support is unchanged.

| Warm scalar representation | Warm gain at fixed scale | Warm gain with prescribed redshift |
| --- | ---: | ---: |
| q=.1–5, native markers | 5.4968731103e−13 | 5.4968896103e−13 |
| Same eleven shells, distributed q4 | 5.2291217119e−13 | 5.2291374082e−13 |

Redshift changes the warm gain by about +3.002 ppm in both representations in this very short window. Resolving the warm quadrature changes the absolute warm gain by about −4.87% at either clock. The evolving runs independently integrate hard redshift work (about 7.7856e−10); subtracting it from the total comoving-energy increment leaves a defect below 1e−15 in every arm. The terminal hard-number and entropy checks also pass. A difference between independently integrated cold populations is not used to measure tiny warm-to-cold feedback. The prescribed scale comparison does not determine a physical expansion rate, evolving polarization, or a late abundance. One macro width is not a time convergence test.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
a=results/early-pole-joined-fine-216-memory-2026-11-27.json.gz
p=early_pole_live_scalar_exterior_frw_reference.py
python -W error "$p" "$a"
python -W error "$p" "$a" --expands
python -W error "$p" "$a" --warm-order 4
python -W error "$p" "$a" --warm-order 4 --expands
```

Raw: `results/early-pole-live-warm-redshift-2026-11-27.json`.
