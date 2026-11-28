# Both S8/S16 explicit scalar replacement orders on the same original display

On the original S8 state, directly evaluate the two children of each of the
58 retained cold q parents crossing an S16 scalar edge. Keep the native
S16 calculation unchanged. Both complete joined S8 masks have 158 panels.
No child rate is inherited. Both replacement directions compare positive
panel sums on the original 17-bin display.

| Observable | Native TV | Explicit population TV | Combined kernel TV | Explicit share of integrated log change |
|---|---:|---:|---:|---:|
| forward number | 0.1108233% | 0.1108232% | 7.845988e-07% | 98.31899% |
| inverse number | 0.1107019% | 0.1107018% | 1.525291e-06% | 98.72704% |
| forward excess | 0.1085988% | 0.1085986% | 1.06184e-06% | 99.75427% |
| inverse excess | 0.1087766% | 0.1087769% | 1.915355e-06% | 99.79093% |

The explicit and combined-kernel terms add as **signed** profiles,
not as TV magnitudes. In the forward-number channel the S16/S8 integrated
relative change is about −0.02057%; the fixed-population TV stays below
7.9e−9 as a fraction in either substitution order. The four conditional
scalar-energy versus number neutrality separations lie between 3.350e−8
and 3.374e−8. No physical photon occupation is assigned.

The two meshes reproduce the four-channel total rectangle within
2.03e−13 relative and normalized profiles within 2.08e−14 absolute.
The largest native-bin relative recomposition error across all four
channels is below 1.77e−13 on the reported mesh. The combined kernel has
not been further decomposed into optical and charged contributions.

```sh
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
d=2026-11-28; p=results/fine-charge2-scalar; s=results/early-pole-fine-charge2-scalar
python -W error early_pole_endpoint_cold_spatial_population_exact_rectangle.py --levels 8 16 \
 --left-baseline "$p"8-endpoint-continuum-cold-{1,2}.jsonl.gz \
 --left-refined "$p"8-s16-population-edge-endpoint-continuum-cold-{1,2}.jsonl.gz \
 --right "$p"16-endpoint-continuum-cold-{1,2}.jsonl.gz --states "$s"{8,16}-x05-state-$d.json.gz \
 --archive results/finite-mass-qed-fixed-pole-cold-scalar-band-$d.json \
 --plan "results/early-pole-s8-s16-population-edge-split-plan-$d.json" \
 --display-fixture "results/early-pole-cold-spatial-population-factor-bounds-$d.json"
```
