# Both S4/S8 explicit scalar replacement orders on the same original display

On the original S4 state, directly re-integrate the two children of each of the 29 cold q panels crossing an S8 scalar edge. Keep the native S8 calculation unchanged. The 100 joined S4 panels pass both complete band gates (mesh discrepancy 1.49e-13). The largest nine-coefficient original panel recomposition defects are 2.76e-13 and 1.75e-13 on the two meshes. No child rate was inherited from a parent. Both replacement paths then use known positive panel sums; the normalized display uses the same 17 q bins as the earlier spatial ladder.

| Observable | Native TV | Population contrast TV, averaged as signed profiles | Bare-kernel contrast TV, averaged as signed profiles | Explicit population share of integrated log change |
| --- | ---: | ---: | ---: | ---: |
| forward number | 0.4516501% | 0.4516477% | 1.267184e-05% | 93.56826% |
| inverse number | 0.4511404% | 0.4511382% | 2.222951e-05% | 94.40285% |
| forward excess | 0.4368845% | 0.436882% | 1.719227e-05% | 99.01708% |
| inverse excess | 0.4376052% | 0.437608% | 2.81216e-05% | 99.04221% |

The counterfactual kernel-only forward TV is 1.281e-07 at the S4 occupation and 1.253e-07 at the S8 occupation. The symmetric terms add as **signed** normalized profile vectors; their TV magnitudes need not add. The whole-band amplitude columns use averaged log contrasts, not an attribution of a signed physical flux.

The relative scalar-energy vs number neutrality split on the four crossings spans 3.374e-08 to 3.472e-08; this is a conditional diagonal and counterfactual ledger, with no physical photon occupation assigned. Two-mesh gross and profile defects are 1.51e-13 relative and 2.39e-14 absolute. No separate decomposition of the bare optical/charge kernel is inferred.

To recompute from the published complete raw bands and split plan:

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
d=2026-11-28; p=results/fine-charge2-scalar; s=results/early-pole-fine-charge2-scalar
python -W error early_pole_endpoint_cold_spatial_population_exact_rectangle.py --levels 4 8 \
 --left-baseline "$p"4-endpoint-continuum-cold-{1,2}.jsonl.gz \
 --left-refined "$p"4-s8-population-edge-endpoint-continuum-cold-{1,2}.jsonl.gz \
 --right "$p"8-endpoint-continuum-cold-{1,2}.jsonl.gz --states "$s"{4,8}-x05-state-$d.json.gz \
 --archive results/finite-mass-qed-fixed-pole-cold-scalar-band-$d.json \
 --plan results/early-pole-s4-s8-population-edge-split-plan-$d.json \
 --display-fixture results/early-pole-cold-spatial-population-factor-bounds-$d.json
```
