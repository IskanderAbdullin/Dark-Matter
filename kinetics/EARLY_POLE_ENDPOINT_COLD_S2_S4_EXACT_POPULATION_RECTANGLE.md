# Both S2/S4 explicit scalar replacement orders on the same original display

On the original S2 state, directly re-integrate the two children of each of the 14 cold q panels crossing an S4 scalar edge. Keep the native S4 calculation unchanged. The 71 joined S2 panels pass both complete band gates (mesh discrepancy 7.62e-14). The largest nine-coefficient original panel recomposition defects are 4.24e-13 and 1.46e-13 on the two meshes. No child rate was inherited from a parent. Both replacement paths then use known positive panel sums; the normalized display uses the same 17 q bins as the earlier spatial ladder.

| Observable | Native TV | Population contrast TV, averaged as signed profiles | Bare-kernel contrast TV, averaged as signed profiles | Explicit population share of integrated log change |
| --- | ---: | ---: | ---: | ---: |
| forward number | 2.02402% | 2.024015% | 2.010987e-05% | 97.50343% |
| inverse number | 2.022563% | 2.02256% | 3.519135e-05% | 98.092% |
| forward excess | 1.907637% | 1.907623% | 2.757828e-05% | 99.6284% |
| inverse excess | 1.890783% | 1.890787% | 4.504531e-05% | 99.68145% |

The counterfactual kernel-only forward TV is 2.1e-07 at the S2 occupation and 1.922e-07 at the S4 occupation. The symmetric terms add as **signed** normalized profile vectors; their TV magnitudes need not add. The whole-band amplitude columns use averaged log contrasts, not an attribution of a signed physical flux.

The relative scalar-energy vs number neutrality split on the four crossings spans 3.472e-08 to 3.866e-08; this is a conditional diagonal and counterfactual ledger, with no physical photon occupation assigned. Two-mesh gross and profile defects are 1.58e-13 relative and 4.08e-14 absolute. No separate decomposition of the bare optical/charge kernel is inferred.

To recompute from the published complete raw bands and split plan:

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
d=2026-11-28; p=results/fine-charge2-scalar; s=results/early-pole-fine-charge2-scalar
python -W error early_pole_endpoint_cold_spatial_population_exact_rectangle.py --levels 2 4 \
 --left-baseline "$p"2-endpoint-continuum-cold-{1,2}.jsonl.gz \
 --left-refined "$p"2-s4-population-edge-endpoint-continuum-cold-{1,2}.jsonl.gz \
 --right "$p"4-endpoint-continuum-cold-{1,2}.jsonl.gz --states "$s"{2,4}-x05-state-$d.json.gz \
 --archive results/finite-mass-qed-fixed-pole-cold-scalar-band-$d.json \
 --plan results/early-pole-s2-s4-population-edge-split-plan-$d.json \
 --display-fixture results/early-pole-cold-spatial-population-factor-bounds-$d.json
```
