# A mixed warm order selected at x=.01 improves the held x=.05 source at a fixed node budget

At the stored S4 x=.01 state, use the previously evaluated G4–G8 cell rates and minimize the sum of absolute local net residuals relative to Volume64. Keep a fixed charge grid and original cold nodes. An integer dynamic program selects a single positive Gaussian order per warm cell at each maximum scalar-node budget. Exhaustive small tests verify the optimizer. The x=.05 source is not used to choose either rule.

Independently re-evaluate each chosen mixed chart directly on both radial meshes at both x=.01 and .05. The calculation preserves the original cell counts and native edges. Direct evaluations agree with independent per-cell assembly within 1.2e−16 of the full normalized source. At the two recorded frozen states:

| Direct rule | Scalar nodes | Local L1 at x=.01 | Local L1 at x=.05 |
| --- | ---: | ---: | ---: |
| Uniform G6 | 255 | 0.0422235% | 0.0529980% |
| Uniform G7 | 277 | 0.0263033% | 0.0342180% |
| Uniform G8 | 299 | 0.0244703% | 0.0286017% |
| Mixed 219 | 219 | 0.0213752% | 0.0316515% |
| Mixed 255 | 255 | 0.0205348% | 0.0258389% |

The mixed 255-node rule beats uniform G8 on this finite reference measure while using 44 fewer nodes. The mixed 219-node rule improves over G6 and G7 here, but not over G8 at x=.05. This tests transfer of a frozen selection between two snapshots on the same underlying trajectory, not an independent sample of production histories or evolving stability. Volume64 is a compared numerical reference, not an exact continuum source.

```sh
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
d=2026-11-28; p=results
python -W error early_pole_endpoint_warm_budgeted_gauss_probe.py \
 --states "$p/early-pole-fine-charge2-scalar4-x"{01,05}-state-$d.json.gz \
 --reference "$p/early-pole-endpoint-warm-source-cells-x"{01,05}-$d.json \
 --comparisons "$p/early-pole-warm-frozen-gauss8-source-x"{01,05}-$d.json
```
