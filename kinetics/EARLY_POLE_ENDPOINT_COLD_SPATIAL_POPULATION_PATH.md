# One exact replacement path separates the cold spatial change

In each adjacent spatial pair, start at the original coarse scalar grid. First evaluate the fine-grid bare kernel with the coarse grid's explicit scalar population, then switch that explicit population on the fine-grid kernel. These replacements concern (1+f) in forward number and forward excess, and (f) in inverse number. They keep all other differences (scalar gap, pole and charged state included) together in the bare response. Every fine q panel carries only one coarse scalar occupation, so both transitions use the stored quadratures directly without subpanel interpolation.

On the same 17-bin display, retaining the coarse explicit population makes the normalized difference between the fine and coarse bare kernels very small:

| Pair | Native forward-number TV | At fixed coarse explicit population | Explicit population share of total forward log change on this path |
| --- | ---: | ---: | ---: |
| S2/S4 | 0.02024020 | 2.10e−7 | 97.50% |
| S4/S8 | 0.00451650 | 1.28e−7 | 93.57% |
| S8/S16 | 0.00110823 | 7.87e−9 | 98.32% |

All three fields follow this pattern on this path; the S8/S16 integrated log shares for inverse number and forward excess are 98.73% and 99.75%, respectively. This is different from changing the **final time step on the same spatial chart**, for which the previously tested fixed-geometry factor split was dominated by the bare kernel. The two controls change different axes. The exact spatial path is not an order-independent attribution: some original coarse q panels straddle fine population boundaries (58 for S8/S16), and their non-interpolated positivity bounds leave the opposite path poorly determined. Refined coefficients are required to resolve it.

Both meshes agree to 2.02e−13 relative for the response bounds and 4.05e−14 absolute for the log/TV diagnostics. The displayed TV here uses actual response profiles after factor replacement, not the separately optimized minimum with moments held and not a physical spectrum.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
d=2026-11-28; p=results/fine-charge2-scalar; s=results/early-pole-fine-charge2-scalar
arc=results/finite-mass-qed-fixed-pole-cold-scalar-band-$d.json
python -W error early_pole_endpoint_cold_spatial_population_factor_bounds.py \
 --s2 "$p"2-endpoint-continuum-cold-{1,2}.jsonl.gz --s4 "$p"4-endpoint-continuum-cold-{1,2}.jsonl.gz \
 --s8 "$p"8-endpoint-continuum-cold-{1,2}.jsonl.gz --s16 "$p"16-endpoint-continuum-cold-{1,2}.jsonl.gz \
 --states "$s"{2,4,8,16}-x05-state-$d.json.gz --archive "$arc"
```

The accompanying population-edge plan adds the 58 missing **actual S16 scalar boundaries** to the S8 q grid, requiring both children to be re-integrated for each boundary in each mesh. This differs from the 41 geometry splits used for shared-moment feasibility. The plan is reproducible with `early_pole_endpoint_cold_population_edge_plan.py` using the same S8/S16 bands and states. Once those 116 new coefficients per mesh are available, the optional `--refined-s8` and `--s8-split-plan` arguments can check the opposite replacement path on the same original 17 display bins.
