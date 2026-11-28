# Native S8/S16 shared reference: independent three-field loss

Using a single S8 optical mass for the reference mean excess and fixed original S8 scales for the three response fields, the native S8/S16 charts have joint shared-count-only minimum L1 between 4.11400e−6 and 4.11469e−6 over five reference blends. With count and mean both fixed, on all 60 common parents the pure S8 reference gives 8.79706e−6; the four nonzero S16 blends are infeasible on those parents before S8 reintegration. The maximal shared contiguous fallback has 36 conditioning parents but retains the 60 display bins; its joint count-and-mean loss over all five blends is 8.29473e−6 to 8.30064e−6.

The unscaled LP failed the independent bracket on at least one coarsened case. Scaling the response-space LP proposal (by 256 for these cards) does not change the optimization's feasible profiles or its unscaled objective. Both support and witness are evaluated on the original coefficients. Their width is no more than 5.31e−12 on the reported coarsened results; the optimum differences across q meshes are at most 3.53e−13. The independent acceptance width is unchanged. These are numerical bracketing checks of the native held panel response, not a prediction for the physical reaction or the refined S8 panels.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
d=2026-11-28; p=results/fine-charge2-scalar; s=results/early-pole-fine-charge2-scalar
arc=results/finite-mass-qed-fixed-pole-cold-scalar-band-$d.json
python -W error early_pole_endpoint_cold_shared_reference_scan.py --levels 8 16 \
 --left "$p"8-endpoint-continuum-cold-{1,2}.jsonl.gz --right "$p"16-endpoint-continuum-cold-{1,2}.jsonl.gz \
 --states "$s"{8,16}-x05-state-$d.json.gz --archive "$arc"
```
