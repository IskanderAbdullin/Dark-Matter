# Higher volume orders verify which diagnostic warm cells benefit from splitting

Keep the archived S4 occupations, all charged nodes, fitted mass and collision clock unchanged. Integrate four diagnostic warm cells, including the two proposed refinement cells and two failure controls, with local positive volume Gauss–Legendre orders 32, 64, 128 and 256. All other cells stay unchanged. Selective integration reproduces the published order-32 and order-64 cell rates; each new volume order is paired with both internal momentum quadratures. Order 128→256 changes the largest selected net contribution by at most 7.20e-7 of the original total net warm source across all four snapshots. This is a convergence comparison, not a strict bound.

Against the order-256 results, splitting reduces error in both proposed cells at every snapshot. In total-source fractions, its improvement across the four snapshots ranges 1.80e-4 to 2.48e-4 for q=[3.96368,5.72884], and 3.08e-4 to 3.37e-4 for q=[5.72884,8.28009]. In the next higher cell q=[8.28009,11.96750], the very same split instead increases error by 1.74e-5 to 1.94e-5. The lower control q=[0.90829,1.31278] improves at x=.04 but not at .01, .02 or .05. At x=.04 its unsplit error is 8.90e-5, still larger than the split-discrepancy estimate 6.96e-5. These controls support selective refinement rather than a universal split stopping criterion. None measures a new trajectory.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
p=results/early-pole-fine-charge2-scalar4; w=results/early-pole-endpoint-warm-source-cells; d=2026-11-28
for i in 01 02 04 05; do
 python -W error early_pole_endpoint_warm_selected_volume_audit.py "$p-x$i-state-$d.json.gz" "$w-x$i-$d.json"
done
```
