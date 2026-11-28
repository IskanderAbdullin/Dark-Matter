# Native fine/third D2 cold comparison without interpolation

Both independently evolved fine and third D2 endpoints use 2048 final steps. Compare the two complete native held cold bands with 49 and 81 panels per mesh, respectively, over the same retained support. The chart, scalar gap, and fitted optical pole change together. The direct measurement is not a one-factor test.

At every union edge, bound the cumulative response between the sum of complete earlier panels and that sum plus any panel crossing the cut. This does not assume a within-panel interpolation. The exact shared panel edges alone are 1e-06, 1e-05, 1.4427e-05, 1.5e-05. A gap between the two cumulative intervals is a lower bound on normalized TV, not its exact value.

| Response | Third / fine retained gross response minus one | Lower bound from one cumulative cut | Minimum TV compatible with both native ledgers |
| --- | ---: | ---: | ---: |
| forward number | -83.07944% | 17.85102% | 26.23343% |
| inverse number | -85.8182% | 17.88863% | 26.40743% |
| forward excess | -84.54903% | 19.66131% | 31.49369% |

The sharper minimum uses maximum common mass on intersections of the native intervals. Its feasibility witness assigns nonnegative shared mass to union subintervals with both panel marginals as caps; independent nonnegative panel prices cap the overlap above. The compatible minimum is one minus maximum overlap, not the true profile TV: no distribution within a panel is presumed. For densities on positive-length intervals, the extremal assignments can be realized without atoms.

The maximum difference between the two independent q assessments is 3.18e-14 in CDF enclosure endpoints and 1.35e-13 relative in band ratio. This measures finite held response differences; no physical whole-spectrum or spatial continuum inference is made.

From `kinetics/`:

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
d=2026-11-28; p=results
python -W error early_pole_endpoint_cold_non_nested_cdf_contrast.py \
 --fine "$p/fine-d2-endpoint-continuum-cold"-{1,2}.jsonl.gz \
 --third "$p/third-d2-endpoint-continuum-cold"-{1,2}.jsonl.gz \
 --states "$p/early-pole-cell-endpoint"-{fine,third}-d2-x05-state-$d.json.gz \
 --archive "$p/finite-mass-qed-fixed-pole-cold-scalar-band-$d.json"
```
