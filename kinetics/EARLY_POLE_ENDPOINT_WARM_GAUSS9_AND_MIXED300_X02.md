# Actual G9 and mixed300 through x=.02, with exact paired-clock pressure

Both independently evolved charts continue the exact accepted x=.01 states on 1024/2048 steps. G9 has scalar-yield drift 2.185e-6, maximum absolute continuity defect 1.69e-13, and paired-number residual below 9.33e-17. At the common original S4 measurement mass, pressure relative to G8 is higher by 0.0272498% for G9 and 0.0230732% for mixed300. The single extra node retains 84.67% of the uniform G9 pressure change at this clock.

In the evolved cell comparison, the summed absolute off-target pressure residual grows to 1.70e-9 of total G8 pressure, or 7.38e-6 of the mixed300 absolute cellwise pressure change. The corresponding off-target number residual accounts for about 8.8% of the very small cellwise number change. Pressure localization thus does not imply absence of count or charged-distribution feedback.

The complete 1024-step checkpoint from each chart binds its exact clock, restart, geometry and recorded original coarse trajectory. Their within-chart coarse-minus-fine pressure changes are -1.2475e-7 and -1.2481e-7 in fractional units. The actual mixed300/G9 pressure contrast changes only 6.00e-11 between the matched 1024 and 2048 clocks. This is time-resolution evidence for this segment, not a full-history global time-error estimate or a spatial convergence result.

```sh
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
d=2026-11-28; p=results/early-pole-fine-charge2-scalar4; ctr="results/early-pole-warm-single-g9-upgrade-$d.json"
python -W error early_pole_endpoint_warm_300_held_clock_contrast.py \
 --control "$ctr" --native-chain "$p-chain-$d.json" --native "$p"-{start,x01,x02}-state-$d.json.gz \
 --g8-reports "$p"-gauss8warm-{start,x01,x02}-$d.json.gz --g8-states "$p"-gauss8warm-{start,x01,x02}-state-$d.json.gz \
 --g9-reports "$p"-gauss9warm-{start,x01,x02}-$d.json.gz --g9-states "$p"-gauss9warm-{start,x01,x02}-state-$d.json.gz \
 --mixed-reports "$p"-mixed300warm-{start,x01,x02}-$d.json.gz --mixed-states "$p"-mixed300warm-{start,x01,x02}-state-$d.json.gz
python -W error early_pole_endpoint_warm_checkpoint_clock_comparison.py \
 --control "$ctr" --accepted "results/early-pole-warm-mixed300-actual-prefix-x02-$d.json" \
 --g9 "$p-gauss9warm-x01-state-$d.json.gz" "$p-gauss9warm-x02-$d.json.gz" "$p-gauss9warm-x02-state-$d.json.gz" \
  "results/early-pole-warm-gauss9-x02-1024-checkpoint-$d.json.gz" "results/early-pole-warm-gauss9-x02-provenance-$d.sha256" \
 --mixed "$p-mixed300warm-x01-state-$d.json.gz" "$p-mixed300warm-x02-$d.json.gz" "$p-mixed300warm-x02-state-$d.json.gz" \
  "results/early-pole-warm-mixed300-x02-1024-checkpoint-$d.json.gz" "results/early-pole-warm-mixed300-x02-provenance-$d.sha256"
```

This is an actual numerical prefix, not a continuum or physical outcome.
