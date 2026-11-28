# Direct comparison of the complete G8 and mixed255 warm histories

Recompute all particle count and pressure measurements from the five actual mixed255 and G8 snapshots, using exactly the same accepted original S4 mass and paired chronology. Verify the exact saved restart traces and independently close the ratios through the two earlier G7 comparisons. Closure is within 2e-14 across all measurements.

The 255-node mixed chart has 0.03313% to 0.03370% higher pressure than the 299-node G8 chart across the five measured clocks; its count differs by +8.68e-10 in relative units at x=.05. The same native cell q=[3.96367894,5.72884111] accounts for 91.15% of the terminal pressure difference. Redistribution among cells adds +0.0345945 percentage points, partly offset by -0.000896533 points within cells. These finite matched histories do not identify which chart is closer to a continuum or physical answer. The reported pressure envelope is conditional on original-cell particle and energy moments. It does not estimate solver time-step error, which requires comparing nodal states recovered at both time steps.

```sh
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
d=2026-11-28; p=results/early-pole-fine-charge2-scalar4
python -W error early_pole_endpoint_warm_three_chart_closure.py \
 --control "results/early-pole-warm-budgeted-gauss-probe-$d.json" \
 --native-chain "$p-chain-$d.json" \
 --mixed-g7-lock "results/early-pole-warm-mixed255-full-prefix-$d.json" \
 --g8-g7-lock "results/early-pole-warm-gauss8-actual-prefix-x05-$d.json" \
 --native-states "$p"-{start,x01,x02,x04,x05}-state-$d.json.gz \
 --mixed-reports "$p"-mixed255warm-{start,x01,x02,x04,x05}-$d.json.gz \
 --mixed-states "$p"-mixed255warm-{start,x01,x02,x04,x05}-state-$d.json.gz \
 --g8-reports "$p"-gauss8warm-{start,x01,x02,x04,x05}-$d.json.gz \
 --g8-states "$p"-gauss8warm-{start,x01,x02,x04,x05}-state-$d.json.gz
```
