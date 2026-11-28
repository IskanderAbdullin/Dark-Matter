# Actual startup localizes the mixed300 pressure effect to the selected cell

Use the three exact independently evolved startup states and the same native measurement mass. From their resolved per-original-cell count and pressure moments, form a diagnostic reference with G8 moments everywhere except the selected [3.96368,5.72884] cell, which uses G9. It is a numerical counterfactual, not a separately integrated mixed trajectory.

The observed mixed300 total pressure increase relative to G8 is 2.31766e-4. The selected cell differs from its uniform G9 counterpart by only 1.42e-12 of total G8 pressure. Across all other cells, the sum of absolute pressure deviations from G8 is 6.45e-11 of total G8 pressure, or about 2.78e-7 of the absolute mixed300 pressure change. Charged occupations do change, so this is near-localization of a pressure diagnostic at startup, not absence of coupled feedback or a universal locality theorem.

```sh
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
d=2026-11-28; p=results/early-pole-fine-charge2-scalar4
python -W error early_pole_endpoint_warm_evolved_single_cell_closure.py \
 --control "results/early-pole-warm-single-g9-upgrade-$d.json" \
 --accepted "results/early-pole-warm-mixed300-actual-startup-$d.json" \
 --states "$p"-{gauss8,gauss9,mixed300}warm-start-state-$d.json.gz
```

The comparison is limited to the held measurement at the first archived clock of these evolved histories; it does not determine a continuum result.
