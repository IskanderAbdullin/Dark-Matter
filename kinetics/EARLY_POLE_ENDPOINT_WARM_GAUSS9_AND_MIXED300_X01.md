# Actual G9 and mixed300 continuations retain the single-cell effect to x=.01

Advance the actual accepted G9 and mixed300 startup states through an exact 512/1024-step continuation. Both validate their exact previous state hashes and report scalar-yield clock drift below 1.253e-7, absolute continuity defects below 1.07e-16 and maximum paired-number residuals below 3.10e-17.

At the common original S4 mass, relative G9 pressure is 0.0273012% higher than G8; mixed300 is 0.0231884% higher. The one-extra-node rule therefore retains about 84.9% of the G9 pressure change after this continuation, with relative count changes only +3.91e-9 and +3.51e-9 respectively. Independently transplanting only the selected G9 cell in measured moments accounts for the change within the small counterfactual residual. At x=.01 the sum of the absolute off-target pressure differences is 1.90e-10 of G8 total pressure. It is about 8.20e-7 of the total mixed300 absolute cellwise pressure change. This is sampled localization, not an absence of dynamical feedback.

```sh
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
d=2026-11-28; p=results/early-pole-fine-charge2-scalar4; ctr="results/early-pole-warm-single-g9-upgrade-$d.json"
python -W error early_pole_endpoint_warm_g7_prefix_check.py --levels 8 9 \
 --native "$p"-{start,x01}-state-$d.json.gz --native-chain "$p-chain-$d.json" \
 --reports6 "$p"-gauss8warm-{start,x01}-$d.json.gz --states6 "$p"-gauss8warm-{start,x01}-state-$d.json.gz \
 --reports7 "$p"-gauss9warm-{start,x01}-$d.json.gz --states7 "$p"-gauss9warm-{start,x01}-state-$d.json.gz
python -W error early_pole_endpoint_warm_300_held_clock_contrast.py \
 --control "$ctr" --native-chain "$p-chain-$d.json" --native "$p"-{start,x01}-state-$d.json.gz \
 --g8-reports "$p"-gauss8warm-{start,x01}-$d.json.gz --g8-states "$p"-gauss8warm-{start,x01}-state-$d.json.gz \
 --g9-reports "$p"-gauss9warm-{start,x01}-$d.json.gz --g9-states "$p"-gauss9warm-{start,x01}-state-$d.json.gz \
 --mixed-reports "$p"-mixed300warm-{start,x01}-$d.json.gz --mixed-states "$p"-mixed300warm-{start,x01}-state-$d.json.gz
python -W error early_pole_endpoint_warm_evolved_single_cell_closure.py --control "$ctr" \
 --accepted "results/early-pole-warm-mixed300-actual-prefix-x01-$d.json" \
 --states "$p"-{gauss8,gauss9,mixed300}warm-x01-state-$d.json.gz
```

This is an actual early numerical prefix only, not a statement of convergence or the final history.
