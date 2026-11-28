# One extra warm node recovers most of the actual G9 startup pressure change

The fixed trial chart upgrades one original cell [3.96368,5.72884] from G8 to G9 and leaves the other cells unchanged (300 total scalar nodes). It independently evolves from the same empty-daughter startup, with the 2048/4096 pair reaching x=.0075. The paired scalar-yield drift is 5.116e-5; the maximum absolute continuity defect is 5.64e-16 and maximum paired-number residual 1.74e-16. The chart identity is bound to the exact earlier held-source control, not an adaptively chosen later state.

The direct three-way comparison uses the same native S4 measurement mass. Relative to G8, pressure rises by 0.0273012% for uniform G9 (321 nodes) and 0.0231766% for mixed300. Thus one additional node accounts for about 84.9% of the 22-node actual G8→G9 pressure change at this clock. The effect matches the independently measured dominant native cell at startup, but this is not evidence of improved continuum accuracy or performance at later clocks.

```sh
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
d=2026-11-28; p=results/early-pole-fine-charge2-scalar4
python -W error -m unittest -q test_early_pole_endpoint_warm_single_upgrade_history
python -W error early_pole_endpoint_warm_300_held_clock_contrast.py \
 --control "results/early-pole-warm-single-g9-upgrade-$d.json" \
 --native-chain "$p-chain-$d.json" --native "$p-start-state-$d.json.gz" \
 --g8-reports "$p-gauss8warm-start-$d.json.gz" --g8-states "$p-gauss8warm-start-state-$d.json.gz" \
 --g9-reports "$p-gauss9warm-start-$d.json.gz" --g9-states "$p-gauss9warm-start-state-$d.json.gz" \
 --mixed-reports "$p-mixed300warm-start-$d.json.gz" --mixed-states "$p-mixed300warm-start-state-$d.json.gz"
```
