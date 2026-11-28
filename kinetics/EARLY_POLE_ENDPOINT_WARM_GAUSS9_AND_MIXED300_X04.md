# Actual G9 and mixed300 through x=.04: pressure contrast and common-mode count drift

Starting at the exact accepted x=.02 occupation in each chart, both restarts reached x=.04 using 1024/2048 steps. Measured at the common native mass, G9 and the one-cell upgrade (mixed300) give pressure excesses of 0.0274928% and 0.0235095% relative to G8. Mixed300 retains 85.51% of the full G9 increment (84.67% at x=.02). This is a chart contrast at the same endpoint, not a continuum estimate.

The mixed300 off-target absolute pressure change is 1.34849e-8 of total G8 pressure and 5.73530e-5 of the mixed300 absolute pressure change, larger than at x=.02. The corresponding number fraction is 1.988%; this smaller fraction does not mean less absolute number feedback. The total absolute off-target number residual, normalized to G8, is 1.03216e-11. The moved cell spans q in [3.96368, 5.72884], and the largest remaining off-target pressure residual is in [1.89741, 2.74240]. The diagnostic substitutes the one cell in measurement only, not in a newly evolved trajectory.

The exact coarse checkpoints bind the initial states, nodal chart, step rule, interval, original recorded trajectory, and complete occupations. Within-chart coarse-minus-fine fractional pressure shifts are +3.01524e-7 (G9) and +2.95706e-7 (mixed300); the mixed300/G9 pressure contrast shifts only +5.81803e-9 going fine. Unlike pressure, the two individual coarse counts shift by about -1.27671e-4 relative to their fine counterparts; the count contrast itself shifts by only -3.79e-13. The latter is cancellation, not absolute convergence of either count.

Recalculate the audits from the archived restart and report bytes:

```bash
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
d=2026-11-28; p=results/early-pole-fine-charge2-scalar4; w=results/early-pole-warm
c="$w-single-g9-upgrade-$d.json"; a="$w-mixed300-actual-prefix-x04-$d.json"
native=("$p"-{start,x01,x02,x04}-state-$d.json.gz)
g8r=("$p"-gauss8warm-{start,x01,x02,x04}-$d.json.gz)
g8s=("$p"-gauss8warm-{start,x01,x02,x04}-state-$d.json.gz)
g9r=("$p"-gauss9warm-{start,x01,x02,x04}-$d.json.gz)
g9s=("$p"-gauss9warm-{start,x01,x02,x04}-state-$d.json.gz)
m3r=("$p"-mixed300warm-{start,x01,x02,x04}-$d.json.gz)
m3s=("$p"-mixed300warm-{start,x01,x02,x04}-state-$d.json.gz)
python -W error early_pole_endpoint_warm_g7_prefix_check.py --levels 8 9 \
 --native-chain "$p-chain-$d.json" --native "${native[@]}" \
 --reports6 "${g8r[@]}" --states6 "${g8s[@]}" --reports7 "${g9r[@]}" --states7 "${g9s[@]}"
python -W error early_pole_endpoint_warm_300_held_clock_contrast.py --control "$c" \
 --native-chain "$p-chain-$d.json" --native "${native[@]}" \
 --g8-reports "${g8r[@]}" --g8-states "${g8s[@]}" --g9-reports "${g9r[@]}" --g9-states "${g9s[@]}" \
 --mixed-reports "${m3r[@]}" --mixed-states "${m3s[@]}"
python -W error early_pole_endpoint_warm_evolved_single_cell_closure.py --control "$c" \
 --accepted "$a" --states "${g8s[3]}" "${g9s[3]}" "${m3s[3]}"
python -W error early_pole_endpoint_warm_checkpoint_clock_comparison.py --control "$c" \
 --accepted "$a" --g9 "${g9s[2]}" "${g9r[3]}" "${g9s[3]}" \
 "$w-gauss9-x04-1024-checkpoint-$d.json.gz" "$w-gauss9-x04-provenance-$d.sha256" \
 --mixed "${m3s[2]}" "${m3r[3]}" "${m3s[3]}" \
 "$w-mixed300-x04-1024-checkpoint-$d.json.gz" "$w-mixed300-x04-provenance-$d.sha256"
```

These are fixed same-start segment comparisons, not global time-error bounds, spatial convergence, or a full physical prediction.
