# Actual G9 and mixed300 through x=.05 with exact endpoint clocks

Using the exact accepted x=.04 restarts, both charts reach x=.05. At the common native mass, their pressure shifts relative to G8 are +0.02832649% (G9) and +0.02413877% (mixed300). Mixed300 captures 85.21623% of the G9 increment at this endpoint. The comparison changes warm nodal quadrature only, with fixed cold and charged axes.

The evolved absolute off-target pressure change is 2.576553e-08 of G8 total, or 1.067164e-04 of total absolute cellwise mixed300 pressure change. The corresponding number measures are 5.635411e-11 and 1.075747e-01. These quantities evaluate actual stored states through a one-cell replacement at measurement time, not a separately evolved replaced trajectory.

For the independent within-segment clocks 1024 and 2048, coarse-minus-fine relative count shifts are -3.180096e-05 and -3.180096e-05; pressure shifts are +2.524458e-07 and +2.489612e-07, respectively. The fractional mixed300/G9 contrast changes by -1.985079e-13 in count and +3.484454e-09 in pressure going from coarse to fine. A small change in the contrast does not bound each chart's absolute count error.

```bash
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
d=2026-11-28; p=results/early-pole-fine-charge2-scalar4; w=results/early-pole-warm
c="$w-single-g9-upgrade-$d.json"; a="$w-mixed300-actual-prefix-x05-$d.json"
native=("$p"-{start,x01,x02,x04,x05}-state-$d.json.gz)
g8r=("$p"-gauss8warm-{start,x01,x02,x04,x05}-$d.json.gz)
g8s=("$p"-gauss8warm-{start,x01,x02,x04,x05}-state-$d.json.gz)
g9r=("$p"-gauss9warm-{start,x01,x02,x04,x05}-$d.json.gz)
g9s=("$p"-gauss9warm-{start,x01,x02,x04,x05}-state-$d.json.gz)
m3r=("$p"-mixed300warm-{start,x01,x02,x04,x05}-$d.json.gz)
m3s=("$p"-mixed300warm-{start,x01,x02,x04,x05}-state-$d.json.gz)
python -W error early_pole_endpoint_warm_g7_prefix_check.py --levels 8 9 \
 --native-chain "$p-chain-$d.json" --native "${native[@]}" \
 --reports6 "${g8r[@]}" --states6 "${g8s[@]}" --reports7 "${g9r[@]}" --states7 "${g9s[@]}"
python -W error early_pole_endpoint_warm_300_held_clock_contrast.py --control "$c" \
 --native-chain "$p-chain-$d.json" --native "${native[@]}" \
 --g8-reports "${g8r[@]}" --g8-states "${g8s[@]}" --g9-reports "${g9r[@]}" --g9-states "${g9s[@]}" \
 --mixed-reports "${m3r[@]}" --mixed-states "${m3s[@]}"
python -W error early_pole_endpoint_warm_evolved_single_cell_closure.py --control "$c" \
 --accepted "$a" --states "${g8s[4]}" "${g9s[4]}" "${m3s[4]}"
python -W error early_pole_endpoint_warm_checkpoint_clock_comparison.py --control "$c" \
 --accepted "$a" --g9 "${g9s[3]}" "${g9r[4]}" "${g9s[4]}" \
 "$w-gauss9-x05-1024-checkpoint-$d.json.gz" "$w-gauss9-x05-provenance-$d.sha256" \
 --mixed "${m3s[3]}" "${m3r[4]}" "${m3s[4]}" \
 "$w-mixed300-x05-1024-checkpoint-$d.json.gz" "$w-mixed300-x05-provenance-$d.sha256"
```

These are same-start local segment results. They neither bound total time error back to startup nor establish spatial convergence or a physical abundance.
