# Absolute G9 count and pressure at three fixed-start clocks

All three discretizations start at the exact same actual G9 x=.02 occupation and use the same geometry, capacity-limited step kernel and common mass when measuring x=.04. The coarser archives are replayed at their original clocks and the 4096-step arm is integrated directly. The compressed final 4096-step checkpoint binds the full recovered occupation, startup, time-step rule, geometry and ledgers; replaying it does not substitute a projected spectrum.

The adjacent fractional shifts (coarse divided by fine minus one) are -1.27670937e-04 and -3.17156958e-05 in scalar count and +3.01524206e-07 and +8.17095240e-08 in pressure. If the two signed logarithmic increments agree in sign, their magnitude ratio gives an observed refinement order, 2.0092 for count and 1.8837 for pressure. No asymptotic Richardson bound is imposed. Maximum absolute continuity and paired-number defects on the new arm are 4.343e-12 and 7.893e-17.

```bash
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
d=2026-11-28; p=results/early-pole-fine-charge2-scalar4-gauss9warm; w=results/early-pole-warm
tmp=$(mktemp); trap 'rm -f "$tmp"' EXIT
gzip -dc "$w-gauss9-x04-4096-checkpoint-$d.json.gz" > "$tmp"
python -W error early_pole_endpoint_warm_x04_g9_third_clock.py  "$p-x02-state-$d.json.gz" "$p-x04-state-$d.json.gz"  "$w-mixed300-actual-prefix-x04-$d.json" "$w-checkpoint-clock-comparison-x04-$d.json" "$tmp"
```

This quantifies one chart's absolute time sensitivity on .02 to .04. It does not bound error before the shared restart, prove any cross-chart cancellation at 4096 steps, or validate a physical production or transfer result.

## Origin of the clock change

The absolute cellwise number changes come almost entirely from the unrefined cold cells: more than 99.99999% at both pairs of clocks. The largest cell for 2048 versus 4096 is q in [3.306e-6, 3.625e-6], accounting for 27.06% of absolute number change. In contrast, refined warm cells account for 71.31% of absolute pressure change on that fine pair. The signed pressure contributions partially cancel: -5.49649e-8 from cold and +1.36674e-7 from warm, in units of the total 4096-step pressure. The same separation holds for 1024 versus 2048. This identifies which cells numerically change, not a spatial truncation error or an explanation of the full production history.

Reproduce the attribution after the previous check (same directory and variables):
```bash
python -W error early_pole_endpoint_warm_x04_g9_clock_cell_attribution.py \
 "$p-x04-state-$d.json.gz" "$w-gauss9-x04-1024-checkpoint-$d.json.gz" \
 "$w-gauss9-x04-4096-checkpoint-$d.json.gz" "$w-gauss9-x04-third-clock-$d.json"
```

The signed scalar-yield change is negative and increases in magnitude at all eight matched outputs, not just at the endpoint. For 1024 versus 2048 steps it grows from -4.8340e-6 at x=0.02181015 to -1.2767e-4 at x=.04. The three-clock observed count orders at these outputs range from 1.9737 to 2.0423; the endpoint order is 2.0092. This is an empirical result for the shared restart on this same chart.

```bash
python -W error early_pole_endpoint_g9_third_clock_count_history.py \
 "$p-x04-$d.json.gz" "$w-checkpoint-clock-comparison-x04-$d.json" \
 "$w-gauss9-x04-third-clock-$d.json"
```
