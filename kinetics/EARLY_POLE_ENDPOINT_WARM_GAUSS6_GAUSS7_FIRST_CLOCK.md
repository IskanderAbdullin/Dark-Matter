# G7/G6 pressure contrast survives a matched first-continuation time replay

Recover the G7 512-step endpoint from its accepted fine x=.0075 startup.
The independently stored coarse time trace agrees exactly on the tracked
number, continuity, and time entries. Together with the existing G6
recovery, both orders use their own original startup and the same
native S4 mass for pressure evaluation at x=.01.

The G7/G6 relative pressure excess is 0.00059080241604 on 512 steps
and 0.00059080241466 on 1024 steps. Their difference is −1.38e−12,
as a relative fraction, against an order effect of about +0.05908%.
This excludes a clock effect of that magnitude **on the matched first
continuation only**. It gives no bound on startup time error, on
subsequent epochs, or on physical pressure. The independent comparison
uses the numeric clock fields of the earlier G6 replay; its original
legacy scope caption described the default final clock.

```sh
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
d=2026-11-28; p=results/early-pole-fine-charge2-scalar4
python -W error early_pole_endpoint_warm_g7_first_clock_check.py \
 --native "$p-x01-state-$d.json.gz" --prefix "results/early-pole-warm-gauss7-prefix-x01-$d.json" \
 --raw-reports "$p"-gauss{6,7}warm-x01-$d.json.gz \
 --starts "$p"-gauss{6,7}warm-start-state-$d.json.gz \
 --coarse "$p"-gauss{6,7}warm-coarse-x01-state-$d.json.gz \
 --fine "$p"-gauss{6,7}warm-x01-state-$d.json.gz \
 --coarse-replays "$p"-gauss{6,7}warm-first-pressure-clock-$d.json
```
