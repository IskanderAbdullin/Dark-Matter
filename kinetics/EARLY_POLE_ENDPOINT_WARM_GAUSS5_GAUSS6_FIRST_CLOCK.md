# The G5/G6 pressure contrast survives both first-continuation clocks

From each order's own archived x=.0075 fine startup state, independently recover its x=.01 coarse 512-step state. The exact input digest is checked against both the paired evolutionary report and its archived 1024-step endpoint. In both cases the four audited columns at all eight stored coarse trace times replay exactly. Pressure is evaluated throughout using the same native S4 mass at x=.01 and each chart's moment-conditioned reconstruction.

| Native S4 mass diagnostic | Value |
| --- | ---: |
| G5 coarse/fine pressure minus one | −5.87234e−8 |
| G6 coarse/fine pressure minus one | −5.87265e−8 |
| G6/G5 pressure minus one, both 512 steps | −0.000801522888414 |
| G6/G5 pressure minus one, both 1024 steps | −0.000801522885284 |
| Difference of those two order contrasts | +3.13e−12 |

The local clock sensitivity is much smaller than the warm order effect on this first continuation. Earlier startup states are left unchanged, so this result does not bound the accumulated pressure error from startup or subsequent stages.

```sh
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
d=2026-11-28; p=results/early-pole-fine-charge2-scalar4
for n in 5 6; do
 python -W error early_pole_endpoint_gauss_warm_pressure_clock.py --first-continuation \
  "$p-x01-state-$d.json.gz" "$p-gauss${n}warm-start-state-$d.json.gz" \
  "$p-gauss${n}warm-x01-state-$d.json.gz" "$p-gauss${n}warm-x01-$d.json.gz" \
  --coarse "/tmp/G$n-first-state.json" >"/tmp/G$n-first-clock.json"
done
python -W error early_pole_endpoint_warm_first_clock_order_contrast.py \
 --replay /tmp/G{5,6}-first-clock.json \
 --comparison "$p"-gauss{5,6}warm-comparison-x01-$d.json
```
