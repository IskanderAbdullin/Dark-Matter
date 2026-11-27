# Tiny cold feedback stays negative, but the warm-shape ordering is not resolved

Refine the [near-null reciprocal cold tangent](INDEPENDENT-EARLY-POLE-PAIRED-COLD-POLE-TANGENT-2026-11-27.md) on the same native hard snapshot. At fixed exposure 2e−6 both histories use a common baseline pole, exact anchored active gap and matched hard and scalar tangent deposits. The twelve-step replays also enforce the 5e−15 local elastic energy budget; the six-step records use ordinary acceptance.

| Joined macro steps | Native warm markers: cold feedback | Distributed q4: cold feedback |
| --- | ---: | ---: |
| 6 | −7.41416e−19 | −7.10356e−19 |
| 12 | −7.30187e−19 | −7.54698e−19 |

Both cases remain of order −3e−12 of their native cold gain, and the 12-step energy errors are below 7e−16. The *ranking by warm representation* reverses, however. The six-step difference cannot be read as a precise shape effect on that feedback. This does not undermine the much larger and separately measured direct warm-production difference. It also does not separate macro-width effects from the additionally tightened nonlinear acceptance.

An independent initial warm-direction check evaluates the cold source derivative under the added warm hard and scalar fluxes at the common snapshot. It yields leading quadratic coefficients for short-time cold feedback of about −1.94e−7 (markers) and −1.92e−7 (q4), similar in magnitude but not a precision ranking or a time-dependent integration. The hold-fixed comparison excludes expansion.

From `kinetics/`, with `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1`:
```sh
a=results/early-pole-joined-fine-216-memory-2026-11-27.json.gz
p=early_pole_live_scalar_exterior_paired_difference_reference.py
python -W error "$p" "$a" --joined --duration 2e-6 --steps 12 --order 8 --anchor-pole --linearize-cold-pole --elastic-energy-budget 5e-15
python -W error "$p" "$a" --joined --duration 2e-6 --steps 12 --order 8 --anchor-pole --linearize-cold-pole --elastic-energy-budget 5e-15 --warm-order 4
p=early_pole_live_scalar_exterior_decay_reference.py
python -W error "$p" "$a" --duration 2e-8 --panels 1 --linear-probe
python -W error "$p" "$a" --duration 2e-8 --panels 1 --linear-probe --warm-order 4
```

Raw: `results/early-pole-paired-cold-tangent-time-2026-11-27.json`.
