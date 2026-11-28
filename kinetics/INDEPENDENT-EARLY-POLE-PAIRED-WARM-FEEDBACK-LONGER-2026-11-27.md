# The joined warm-marker feedback remains tiny on a longer common grid

Continue the [shared-grid warm-marker comparison](INDEPENDENT-EARLY-POLE-JOINED-EXTERIOR-PAIRED-COLD-2026-11-27.md) from the same fine final n216 hard state, holding scale fixed, with the q=.1–5 marker selection. Increase the exposure from 2e−6 to 2.4e−5. Both joined resolutions use refreshed unlike screening and the same stricter 5e−15 per-step elastic acceptance, without changing the terminal conservation tests. The third row suppresses the hard channels and leaves reciprocal decay only.

| Run at 2.4e−5 | Baseline cold gain | Warm gain | Warm-to-cold feedback |
| --- | ---: | ---: | ---: |
| Joined, 72 macro steps | 2.845455205e−6 | 6.597486966e−12 | −1.181883874e−16 |
| Joined, 144 macro steps | 2.845542970e−6 | 6.597537357e−12 | −1.185262206e−16 |
| Decay only, 72 RK4 steps | 2.818961451e−6 | 6.595892078e−12 | −1.102650200e−16 |

The joined cold feedback is about −4.17e−11 of the baseline cold gain, and its end value changes by .286% between the two joined widths. The last two sampled values of feedback divided by exposure squared are near −2.05e−7; the earlier samples are less precise. The matched decay-only result has the same sign and order. Both joined runs close terminal energy to 8.9e−16 or better. These are finite fixed-scale histories, with no claim at later times or about abundance.

There are two additional checks on a separate short calculation. The [n6 reciprocal pole-tangent feedback](INDEPENDENT-EARLY-POLE-PAIRED-COLD-TANGENT-TIME-2026-11-27.md) is unchanged in both representations when imposing the same 5e−15 acceptance used at n12, so that budget difference did not cause the reported ordering reversal. In the manufactured n12 experiment multiplying reciprocal warm fluxes by 32, the cold differential per unit multiplier is −7.93e−19 for markers and −7.76e−19 for q4. Compared with the unscaled differential, the magnitudes increase by about 9% and 3%, respectively. This strengthens the sign and magnitude check but does not restore a precision ranking between shapes.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
a=results/early-pole-joined-fine-216-memory-2026-11-27.json.gz
p=early_pole_live_scalar_exterior_paired_difference_reference.py
for n in 72 144; do
  python -W error "$p" "$a" --duration 2.4e-5 --steps "$n" --order 8 --joined --elastic-energy-budget 5e-15
done
python -W error "$p" "$a" --duration 2.4e-5 --steps 72 --order 8
# For the separate short checks, use q4 with --warm-order 4 or omit it for markers.
python -W error "$p" "$a" --duration 2e-6 --steps 6 --order 8 --joined --anchor-pole --linearize-cold-pole --elastic-energy-budget 5e-15
python -W error "$p" "$a" --duration 2e-6 --steps 12 --order 8 --joined --anchor-pole --linearize-cold-pole --elastic-energy-budget 5e-15 --warm-multiplier 32
```

Raw: `results/early-pole-paired-warm-feedback-longer-2026-11-27.json`.
