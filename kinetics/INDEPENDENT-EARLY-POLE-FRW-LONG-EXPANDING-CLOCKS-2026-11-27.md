# Matched expanding five-channel clock comparison at 7.2e−5

Run both expanding fine-archive trajectories through the same log-scale interval 7.2e−5 and 96 macro steps. Both dielectric responses are held fixed, both occupations use the native step reconstruction, and the modeled channels and energy ledger are the same. The only change is whether the represented dark energy enters H. The checkpointed solver [recovered](INDEPENDENT-EARLY-POLE-FRW-SAME-PATH-CLOCK-2026-11-27.md) a few near-tolerance elastic midpoint solves by bisection; it records this independently for each arm.

| Expanding clock | Integrated scalar gain | Final interval mean source | Macro steps with elastic bisection |
| --- | ---: | ---: | ---: |
| Radiation reference | 7.840596431e−6 | .117705872 | 86 |
| Conditional full H | 7.631411815e−6 | .114642390 | 80, 83, 91 |

Including the modeled dark energy in H reduces the accumulated gain by 2.668%. Both mean rates remain positive and increasing in the last interval; this is not an observation of freeze-out. The absolute trace-subtracted comoving energy errors are below 7.3e−13, and the entropy ledger residuals below 2.8e−12. This comparison isolates clock choice **at the same prescribed expansion**, not the effect of redshift relative to a fixed-scale solve, and it is not a temporal convergence study. Since the bisections differ and the same-path constant/shortened reference controls are separate, no claim is made about the much smaller dynamic clock correction.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1
python -W error early_pole_joined_screened_frw_reference.py fine --duration 7.2e-5 --steps 96 --expands
python -W error early_pole_joined_screened_frw_reference.py fine --duration 7.2e-5 --steps 96 --expands --full-H
```
The same macro-boundary recovery also works with a fresh `--checkpoint` path for each trajectory. Raw: `results/early-pole-joined-frw-long-expanding-clocks-2026-11-27.json.gz`.
