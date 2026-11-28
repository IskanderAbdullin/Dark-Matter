# Couple pair conversion to prompt and elastic evolution after the 12 MeV state

Restart from the [already paired, signed-source 50 MeV history](INDEPENDENT-50MEV-PAIR-FROM-EMPTY-2026-11-26.md)
at photon temperature 12 MeV. Split signed pair collisions around
the conservative implicit prompt-plus-elastic step. The two
conservative redshift half-steps and the independent exterior
cohort ledger record the expansion and outgoing tails. Use the
tabulated radiation with the included dark energy in the Hubble rate.
Prior unresolved cohorts are not returned as incident particles.

| Final relative expansion | Intervals | Included late pair conversion | Inclusive stable count / starting retained count |
| --- | ---: | --- | ---: |
| 1.2 | 4 | Yes | .95733 |
| 1.2 | 8 | Yes | .95641 |
| 1.2 | 8 | No | 1 (roundoff) |
| 1.5 | 12 | Yes | .89500 |

Thus the late loss is 4.36% by about 10.01 MeV and 10.5% by about
8.01 MeV, within this finite-support omitted-source comparator.
The scalar inclusive count at the latter point is 1.2623 times
its retained starting count. Its equivalent thermal saturation
ratio is about 3.87e-4, not near unity. The first-bin absolute
count is 2.037e-8 of the starting scalar number in the matched
8-step paired case, versus 1.980e-8 without conversion. A finite
first cell is neither a condensate nor a complete cold spectrum.

The signed reverse rate reaches about 84% of the forward rate
on this interval. Carrying the earlier external counts without
new reactions, the conditional fixed-count stable Ωh² at the
1.5 snapshot would still be about 78, not .120. Neither the
rates nor this intermediate count are an asymptotic relic
prediction. Visible residual injection, heavy recycling,
quartic multiplication, and re-entry are omitted. Only the
shorter endpoint has a step-refined partner here.

```sh
export OPENBLAS_NUM_THREADS=1
cmd=(python -W error coupled_warm_pair_continuation_reference.py EOS.late.dat results/50MeV-from-empty-pair-quantum-shape-2026-11-26.json)
"${cmd[@]}" --stop 1.2 --steps 4 8
"${cmd[@]}" --stop 1.2 --steps 4 8 --no-pair
"${cmd[@]}" --stop 1.5 --steps 12
```
