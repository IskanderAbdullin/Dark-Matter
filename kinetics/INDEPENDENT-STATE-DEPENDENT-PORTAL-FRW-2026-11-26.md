# Test of state-dependent visible portal injection during FRW evolution

This control replaces the continuing vacuum source by a signed source that
is evaluated at the current represented FD/BE state before each midpoint
source half step. It retains the prescribed visible entropy history and the
same collision operators as the previously benchmarked vacuum evolution.
The initial physical pulse is from an actually empty sector in both cases.

The direct light portal uses one positive visible weight for both forward
and inverse reactions and splits on both outgoing shell boundaries,
including the exterior. The heavy portal uses the locally adiabatic parent
occupation and the full lab angular response proportional to A+B. The
resulting positive and negative rates are kept separately up to the signed
moment update. A fixed-state check takes both channels from an
overpopulated equal-temperature bath and verifies that the output records
removal rather than injection. The evolving calculation rejects
inadmissible shell moments; passing the displayed controls does not ensure
that explicitly split source kicks are stable for arbitrary parameters.

The matched vacuum benchmark uses the same invariant and lab mesh orders
for portal production. Separate outgoing pair and daughter checks recover
the previously validated vacuum shell normalization. Dark feedback on
visible entropy, propagation of the eliminated heavy distribution, and
the exterior occupation history are not solved in this control. Treat the
quantum heavy source as a leading fixed-bath adiabatic response, not as
a resummation when the parent source itself dominates the decay scale.

A separate 50 MeV evolved snapshot audit (6 steps with the original
vacuum source) tests the published heavy kernel on entropy-variable
shells, not just the collisionless production proxy. At 45 MeV the
net heavy number correction is about 5.1e-10; the scalar injection
share is 0.479992. About 7.0e-5 of events send at least one daughter
outside the finite occupation grid. This fraction includes both low
and high exterior, without conditioning the parent branch on the grid.

## Matched 50 and 100 MeV controls

Both start empty, end at 0.9 T_R, and retain 0.05 <= q/T_R <= 16. The
12-interval results are:

| T_R | Final scalar energy fraction, vacuum | Final scalar energy fraction, state-dependent | Change in scalar collision heat / all retained dark input |
| ---: | ---: | ---: | ---: |
| 50 MeV | 0.299949463 | 0.299949756 | +1.85e-6 |
| 100 MeV | 0.325384621 | 0.325386966 | +3.92e-5 |

For stable fermion number the state-dependent correction is -1.98e-6
and -2.13e-5 relative to the respective matched vacuum run. Number
and energy ledgers close to 9.2e-14 or better. At both temperatures the
change in the final scalar share from 6 to 12 time steps is much larger
than the matched quantum correction. The comparison therefore checks
that these source effects are small for the tested bulk observables;
it does not resolve a continuum-limit correction at the shown digits.

## Reproduce

From `kinetics`, with `EOS.dat` the external central EOS and `empty.json`
the verified standalone vacuum production spectra:

```sh
export OPENBLAS_NUM_THREADS=1
python -W error joint_visible_quantum_source_reference.py
for TR in .05 .10; do
  for MODE in matched-vacuum-order quantum-visible-source; do
    python -W error joint_resonant_frw_split_reference.py EOS.dat --physical-amplitude --visible-source --"$MODE" --with-elastic --empty-pulse 1e-8 --empty-reheat "$TR" --support-low .05 --steps 6 12 --vacuum-scaled-residual
  done
done
python -W error joint_resonant_frw_split_reference.py EOS.dat --physical-amplitude --visible-source --with-elastic --empty-pulse 1e-8 --support-low .05 --steps 6 --vacuum-scaled-residual --save-snapshots snapshot.json
python -W error joint_visible_heavy_adiabatic_reference.py EOS.dat empty.json --reheat .05 --evolved-snapshots snapshot.json
```
