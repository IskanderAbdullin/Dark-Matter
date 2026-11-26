# Resolve the driven heavy lag inside two corrected fine q64 intervals

The [earlier heavy-lag check](INDEPENDENT-WARM-HEAVY-RELAXATION-LAG-2026-11-26.md)
interpolated the stationary heavy population from saved endpoints.
Here instead first interpolate the daughter entropy coefficients
in ln(a), at fixed normalized kinetic coordinate on the
redshifting momentum shells. Independently recalculate the
angular blocking, inverse source and instantaneous stationary
parent occupation at each interior sample. The selected visible
radiation follows the tabulated clock, the retained dark energy
comes from the interpolated phase, and the physical exterior
energy is linearly interpolated from its endpoints for gravity.
No exterior particles participate in the driven collision equation.

On the original fine q64 step 24 (near the native cold spike)
and step 48 (the terminal interval), solve the driven heavy equation
on fixed comoving parent momenta from 0 to 1 GeV in the local
step-start units. Approximate the independently computed equilibrium
linearly and use the lesser of the two sampled relaxation
coefficients in each small cell. This coefficient is a prescribed
approximation, not a lower enclosure between the samples.
Start at the instantaneous equilibrium for that cell sequence.

The update propagates the **lag** directly, avoiding a subtraction
of almost equal populations. At 16 subdivisions the least
integrated coefficient is at least 1.27e9 in either interval.
The maximum sampled absolute energy lag in units of retained
dark energy is 1.283e-20 (step 24) and 7.391e-21 (step 48).
The terminal quasi-static heavy energy fractions are 1.070e-7
and 6.469e-8. Refining the time interpolation from 8 to 16
moves the terminal lag by under 0.1%. Separate angular and radial
controls leave the displayed energy-lag conclusion unchanged.

These are driven comparisons on two smooth specified paths, not
bounds for a dynamically coupled heavy spectrum, unknown
outside daughter phases, or a short source burst between samples.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
d=ladder/wide48/q64
for i in 24 48; do
  python -W error warm_corrected_heavy_lag_entropy_surrogate_reference.py EOS.dat \
    "$d/step-$((i-1)).json" "$d/step-$i.json"
done
~~~

[Crossing interval](results/warm-q64-driven-heavy-lag-step24-2026-11-26.json.gz);
[terminal interval](results/warm-q64-driven-heavy-lag-step48-2026-11-26.json.gz).
