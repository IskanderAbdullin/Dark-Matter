# Late continuation from archived shells with no thermal restart

This source-free calculation uses the actual scalar and stable-fermion
entropy variables at the end of the 48-step 12 MeV run. It does not
reinitialize either species to a thermal shape. The comoving momentum
grid is fixed; both particle counts and total energy including pressure
work are tracked. Below-pole elastic scattering is included throughout.
The background uses either a freely redshifted photon-only comparison or
the [tabulated radiation clock](INDEPENDENT-LATE-RADIATION-AND-RESIDUAL-SOURCE-2026-11-26.md)
plus the explicitly retained dark density.

The on-shell parent is optional in the interval where the numerical
support permits it. For both daughter momenta below p_max, no parent
with mass M is kinematically possible once p_max drops below the
rest-frame decay momentum P. This follows from the maximal invariant
energy at equal upper momentum and opposite directions. The numerical
grid inserts the endpoint of the retained on-shell interval in log(b),
and checks that the parent has no sampled contributions beyond it.
For the archived upper grid edge the geometric end is b=3.63867.
This is not a bound on reactions involving particles outside the mesh.

The continuum common state factors the stiff prompt affinity. If the
implicit solve fails, retry from the other controlled initial guess or
subdivide the fixed-background collision update. The published residual
criterion is unchanged. For momentum refinement, the original linear
kinetic entropy variables are restricted exactly to subdivided cells
before any evolution, so no extra thermalization is assumed.

## Result

At b=410 with tabulated radiation and retained dark gravity the nominal
24-step run gives scalar and fermion effective temperatures of 13.90
and 0.403 keV, compared with 11.77 and 0.470 keV in the photon-only
comparator. The scalar count divided by the excited-state capacity at
its moment-equivalent temperature is 0.00221, versus 0.00363 for that
photon-only comparator. The coarse/fine results test numerical
sensitivity, not the accuracy of a late cosmological history.
Dropping the retained on-shell channel on the same radiation clock
changes the 24-step scalar ratio from 0.00221 to 0.00219.
At 36 steps the radiation-clock value is 0.00233, and tightening the
angular quadrature at 24 steps gives 0.00230. The highest resolved
scalar occupation at the final time is only about 0.0011 in the baseline
36-step case. Separate count and pressure/energy ledgers close within
6e-12 relative to their initial reference scales. No ground-state
degree of freedom is included.

The scalar count below the original q/T_R=0.005 edge in the radiative
run is about 4.1e-8 of all scalars, against 1.9e-8 for pure redshift
from the same archived shape. This retains strong memory of the input
shape; the earlier calculation reinitialized at b=17.7 gave about
2.0e-7 after that different preparation. The first-bin count is not a
singular occupation at zero momentum. These are finite-domain, source-free
late benchmarks, not an exclusion of condensation for other histories.

## Reproduce

Prepare the [external EOS inputs](INDEPENDENT-BENCHMARK-SETUP.md), then
from `kinetics` run:

```sh
export OPENBLAS_NUM_THREADS=1
BASE=(EOS.dat results/12-quantum-48-2026-11-26.json results/12-quantum-48-shape-2026-11-26.json)
python -W error late_archived_shape_continuation_reference.py "${BASE[@]}" --steps 24 --prompt-until 8 --radiation-table EOS.late.dat --dark-gravity
python -W error late_archived_shape_continuation_reference.py "${BASE[@]}" --steps 24 --prompt-until 8
python -W error late_archived_shape_continuation_reference.py "${BASE[@]}" --steps 12 --prompt-until 0 --off
```
