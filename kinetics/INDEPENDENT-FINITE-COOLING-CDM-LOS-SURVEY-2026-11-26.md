# A split kinetic check and a CDM-metric tracer survey

This follows the [memory-bounded relativistic LOS calculation](INDEPENDENT-FINITE-COOLING-RELATIVISTIC-LOS-2026-11-26.md).
The 86 eV scalar has the same very small forward-only, first-order
prompt-plus-recoil test abundance and positive distribution.
The purpose is to test numerical grouping and the scalar tracer,
not to compute a self-gravitating scalar cosmology.

## Split hierarchy without a fluid closure

Use the two positive prompt and recoil inputs from the
[same-sum distribution check](INDEPENDENT-FINITE-COOLING-TRACER-SUPERPOSITION-2026-11-26.md).
Turn off the fluid closure and solve all six focused native modes
using the explicit RK evolver at lmax=36. The matched wavenumbers
are exactly k=4 and 10 h/Mpc. At z=20:

| k [h/Mpc] | Native joined, no fluid | Native split RK, no fluid, energy weighted | Difference (split − joined) |
| --- | ---: | ---: | ---: |
| 4 | 0.88129744 | 0.88129773 | +2.87×10⁻⁷ |
| 10 | 0.75966037 | 0.75966023 | −1.41×10⁻⁷ |

This recovers grouping consistency at the tested points without
reducing the momentum grids or multipole order to fit memory.
It is not an all-k hierarchy-convergence test. The individually
computed recoil ratios are positive, 0.59911 and 0.34893
at these modes. The difference from the independent fixed-metric
LOS control is of order 10⁻⁵ or smaller for the combined value
depending on whether the transfer or time-trajectory probe is
taken. This is decisively smaller than the earlier ~0.215
difference under the default one-fluid grouping at k=10,z=20.

## Survey using a separately generated CDM-only driver

Generate the metric in native CLASS with **no** scalar component, then
use the same metric in the relativistic collisionless LOS
integration for each scalar distribution. The scalar therefore does
not source the gravitational metric here. At the previous exact
k=4 and k=10 modes at z=20, this calculation agrees with the
matched native-recorded-metric LOS to within 1.1×10⁻⁷ absolute.
The joined continuous-input quadrature and the energy-weighted
separate responses agree within 7.1×10⁻⁹ at the selected probes.

For the joined scalar tracer, the fixed-metric values at k=22 h/Mpc
are 0.856585, 0.736954, and 0.604525 at z=0,5,20, respectively
(default wider *fluid* grouping: 0.725052, 0.470620, 0.260809).
None of the seven tested joined LOS samples from k=1 to 22 falls
below half at these redshifts. This deliberately does not certify
every intervening k or locate a continuous cutoff. The solver
keeps the approximate early adiabatic monopole, dipole and quadrupole,
checks its time/momentum refinement, and records its metric hashes.
A separate z=20 extension at k=26,30,34,40 gives 0.566204,
0.532253, 0.501896, and 0.461851 for the joined shape.
This brackets an amplitude half level between the last two
sampled modes in this conditional calculation, rather than at
the default-fluid numerical crossing near 10.7. There is no
continuous root solver or claim of an all-scale cutoff.
Total matter, physical abundances, and observations remain outside
this test.

## Reproduction

Run the focused split case generated with the usual baseline
inputs, no fluid closure and explicit RK:

~~~sh
OPENBLAS_NUM_THREADS=1 python -W error finite_cooling_focused_closure_reference.py native --prepare native
# Run the generated split_nfaoff_l36_rk_focus_sync.ini with native CLASS.
OPENBLAS_NUM_THREADS=1 python -W error finite_cooling_focused_closure_reference.py native --cases joined_on_l36 joined_off_l36 split_off_l36_rk
~~~

For the survey, generate a small CDM-only driver from the focused
joined-on input; its archive contains the exact non-scalar settings.

~~~sh
OPENBLAS_NUM_THREADS=1 python -W error finite_cooling_cdm_los_survey_reference.py \
  --seed native/finite_cooling_pairmix_tracer_ratio5_s14_nfaon_l36_focus_sync.ini --out cdm_native
# Run cdm_native/cdm_los.ini with native CLASS and retain its perturbation files.
for z in 20 5 0; do
  OPENBLAS_NUM_THREADS=1 python -W error finite_cooling_relativistic_los_reference.py \
    EOS.dat MIXED.dat PROMPT.dat RECOIL.dat \
    cdm_native/cdm_los_00_perturbations_k*_s.dat --cdm-only --z "$z"
done
~~~

The non-cold distributions are the same **unremapped** sigma=0.14
files used in the earlier relativistic LOS test. No new stiff
native massive-species evolution is needed for the broad driver.

For the separate high-k driver, give the same preparation command
a fresh output directory and add --modes 26 30 34 40. Run the
same LOS command on those outputs with --cdm-only --z 20.
