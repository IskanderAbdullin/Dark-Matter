# Visible energy feedback and exterior source ledger in the coupled FRW test

This control replaces the imposed adiabatic visible thermal history after
each source kick. The visible equation of state still enforces local
thermal equilibrium. The expanding dark sector uses the tested signed
state-dependent portal and the same two nonlinear dark collision terms.
The endpoint is the **same scale factor** as the frozen-background control,
not a retuned common final temperature.

During redshift substeps the visible comoving entropy is conserved. At a
source half kick, subtract all net portal energy from the visible density
and invert the monotone equation of state to obtain the new temperature.
This subtracts energy carried by omitted source daughters as well as
energy carried by retained ones. The next half kick uses the updated
temperature; both kicks and collisions use the Hubble rate evaluated
before the midpoint kicks, by design. The generated energy in exterior
daughter bins is recorded separately. No flux across the momentum cut
from the restricted collision kernel is assumed.

The interpolated entropy data are unchanged. For a strict pressure audit,
integrate dP/dT = s over the same interpolant on each adjacent EOS knot
pair, retaining the pressure anchor. Between kicks independently integrate
the visible trace term a⁴(rho-3P) in log(a). A manufactured source with
both signs verifies the implementation against an independent temperature
ODE: the final temperature discrepancy drops from 1.50e-5 to 9.05e-7
under 6-to-24 step refinement. The pressure/energy ledger fits to 6e-15
of initial visible energy in that test. The source-free calculation
recovers comoving entropy conservation.

## Exterior energy

An optional collisionless tracer groups the positive net omitted source
by species, tail and injection time. It uses exactly the measured count
and injected energy. Redshifting a mean injected lab energy gives a
lower bound to the later energy because the redshift map is convex in
injected energy. The finite lower-tail support gives a chord upper
bound; on the high tail use the asymptotic slope and minimum support.
An independent distributed sample falls within these bounds. The traced
exterior contributes its lower-bound estimate to H, and its pressure
trace and energy are reconciled separately. This approximation does not
predict the exterior spectrum or model re-entry. A signed net removal
from an exterior bin is rejected by this positive-cohort control.

## Physical benchmark

| Reheating | Relative visible temperature change at fixed a | Energy withdrawn / initial visible energy | Final scalar energy share, frozen to coupled (12 intervals) |
| ---: | ---: | ---: | ---: |
| 50 MeV | -6.94e-6 | 2.98e-5 | 0.299949756 to 0.299949576 |
| 100 MeV | -5.08e-5 | 2.22e-4 | 0.325386966 to 0.325387264 |

The lower/upper final exterior energy estimates differ by less than
3e-6 of the exterior energy in these 12-interval cases, and including its
gravity in H changes the predicted visible temperature shift by at most
7e-13 in absolute relative-temperature units. The collective
visible/dark/exterior ledger closes to 4.2e-14 of initial visible energy, including separately integrated pressure trace terms.
These numerical checks do not validate the portal amplitude or prescribe
the treatment of out-of-domain states in a broader model.

## Reproduce

From `kinetics` with the same external central EOS:

```sh
export OPENBLAS_NUM_THREADS=1
python -W error visible_frw_energy_ledger_reference.py EOS.dat
python -W error frw_exterior_cohort_reference.py
for TR in .05 .10; do
  for MODE in accurate-visible-eos visible-backreaction; do
    python -W error joint_resonant_frw_split_reference.py EOS.dat --physical-amplitude --visible-source --quantum-visible-source --"$MODE" --with-elastic --empty-pulse 1e-8 --empty-reheat "$TR" --support-low .05 --steps 6 12 --vacuum-scaled-residual
  done
  python -W error joint_resonant_frw_split_reference.py EOS.dat --physical-amplitude --visible-source --quantum-visible-source --visible-backreaction --trace-exterior --with-elastic --empty-pulse 1e-8 --empty-reheat "$TR" --support-low .05 --steps 6 12 --vacuum-scaled-residual
done
```
