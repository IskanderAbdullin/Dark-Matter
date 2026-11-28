# Heavy visible portal and quantum response of prompt decay

This independent control computes the heavy-pair visible forward and reverse
rates after elimination of the rapidly decaying heavy parent. The stable
fermion and scalar distributions are prescribed for each snapshot.

For a parent of a fixed lab energy, use the exact two-body daughter angular
range with inverse factor A = f_phi f_chi and forward factor
B = (1+f_phi)(1-f_chi). Elimination in the absence of a visible perturbation
gives f_H = <A>/<A+B>. The forward and reverse visible heavy-pair rates use
that same parent occupation for each individual lab leg, with the thermal
visible detailed-balance factors already used for the direct light channel.
For a small external net parent source, the change in the dark decay and
inverse rates at fixed daughter distributions is δf_H (A+B). Thus the
conditional angular response is proportional to A+B, not simply B.
A separate forced two-angle ODE reproduces this response within 1e-12.

Compute the daughter count and kinetic-energy moments using disjoint exact
lab-angle gates for the scalar and fermion shells. The branch denominator
always integrates over the full daughter angle; outside a supplied
occupation grid the occupations are set to vacuum, not the phase-space
measure. Both exterior tails are retained in the projection. Number in
each daughter channel and summed daughter energy balance the parent
in each sign of the net heavy source. In vacuum the shell projection
agrees with the independent earlier analytic overlap integrator within
1.8e-15 normalized L1. Common-temperature thermal inputs with
eta_phi + eta_chi = 0 give a vanishing net visible heavy-pair rate
at floating-point precision.

## Fixed snapshots from 20, 50 and 100 MeV reheating

At the respective 0.9 T_R source-only snapshots, the largest sampled
parent occupations are about 4.7e-15, 2.0e-9 and 3.5e-6. The net heavy
number correction relative to vacuum reaches about 4.0e-7 in the
100 MeV control, where the inverse contribution is about 1.15e-11
of the vacuum number rate. The daughter angular weighting is more
noticeable. In that control, the scalar share of parent injection
energy goes from 0.480000 (vacuum isotropic) to 0.479003 and the
number fraction emitted below p/T=0.1 changes from 0.00468 to 0.00507.
The 50 MeV scalar share changes by about 1.87e-5 in absolute units.
The full-shell scalar-number L1 correction at the high-temperature
snapshot is 0.00359 after normalization to the vacuum scalar count.
The smallest decay response rate on the sampled parent geometry is
above 4.0e9 times visible H in all three controls. This is a local
response-time comparison, not a bound for an evolving parent distribution.

The softest portal tails from the parent quadrature are not precision
predictions; e.g., the relative numerical change in the very tiny
20 MeV snapshot tail can be large without affecting any bulk rates.
The exterior-vacuum convention and prompted, isotropic rest-frame
matrix element must not be taken as a self-consistent treatment of
all possible angular spin correlations or of an evolved out-of-domain bath.

## Reproduce

From `kinetics`, with the external central EOS and the raw source spectrum:

```sh
export OPENBLAS_NUM_THREADS=1
python -W error joint_visible_heavy_adiabatic_reference.py EOS.dat empty.json --shell-test
```
