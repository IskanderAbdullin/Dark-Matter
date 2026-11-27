# Flavor-mass shifts on the earlier unshifted fermion trajectories

The [frozen unequal mass check](INDEPENDENT-EARLY-POLE-FROZEN-FLAVOR-ASYMPTOTIC-MASS-2026-11-26.md)
shows that unequal occupancies displace the aperture at fixed x. A
shift in that aperture need not cause a comparable shift in the
time when the stimulated sign changes. Here we rerun the unshifted
separate-temperature hard pair/decay trajectory and evaluate alternative
instantaneous kinematics on its same fitted populations and scalar mass.
No shifted mass is included in the evolution, expansion, pair channels,
or heat exchange.

The earlier model's three baseline roots are reproduced. A common
photon contribution to the asymptotic mass is kept fixed at the
reference value to isolate the flavor-dependent fermion term.
Each FD occupation is transported as a function of the original
physical momentum: evaluate it at the bare energy associated with
the shifted pole momentum. The fixed-energy control instead evaluates
the FD fit at the shifted energy without adjusting the fit; it
is not the same preserved population.

At g=.005 and scalar q/T=4.1e−6:

| Imposed heat exchange / H_ref | Bare gate x | Shift of x, flavor masses and fixed momenta |
| ---: | ---: | ---: |
| 0 | .05019097817 | +5.23e−8 |
| 100 | .05031727257 | +8.62e−8 |
| 1500 | .05096196750 | −1.28e−10 |

The different imposed heat rates move the baseline by about 7.7e−4
across this set, much more than the post-processed kinematic shift.
This is a negative sensitivity result for the specific fitted
trajectories. It is not evidence that dynamical mass feedback or
an in-medium collision vertex can be neglected.

A separate first-order energy witness shows how to avoid counting
both quasiparticle shifts as independent energies if such feedback
is introduced. Let J_i be the hard moment integral of p f_i(p).
The massless-loop task admits a reciprocal quadratic energy density
g²/π⁴ [ J_gamma sum_i J_i + half sum_i J_i² ], where i runs over
the charged Dirac species and degeneracies are 2 for photons and
4 per Dirac species. Differentiating this functional per state
yields the desired lightlike photon and fermion shifts. Summing
all shifted quasiparticle energies at fixed occupations would count
this interaction term twice. An independent finite-difference check
verifies the first variation and the conformal trace for massless
external lines. A reciprocal illustrative finite-bare-energy
variant replaces J_i by the integral of p² f_i/E_bare,i in both
derivatives, not just one side; its pressure/trace relation matches
the scale derivative. It is not the matched finite-mass field theory.

Reproduce with early_pole_frozen_mass_gate_postprocessing_reference.py
(--exchange 0/100/1500) and early_pole_reciprocal_gauge_energy_functional_reference.py.
Raw files are in results/. No physical abundance or CLASS input follows.
