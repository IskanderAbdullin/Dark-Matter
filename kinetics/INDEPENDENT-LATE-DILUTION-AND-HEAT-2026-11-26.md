# Late dilution limits versus below-pole heat exchange

This independent check starts from the retained continuum-matched moments
and reconstructed shells at the end of the 50-to-12 MeV calculation.
No late visible source or change of either dark species count is assumed.
Use b = a/a_snapshot. This compares two limiting evolutions, rather
than assuming common kinetic equilibrium indefinitely.

1. **Collisionless.** Redshift the archived momentum shells exactly
   with p -> p/b and the physical density factor b^-3. Recompute the
   kinetic energy relative to that required at the continuum common
   normal-state saturation boundary at the same separate counts.
2. **Reversible common-temperature fluid.** Start from the inferred
   common entropy, keep entropy and both counts comoving-conserved,
   and solve the full massive FD/BE equations. The normal branch ends
   when the entropy equals that of the saturated normal state. Below
   that threshold, compute the conditional ideal ground population
   while allowing only the excited scalar and fermion entropy to carry
   the conserved entropy. The all-nonrelativistic analytic limit uses
   the saturated Bose zeta functions and the independently checked
   sparse, Boltzmann fermion asymptotic.

In the 48-step vacuum-source input, freely redshifted kinetic energy
remains at least 7.03 times the common-state threshold over the sampled
range b=1 to 1e9. By contrast, the hypothetical reversible common
branch reaches saturation at b=346.38 with dark temperature 2.108 keV.
The 24-step input gives b=346.73, and the 48-step quantum-source input gives b=346.38. The late ideal reversible continuation
approaches about 96% scalar ground population, but that branch requires
efficient heat exchange well past the onset. Its ground fraction is
not a dynamical prediction.

## Independent microscopic heat-rate check

On the reversible trial trajectory evaluate the exact invariant
below-pole scalar–stable-fermion scattering, including BE stimulation
and FD blocking. Compute both total scalar collision rate and
the temperature-difference eigenvalue at fixed species counts.
The latter weights every event by the squared lab energy transfer,
consistently with the derivative of the net two-temperature heat flow.
Its analytic equilibrium derivative agrees with a finite-difference
check at an unsaturated point. Near saturation the two quadrature
settings differ by less than 0.2%.

| b (48-step vacuum input) | Scalar scattering rate / photon-only H | Heat relaxation / photon-only H |
| ---: | ---: | ---: |
| 12.10 | 17.42 | 1.887 |
| 24.40 | 9.27 | 0.388 |
| 49.18 | 4.73 | 0.0694 |
| 99.16 | 2.40 | 0.01165 |
| 199.91 | 1.27 | 0.00185 |
| 346.38 (conditional saturation) | 1.22 | 0.000298 |

Here photon-only H is computed from T_gamma = T_visible,snapshot / b,
ignoring all other energy and any photon heating. It is a deliberately
small comparison expansion rate in a no-late-drain comparison, not a
modeled thermal history. Counting collisions would mask how little
energy they transfer. These below-pole rates are too small around the
conditional saturation point to justify imposing the reversible
common-temperature branch using this process alone. This does not
exclude other mechanisms, nonthermal transfer, or an actual condensate
once physical late evolution is modeled. No ground-state kinetic
equation is solved here.

## Reproduce

Prepare the external `EOS.dat` with the [shared benchmark setup](INDEPENDENT-BENCHMARK-SETUP.md).

From `kinetics`, using the external central EOS:

```sh
export OPENBLAS_NUM_THREADS=1
python -W error continuum_dark_dilution_limits_reference.py EOS.dat results/12-vacuum-48-2026-11-26.json results/12-vacuum-48-shape-2026-11-26.json > dilution.json
python -W error continuum_dark_normal_capacity_reference.py EOS.dat results/12-vacuum-48-2026-11-26.json --snapshots results/12-vacuum-48-shape-2026-11-26.json > capacity.json
python -W error late_scalar_fermion_heat_rate_reference.py capacity.json dilution.json
```
