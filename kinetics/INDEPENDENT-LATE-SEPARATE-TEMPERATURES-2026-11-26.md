# Finite below-pole heat exchange with separate dark temperatures

This source-free late control projects each species separately to its own
FD or normal BE state. It is a thermal moment closure, not a kinetic
calculation of an evolving spectrum or ground mode. It uses the exact
represented separate counts and kinetic energies of the 48-step 12 MeV
snapshot as its input. The main run prepares a common thermal trial at
b=3 using the inferred initial common entropy. Independently test
separate no-exchange preparation at b=3 and common preparations at
b=1.7 and 5. This short preparation is explicitly prescribed, not an
evolution of the missed hot reactions.

The signed elastic heat uses the same event measure for forward and
reverse terms, independent scalar and fermion thermal integration
variables, a fixed below-pole cut s<0.82 M², and exact boosts. At equal
temperature the result is exactly zero. The direction of heat agrees
with detailed balance for both signs and a temperature ratio of 450.
The equilibrium derivative agrees with the earlier independent kernel
and finite differences. No prompt, number-changing or exterior return
processes are included. The background clock is only freely redshifted
photons, giving relatively long interaction time.

Redshift each species by conserving its own comoving entropy, with an
independent pressure-work quadrature. The collision step solves a signed
backward-Euler energy exchange at fixed counts and exactly opposite
energies; entropy cannot decrease. No exchange conserves the individual
entropies. Instantaneous equilibration with step refinement approaches
the separate reversible entropy calculation at b=140. It remains a
manufactured strong-coupling control, not the assumed physical result.

| Run at b=410 | Scalar T (keV) | Fermion T (keV) | Scalar count / excited-state saturation at its own T |
| --- | ---: | ---: | ---: |
| No heat, common initial preparation | 21.00 | 0.212 | 0.000642 |
| Physical heat, 12 time steps | 11.28 | 0.487 | 0.00412 |
| Physical heat, 24 time steps | 10.85 | 0.481 | 0.00463 |
| Physical heat, 48 time steps | 10.63 | 0.478 | 0.00492 |

Those counts are well below saturation despite the conditional shared
temperature adiabat crossing near b=346. The maximum ratio in the
48-step finite-rate run is 0.00493; the terminal ratios from the
alternate 24-step preparations range from 0.00424 to 0.00485.
Angular refinement at fixed time mesh changes this criterion by about
5e-9 relative, while time step effects are visible. The signed scalar
heat removes 0.0838 of initial total comoving kinetic energy in the
48-step control, separately from pressure work. Kinetic energy and
entropy ledgers close to better than 1e-10. This establishes a normal
state in the stated closure, not a no-condensate theorem for nonthermal
distributions or additional channels.

## Reproduce

```sh
cd kinetics
export OPENBLAS_NUM_THREADS=1
python -W error late_signed_scalar_heat_reference.py
BASE=(EOS.dat results/12-quantum-48-2026-11-26.json results/12-quantum-48-shape-2026-11-26.json)
python -W error late_two_temperature_expansion_reference.py "${BASE[@]}" --steps 48 --orders 6 6 8 4
python -W error late_two_temperature_expansion_reference.py "${BASE[@]}" --steps 24 --orders 9 9 12 6
python -W error late_two_temperature_expansion_reference.py "${BASE[@]}" --steps 18 --mode off
python -W error late_two_temperature_expansion_reference.py "${BASE[@]}" --steps 72 --mode equilibrate --stop 140
```
