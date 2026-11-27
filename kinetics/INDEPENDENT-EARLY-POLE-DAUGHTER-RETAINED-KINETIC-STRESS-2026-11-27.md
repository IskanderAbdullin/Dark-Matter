# Retained quasiparticle velocity is not the scalar fluid equation of state

Evaluate the accepted native/local terminal 32-mode histories
at both hard depths and the original-hard 78-node angular
histories. Keep the chart's scalar gap equation and its frozen
spectator occupations. In the retained window only, evaluate
the kinetic pressure weight (p^2/(3E)) at the final pole.
Separate the local-minus-native pressure difference into a
count change at the native mean weight, a shape correction
at fixed native pole, and a pole correction at fixed local
occupations. Form small differences with stable formulas and
independently check the subtracted gap equation.

| Hard depth × cold modes | Relative retained kinetic pressure of net production increment | Relative retained kinetic pressure of populated terminal state |
|---|---:|---:|
| 1 × 32 | +0.77336% | +6.632e-9 |
| 1 × 78 | +0.77220% | +6.620e-9 |
| 4 × 32 | +0.51547% | +4.435e-9 |

In the populated state the pole increase cancels about a third
of the count-only pressure term. It makes a much smaller relative
difference to the increment. Roughly 87–88% of the added
fixed-pole kinetic pressure is in old scalar groups five and six.
Pressure weighting and count weighting are not identical.

In the full scalar-only *bare-massless conformal Hartree model*
used by these charts, the pressure is (P_{m kin}+m^4/(2lambda)),
including the prescribed spectator component in (P_{m kin}).
The energy contains the opposite counterterm. This yields
(ho-3P=m^2sum n/E-2m^4/lambda=0) by the same scalar
gap equation. An independent finite dilation of the complete
state, re-solving its scalar pole, checks the pressure derivative.
The Hartree term supplies about 98.41% of that scalar-only
pressure in these fixtures. The paired full scalar pressure and
energy shifts are equal as relative changes, about 1.28e-8 on
the original hard mesh and 8.44e-9 on the finer hard mesh.
The represented scalar-only fluid has (P/ho=1/3) in
this conformal closure despite the much lower kinetic
velocities of particles in the selected narrow momentum
range. No physical bare mass, later transition, gravity or
dark-matter abundance follows from this calculation.

Run `early_pole_daughter_retained_kinetic_stress.py` with the
corresponding paired initializer and native/local raw histories.
Pass the common angular mesh via `--angular-mesh` for the
78-node history. The three receipts retain the independent
gap difference and the separate pressure terms.
