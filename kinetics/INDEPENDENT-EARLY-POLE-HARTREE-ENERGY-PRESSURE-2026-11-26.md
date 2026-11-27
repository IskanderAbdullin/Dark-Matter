# Tiny momentum does not make the self-mass-dominated shell pressureless

The [narrow-shell feedback](INDEPENDENT-EARLY-POLE-SHELL-HARTREE-FEEDBACK-2026-11-26.md)
changes the quasiparticle mass with scalar occupation.
Summing the resulting particle energies by itself
double counts part of the quartic energy.
The same simple Hartree closure has a conserving
variational completion. Set $A=m^2-m_0^2$, with
$m_0$ the vacuum mass and
$A=\lambda\int f/E\,d\Pi_p/2$, where
$d\Pi_p=d^3p/(2\pi)^3$ here is the number measure.
At fixed $f$,
$$
 \rho=\int E f\,d\Pi_p-\frac{A^2}{2\lambda},\qquad
 P=\int \frac{p^2}{3E}f\,d\Pi_p+\frac{A^2}{2\lambda}.
$$
Stationarity in $m^2$ is exactly the gap
equation, so changes of the pole cancel in
the first variation $\delta\rho/\delta f=E$.
Redshifting each momentum and diluting its number
at fixed comoving occupation then gives
$d\rho/d\ln a=-3(\rho+P)$.
An independently added source would enter the
energy equation weighted by $E$; this note
does not run a closed source history.

For a very cold distribution in momentum
($p\ll m$) whose pole is mostly generated
by the quartic occupation,
$m^3\simeq\lambda n/2$. The interaction subtraction
is then $mn/4$, giving
$\rho\simeq3mn/4,\ P\simeq mn/4$, even though
the velocities $p/E$ are small.
Thus the pressure-to-density ratio is about one
third during that conditional regime. In the
zero-momentum limit keeping $m_0$ finite, the
adiabatic derivative is
$$
 \frac{dP}{d\rho}
 =\frac{m^2-m_0^2}{3m^2-m_0^2}.
$$
This is a homogeneous adiabatic derivative,
**not** a certified rest-frame sound speed,
perturbation closure or linear transfer result.

The archived shell examples reproduce the
continuity derivative to better than $10^{-6}$
under subsequent redshifting at fixed particle
count, and the source variation is checked
independently by changing one shell's occupation.
At the ends of the stated shell comparisons
the pressure ratio is approximately .33333 and
the double-counting subtraction is close to 25%
of the quasiparticle energy sum. In the **frozen**
control it takes subsequent scale-factor growth
of about $2.0\times10^7$ before the illustrative
mass-squared correction equals the vacuum value.
That ratio is not a prediction of an actual epoch:
the scalar distribution, interactions and visible
degrees of freedom can evolve.

The conclusion is conditional on this repulsive
Hartree description with fixed count after
its diagnostic handoff. A non-Hartree coherent
field, other self energies, gravitational
perturbations and collisional evolution require
their own treatment. Low momentum alone is
insufficient to justify a pressureless fluid
while the adopted interaction mass dominates.

~~~sh
python -W error early_pole_hartree_energy_pressure_reference.py \
 results/early-pole-shell-hartree-feedback-2026-11-26.json
~~~

[Continuity, energy variation and pressure controls](results/early-pole-hartree-energy-pressure-2026-11-26.json).
