# Can a rounded critical core be distinguished from a subcritical fit?

Compare a [clipped critical halo](INDEPENDENT-HALO-BOSE-PRIMORDIAL-CAP-CORE-2026-11-26.md)
at finite $F$ with an ordinary smooth subcritical
Bose halo. Hold the two scalar central densities and
velocity dispersion $V_0$ equal. This compares how
two shapes look as static fits, **not** whether both
fit the same primordial phase-space constraints.

At large $F$, clipping below
$\delta=\ln(1+1/F)$ removes
$(4/3)\sqrt\delta+O(\delta^{3/2})$ from the
unnormalized central density integral.
An ordinary subcritical offset $\delta_s$
removes $\pi\sqrt{\delta_s}+O(\delta_s)$.
Matching the central densities gives
$\delta_s/\delta\to(4/3\pi)^2$,
where the notation means $[4/(3\pi)]^2$.
The subcritical fit therefore uses maximum
occupation approximately $5.55F$; it would
exceed the original cap if imposed on both.
The ratio of their central density gradients
with respect to potential tends to $3\pi^2/16$,
so matching density and dispersion does not
make the two phase-space cores identical.

Solve both Poisson equations with the same
dimensionless central density and (optionally)
the same static spherical Plummer baryon density
$\rho_b/\rho_{\phi,0}=b_0[1+(\xi/.3)^2]^{-5/2}$.
The baryonic term affects the potential as well
as the circular speed in both solutions.
For $F=7000$ and no baryons, the clipped model
has total circular speed about .0712% larger at
$\xi=.03$, .0307% at .1, .0101% at .3 and
.00246% at 1. In the $b_0=5$ case the
difference at .03 shrinks to .00683%.
At the tested $F=100$ without baryons the
difference at .3 is larger, about .591%.
At $F=7000$ a direct quadrature Poisson control
independently verifies the tabulated differences.

These differences use toy spherical components
at fixed matched parameters. They are not a
likelihood, angular-resolution requirement,
baryon constraint or statement about individual
galaxies. No cap is inferred from rotation
curves here. A finite halo and its envelope,
proper forward modeling of observations, and
the cosmological origin of any fitted DF remain
separate problems.

~~~sh
python -W error halo_bose_cap_shift_degeneracy_reference.py
~~~

[Matched density and circular-speed grids](results/halo-bose-cap-shift-degeneracy-2026-11-26.json).
