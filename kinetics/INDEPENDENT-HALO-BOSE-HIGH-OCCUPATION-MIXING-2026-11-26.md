# A primordial cap alone does not make a clipped critical core reachable

For collisionless evolution followed by coarse
phase-space averaging, the [excess measure](https://arxiv.org/abs/astro-ph/0504246)
$\mathcal D(f_*)=\int(f-f_*)_+\,d^3x\,d^3p$
cannot increase. This follows from invariance of
fine-grained occupation along Hamiltonian trajectories
and convexity of the positive part under averaging.
The integral can be taken in canonical momentum
coordinates with the same quantum phase-space
normalization at both times. Use only the inner halo
contribution as a **necessary** mixing test; an
unknown halo envelope can only raise its excess.

The [clipped illustrative core](INDEPENDENT-HALO-BOSE-PRIMORDIAL-CAP-CORE-2026-11-26.md)
attains its finite maximum $F$ on an open region of
six-dimensional phase space. For any distribution
bounded by $F$, provided the total particle count is
finite, $\mathcal D(F-\epsilon)/\epsilon$ tends to
the phase volume at exactly $F$ as $\epsilon\to0^+$.
The high-level support at $F-\epsilon$ has bounded
measure by finite particle count, so dominated
convergence applies. Thus no finite associated
Lagrangian patch whose initial distribution has
no plateau at $F$ can exactly produce the clipped
core at the **same** maximum by collisionless mixing.
A smaller final maximum or a different, rounded
energy distribution is not subject to that specific
obstruction. The initial capped waterbag is a
separate example with a nonzero plateau volume.

## A shifted, nonplateau core has a different tail

For a smooth subcritical Bose core
$f=[\exp(\delta+t+y)-1]^{-1}$ with
$\delta=\ln(1+F^{-1})$, its maximum is attained
only at zero position and momentum in this local
model. Write $\eta=1-f_*/F$. Near the center,
$1-f/F=(F+1)(t+y)+O[(t+y)^2]$, and both the
harmonic potential and nonrelativistic kinetic
energy are quadratic in their three coordinates.
Their six-dimensional phase volume at total
energy below $E$ scales as $E^3$;
$\mathcal D_\mathrm{core}$ consequently scales
as $\eta^4$.
The leading coefficient is checked against the
exact energy DF and the numerically solved potential.

For a spatially uniform finite initial patch, if
$f/F=1-Cu^n+o(u^n)$ at small radial primordial
momentum $u$, then
$\mathcal D_\mathrm{patch}\propto\eta^{1+3/n}$.
This follows by integrating
$u^2(\eta-Cu^n)$ up to $u=(\eta/C)^{1/n}$.
A primordial profile with $n<1$ thus also cannot
reach the **shifted** core at the same maximum
if the associated initial mass is finite. At $n=1$
the coefficients matter, and for $n>1$ this
leading endpoint argument alone does not rule it out.

## Numeric illustrations

The core count below is integrated to one local
radius $R$. Divide its excess per core particle
by the primordial excess per particle to obtain the
minimum associated patch count in core-count units
at a given threshold. For $F=7000$, use the
previous [compact logistic spectrum](INDEPENDENT-SCALAR-COLD-PACKING-NATIVE-CLASS-2026-11-26.md),
rescaled so its *actual peak* is $F$.
At threshold deficits $\eta=10^{-8}$ and
$10^{-10}$, the necessary patch/core count
ratios for the clipped example are roughly
$1.43\times10^{-5}$ and .0302. At $10^{-12}$
the ratio exceeds $10^4$.
This divergence at the endpoint is independent
of a successful CLASS result; it is also at such
an extreme threshold that it is not by itself an
observational test of an approximately clipped core.
The normalized shifted example has a much smaller
excess at the endpoint and no analogous plateau.

Neither passing any finite subset of thresholds
nor clipping a mathematical profile proves
collisionless halo formation. A halo with mergers,
mass loss or several components requires the correct
progenitor patch accounting; real particle collisions,
phase conversion and production can invalidate
occupation transport from the original epoch.

~~~sh
python -W error halo_bose_high_occupation_mixing_reference.py
~~~

[Per-threshold necessary core/patch comparisons](results/halo-bose-high-occupation-mixing-2026-11-26.json).
