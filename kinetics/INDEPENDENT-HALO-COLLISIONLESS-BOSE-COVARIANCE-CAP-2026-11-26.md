# A sharp local phase-space test without choosing a halo distribution

The [finite-cap critical-core control](INDEPENDENT-HALO-BOSE-PRIMORDIAL-CAP-CORE-2026-11-26.md)
chooses a specific energy-dependent halo distribution.
A weaker but more general test uses only the scalar
mass density, its local velocity covariance and a
finite fine-grained occupation cap $F$.
The same pointwise cap is necessary after collisionless
transport and coarse averaging; the full
[excess-mass constraint](https://arxiv.org/abs/astro-ph/0504246)
can still exclude profiles that pass a moment bound.

Let the nonrelativistic scalar velocities have mean
$\bar{\mathbf v}$ and positive-definite covariance
$\Sigma$. Transform to
$\mathbf u=\Sigma^{-1/2}(\mathbf v-\bar{\mathbf v})$,
so the mean squared radius is 3 and the momentum
Jacobian is $m^3\sqrt{\det\Sigma}$.
At fixed count and occupation cap, filling the smallest
$|\mathbf u|$ first minimizes the second moment.
This follows by comparing the signed difference from
the filled ball, subtracting the value of the cost
at its edge. A filled ball with mean square 3 has
radius $\sqrt5$; it also has exactly unit covariance.
The sharp local bound for one scalar state is
$$
 \rho_\phi\le
 \frac{5^{3/2}}{6\pi^2}
 \frac{F m^4\sqrt{\det\Sigma}}{\hbar^3}.
$$
Thus this necessary bound is valid for anisotropic or
nonthermal local velocity distributions. It is
saturated by an ideal filled ellipsoid, not generally
by a Bose thermal distribution or a formation history.

For $m=86$ eV and the illustrative scalar-only
density $0.1\,M_\odot\,{\rm pc}^{-3}$, isotropic
one-axis dispersion $20\,{\rm km\,s}^{-1}$ requires
$F\ge9.51$, while dispersions $20,20,10$ require
$F\ge19.03$. These are examples, not measured
scalar dispersions or densities. In particular a
gas or stellar rotation curve does not determine
the scalar covariance tensor.

A second check uses the moments of the ideal critical
noncondensed Bose core: $n=\zeta(3/2)/\lambda_{\rm dB}^3$
and $\sigma^2=V_0^2\zeta(5/2)/\zeta(3/2)$.
Those moments alone require $F\ge2.3875$, even
though the exact critical Bose DF has **no finite**
cap. The modeled late underoccupied thermal input
has peak occupation about .512. If collisionless
transport preserves that cap it cannot supply even
the ideal critical core's paired density and velocity
moments, independent of how its halo shape is fitted.
A high-occupation alternative may pass this moment
test without proving reachability. Self-interactions,
phase transfer or new production need a different
phase-space accounting.

~~~sh
python -W error halo_collisionless_bose_covariance_cap_reference.py
~~~

[Sharp-bound constants and fiducial density/dispersion ledger](results/halo-collisionless-bose-covariance-cap-2026-11-26.json).
